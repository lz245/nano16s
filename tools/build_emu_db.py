#!/usr/bin/env python3
"""Build an Emu database from the current NCBI 16S RefSeq Targeted Loci set.

This is a portable rewrite of the cluster scripts in ONT-16S/docs/emu_db_update/.
It makes no assumptions about the machine: no wget, no GNU sed, no module
system, no scheduler. Standard library plus `emu` on PATH.

Why not just download a prebuilt database? The Emu project's default database
is built from NCBI 16S RefSeq as of September 2020 and has not been refreshed.
Species described since then are simply absent, so reads from them are either
unclassified or misassigned to a relative. This rebuilds from the current
release.

Two details matter and are easy to get wrong:

1. NCBI's nodes.dmp uses the rank name "domain"; Emu expects "superkingdom".
   Without the substitution, every lineage terminates one level early.

2. NCBI assigns many 16S sequences a strain- or "no rank"-level taxid, but Emu
   requires species-level taxids in its seq2tax map. Each taxid is therefore
   walked up the tree to its species ancestor; sequences with no species
   ancestor are dropped, and the FASTA is filtered to match so the two files
   stay in agreement.

Usage:
    python tools/build_emu_db.py --outdir ~/.nano16s/db/2026.07
    python tools/build_emu_db.py --outdir DIR --cache DIR --keep-intermediates
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from datetime import date
from pathlib import Path

NCBI = "https://ftp.ncbi.nlm.nih.gov"

SOURCES = {
    "bacteria.16SrRNA.fna.gz": f"{NCBI}/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz",
    "archaea.16SrRNA.fna.gz": f"{NCBI}/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz",
    "bacteria.16SrRNA.gbff.gz": f"{NCBI}/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz",
    "archaea.16SrRNA.gbff.gz": f"{NCBI}/refseq/TargetedLoci/Archaea/archaea.16SrRNA.gbff.gz",
    "taxdump.tar.gz": f"{NCBI}/pub/taxonomy/taxdump.tar.gz",
}

# Ranks at or above genus. Hitting one of these while walking up from a
# sequence's taxid means no species-level ancestor exists.
ABOVE_SPECIES = {
    "genus", "family", "order", "class", "phylum", "clade",
    "superkingdom", "domain", "kingdom",
}

_TAXON = re.compile(r'/db_xref="taxon:(\d+)"')


def log(msg: str) -> None:
    print(msg, flush=True)


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def fetch(cache: Path) -> None:
    """Download any source file not already cached."""
    cache.mkdir(parents=True, exist_ok=True)
    for name, url in SOURCES.items():
        dest = cache / name
        if dest.exists() and dest.stat().st_size > 0:
            log(f"  cached   {name} ({dest.stat().st_size / 1e6:.1f} MB)")
            continue
        log(f"  fetching {name} ...")
        tmp = dest.with_suffix(dest.suffix + ".part")
        with urllib.request.urlopen(url) as r, open(tmp, "wb") as out:
            shutil.copyfileobj(r, out)
        tmp.rename(dest)
        log(f"           {dest.stat().st_size / 1e6:.1f} MB")


# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------

def merge_fasta(cache: Path, dest: Path) -> int:
    """Concatenate the bacterial and archaeal FASTAs. Returns sequence count."""
    n = 0
    with open(dest, "w", encoding="utf-8") as out:
        for name in ("bacteria.16SrRNA.fna.gz", "archaea.16SrRNA.fna.gz"):
            with gzip.open(cache / name, "rt", encoding="utf-8") as fh:
                for line in fh:
                    if line.startswith(">"):
                        n += 1
                    out.write(line)
    return n


def parse_seq2taxid(cache: Path) -> dict[str, int]:
    """Map accession.version -> taxid by reading the GenBank flat files.

    The alternative is nucl_gb.accession2taxid.gz, which is ~9 GB and needs
    tens of GB of RAM to hold. The .gbff files carry the same information for
    exactly the sequences we care about, in ~20 MB.
    """
    mapping: dict[str, int] = {}
    for name in ("bacteria.16SrRNA.gbff.gz", "archaea.16SrRNA.gbff.gz"):
        accession = None
        with gzip.open(cache / name, "rt", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("VERSION"):
                    parts = line.split()
                    accession = parts[1] if len(parts) >= 2 else None
                elif accession and '/db_xref="taxon:' in line:
                    m = _TAXON.search(line)
                    if m:
                        mapping[accession] = int(m.group(1))
                        accession = None
    return mapping


def extract_taxonomy(cache: Path, dest: Path) -> None:
    """Unpack names.dmp and nodes.dmp, rewriting NCBI's 'domain' rank.

    Emu's taxonomy builder looks for 'superkingdom'. NCBI renamed that rank to
    'domain'. Left alone, every lineage stops at phylum.
    """
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(cache / "taxdump.tar.gz") as tar:
        for member in ("names.dmp", "nodes.dmp"):
            src = tar.extractfile(member)
            if src is None:
                raise RuntimeError(f"{member} missing from taxdump.tar.gz")
            with open(dest / member, "wb") as out:
                shutil.copyfileobj(src, out)

    nodes = dest / "nodes.dmp"
    text = nodes.read_text(encoding="utf-8")
    patched = text.replace("\tdomain\t", "\tsuperkingdom\t")
    if patched != text:
        nodes.write_text(patched, encoding="utf-8")
        log("  patched nodes.dmp: rank 'domain' -> 'superkingdom'")


def load_nodes(nodes_dmp: Path) -> dict[int, tuple[int, str]]:
    """taxid -> (parent_taxid, rank)."""
    nodes: dict[int, tuple[int, str]] = {}
    with open(nodes_dmp, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("|")
            if len(parts) < 3:
                continue
            nodes[int(parts[0].strip())] = (int(parts[1].strip()), parts[2].strip())
    return nodes


def species_ancestor(taxid: int, nodes: dict[int, tuple[int, str]]) -> int | None:
    """Walk up to the species-level ancestor, or None if there isn't one."""
    seen: set[int] = set()
    current = taxid
    while current in nodes and current not in seen:
        seen.add(current)
        parent, rank = nodes[current]
        if rank == "species":
            return current
        if rank in ABOVE_SPECIES:
            return None
        if current == parent:
            break
        current = parent
    return None


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--outdir", type=Path, required=True,
                    help="where the finished database goes")
    ap.add_argument("--cache", type=Path,
                    default=Path.home() / ".nano16s" / "cache",
                    help="download cache (default: ~/.nano16s/cache)")
    ap.add_argument("--keep-intermediates", action="store_true",
                    help="keep the working files instead of deleting them")
    args = ap.parse_args()

    if shutil.which("emu") is None:
        print("ERROR: 'emu' is not on PATH. Activate the environment first:\n"
              "         conda activate nano16s", file=sys.stderr)
        return 1

    outdir: Path = args.outdir.expanduser()
    work = outdir / "_build"
    work.mkdir(parents=True, exist_ok=True)

    log("[1/6] Downloading NCBI sources")
    fetch(args.cache.expanduser())

    log("[2/6] Merging FASTA")
    combined = work / "ncbi_16s_combined.fna"
    n_seqs = merge_fasta(args.cache.expanduser(), combined)
    log(f"  {n_seqs:,} sequences")

    log("[3/6] Reading taxids from GenBank flat files")
    seq2taxid = parse_seq2taxid(args.cache.expanduser())
    log(f"  {len(seq2taxid):,} accession -> taxid mappings")
    if len(seq2taxid) != n_seqs:
        log(f"  NOTE: {n_seqs - len(seq2taxid):,} sequences have no taxid and "
            f"will be dropped")

    log("[4/6] Extracting taxonomy")
    taxonomy = work / "taxonomy"
    extract_taxonomy(args.cache.expanduser(), taxonomy)
    nodes = load_nodes(taxonomy / "nodes.dmp")
    log(f"  {len(nodes):,} taxonomy nodes")

    log("[5/6] Remapping to species-level taxids")
    kept: dict[str, int] = {}
    remapped = dropped = 0
    for acc, taxid in seq2taxid.items():
        sp = species_ancestor(taxid, nodes)
        if sp is None:
            dropped += 1
            continue
        if sp != taxid:
            remapped += 1
        kept[acc] = sp
    log(f"  kept {len(kept):,}  remapped {remapped:,}  dropped {dropped:,}")

    seq2tax_file = work / "seq2taxid.map"
    with open(seq2tax_file, "w", encoding="utf-8") as out:
        for acc, taxid in kept.items():
            out.write(f"{acc}\t{taxid}\n")

    filtered = work / "ncbi_16s_filtered.fna"
    n_written = 0
    with open(combined, encoding="utf-8") as fin, \
         open(filtered, "w", encoding="utf-8") as out:
        writing = False
        for line in fin:
            if line.startswith(">"):
                writing = line[1:].split()[0] in kept
                if writing:
                    n_written += 1
            if writing:
                out.write(line)
    log(f"  filtered FASTA: {n_written:,} sequences")

    if n_written != len(kept):
        print(f"ERROR: FASTA has {n_written} sequences but the map has "
              f"{len(kept)} entries. Refusing to build a mismatched database.",
              file=sys.stderr)
        return 1

    log("[6/6] Running emu build-database (this shells out to minimap2)")
    dbname = "ncbi_16s"
    target = outdir / dbname
    if target.exists():
        shutil.rmtree(target)

    proc = subprocess.run(
        ["emu", "build-database", dbname,
         "--sequences", str(filtered),
         "--seq2tax", str(seq2tax_file),
         "--ncbi-taxonomy", str(taxonomy)],
        cwd=outdir,
    )
    if proc.returncode != 0:
        print("ERROR: emu build-database failed.", file=sys.stderr)
        return 1

    fasta = target / "species_taxid.fasta"
    taxtsv = target / "taxonomy.tsv"
    if not fasta.exists() or not taxtsv.exists():
        print(f"ERROR: expected outputs missing under {target}", file=sys.stderr)
        return 1

    with open(fasta, encoding="utf-8") as fh:
        n_db = sum(1 for line in fh if line.startswith(">"))
    with open(taxtsv, encoding="utf-8") as fh:
        n_taxa = sum(1 for _ in fh) - 1

    manifest = {
        "built": date.today().isoformat(),
        "source": "NCBI RefSeq Targeted Loci 16S (Bacteria + Archaea)",
        "source_urls": SOURCES,
        "sequences": n_db,
        "taxa": n_taxa,
        "dropped_no_species_ancestor": dropped,
        "remapped_to_species": remapped,
    }
    (target / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if not args.keep_intermediates:
        shutil.rmtree(work)

    log("")
    log(f"Database ready: {target}")
    log(f"  {n_db:,} sequences, {n_taxa:,} taxa")
    log("")
    log("Use it with:")
    log(f"  nano16s -d /path/to/fastq_pass -o results/ --db {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
