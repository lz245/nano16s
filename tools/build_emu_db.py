#!/usr/bin/env python3
"""Build an Emu database from NCBI 16S RefSeq Targeted Loci, optionally
deepened with rrnDB.

This is a portable rewrite of the cluster scripts in ONT-16S/docs/emu_db_update/.
It makes no assumptions about the machine: no wget, no GNU sed, no module
system, no scheduler. Standard library plus `emu` on PATH.

Why not just download a prebuilt database? The Emu project's default database
is built from NCBI 16S RefSeq as of September 2020 and has not been refreshed.
Species described since then are simply absent, so reads from them are either
unclassified or misassigned to a relative. This rebuilds from the current
release.

`--with-rrndb` adds rrnDB's 16S set to the build. This matters because
RefSeq Targeted Loci is a *type strain* collection: one authoritative sequence
per named species, which is broad and shallow by design. Measured on the
2026.08 build, the median species has exactly **one** reference sequence and
80% have only one. With so few references a nanopore read, which carries about
1% error, can sit closer to a sibling species' single reference than to its
own: on the ZymoBIOMICS mock this misnames four of eight species, every time,
with the genus always correct. rrnDB carries 16S genes pulled from complete
genome assemblies, many per species, which is the depth that fixes it.

Three details matter and are easy to get wrong:

1. NCBI's nodes.dmp uses the rank name "domain"; Emu expects "superkingdom".
   Without the substitution, every lineage terminates one level early.

2. NCBI assigns many 16S sequences a strain- or "no rank"-level taxid, but Emu
   requires species-level taxids in its seq2tax map. Each taxid is therefore
   walked up the tree to its species ancestor; sequences with no species
   ancestor are dropped, and the FASTA is filtered to match so the two files
   stay in agreement.

3. rrnDB lists every 16S copy in every genome it catalogues, so the same
   sequence recurs many times -- 265,316 records collapse to roughly a third
   of that. Identical sequences under one species taxid are therefore
   deduplicated; leaving them in inflates the database without adding any
   information, and makes minimap2 slower for nothing.

Usage:
    python tools/build_emu_db.py --outdir ~/.nano16s/db/2026.07
    python tools/build_emu_db.py --outdir DIR --with-rrndb
    python tools/build_emu_db.py --outdir DIR --cache DIR --keep-intermediates
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.error
import urllib.request
import zipfile
from datetime import date
from pathlib import Path

NCBI = "https://ftp.ncbi.nlm.nih.gov"
RRNDB = "https://rrndb.umms.med.umich.edu/downloads"

SOURCES = {
    "bacteria.16SrRNA.fna.gz": f"{NCBI}/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz",
    "archaea.16SrRNA.fna.gz": f"{NCBI}/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz",
    "bacteria.16SrRNA.gbff.gz": f"{NCBI}/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.gbff.gz",
    "archaea.16SrRNA.gbff.gz": f"{NCBI}/refseq/TargetedLoci/Archaea/archaea.16SrRNA.gbff.gz",
    "taxdump.tar.gz": f"{NCBI}/pub/taxonomy/taxdump.tar.gz",
}

# Only fetched with --with-rrndb. Version-pinned rather than "latest": a
# database whose contents depend on the day it was built cannot be compared
# with one built last month, and the manifest has to be able to say what went
# in. Bump these together when a new rrnDB release is adopted.
RRNDB_VERSION = "5.10"
RRNDB_SOURCES = {
    f"rrnDB-{RRNDB_VERSION}_16S_rRNA.fasta.zip":
        f"{RRNDB}/rrnDB-{RRNDB_VERSION}_16S_rRNA.fasta.zip",
    f"rrnDB-{RRNDB_VERSION}.tsv.zip":
        f"{RRNDB}/rrnDB-{RRNDB_VERSION}.tsv.zip",
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

# No timeout means block forever: urlopen defaults to the global socket
# timeout, which is None. An NCBI connection that stalls rather than closing
# left `nano16s db build` hanging with no output and nothing to wait for.
TIMEOUT_S = 120
ATTEMPTS = 3


def _download(url: str, tmp: Path) -> None:
    """One attempt. Verifies the length the server promised."""
    with urllib.request.urlopen(url, timeout=TIMEOUT_S) as r:
        declared = r.headers.get("Content-Length")
        with open(tmp, "wb") as out:
            shutil.copyfileobj(r, out)
    got = tmp.stat().st_size
    if declared is not None and got != int(declared):
        raise OSError(f"truncated: got {got} bytes, expected {int(declared)}")
    if got == 0:
        raise OSError("empty response")


def fetch(cache: Path, sources: dict[str, str] | None = None) -> None:
    """Download any source file not already cached."""
    cache.mkdir(parents=True, exist_ok=True)
    for name, url in (sources if sources is not None else SOURCES).items():
        dest = cache / name
        if dest.exists() and dest.stat().st_size > 0:
            log(f"  cached   {name} ({dest.stat().st_size / 1e6:.1f} MB)")
            continue
        log(f"  fetching {name} ...")
        tmp = dest.with_suffix(dest.suffix + ".part")
        # A connection cut mid-transfer ends cleanly as far as Python is
        # concerned, so an incomplete file used to be renamed into place and
        # then treated as cached forever -- every later build failed in
        # gzip.open with EOFError and nothing suggested clearing the cache.
        # The file is only renamed once its length matches.
        for attempt in range(1, ATTEMPTS + 1):
            try:
                _download(url, tmp)
                break
            except OSError as exc:      # URLError subclasses OSError
                tmp.unlink(missing_ok=True)
                if attempt == ATTEMPTS:
                    hint = ""
                    # rrnDB's server has served an incomplete certificate
                    # chain at times: the leaf names an intermediate it then
                    # does not send, so verification fails for every client
                    # that does not chase the AIA pointer (curl manages, most
                    # library HTTP clients do not). Nothing is wrong with the
                    # certificate itself, and the fix is not ours to make --
                    # but the user can supply the file by hand, and the cache
                    # is checked before anything is downloaded.
                    if isinstance(exc, urllib.error.URLError) and \
                            "CERTIFICATE_VERIFY_FAILED" in str(exc):
                        hint = (
                            f"\n\n"
                            f"       This is a certificate-verification failure, which for\n"
                            f"       rrnDB usually means their server is serving an incomplete\n"
                            f"       chain rather than anything being wrong here.\n\n"
                            f"       To work around it, download the file yourself and put it\n"
                            f"       in the cache, then re-run this command unchanged:\n\n"
                            f"           curl -L -o {cache / name} \\\n"
                            f"               {url}\n")
                    raise SystemExit(
                        f"ERROR: could not download {name} after {ATTEMPTS} "
                        f"attempts: {exc}\n"
                        f"       {url}\n"
                        f"       Nothing was cached, so re-running is safe.{hint}")
                log(f"           attempt {attempt} failed ({exc}); retrying")
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


def _zip_member(path: Path, suffix: str) -> str:
    """The single member of a one-file zip, checked rather than assumed."""
    with zipfile.ZipFile(path) as z:
        names = [n for n in z.namelist() if n.endswith(suffix)
                 and not n.startswith("__MACOSX/")]
    if len(names) != 1:
        raise SystemExit(
            f"ERROR: expected exactly one '*{suffix}' inside {path.name}, "
            f"found {len(names)}: {names[:5]}\n"
            f"       The download may be corrupt; delete it and re-run.")
    return names[0]


def parse_rrndb_taxids(cache: Path) -> dict[str, int]:
    """rrnDB genome record id -> NCBI taxid, from rrnDB-<v>.tsv.

    The FASTA headers carry the genome's record id in their second
    pipe-separated field, but no taxid; the genome table carries both. Joining
    the two is the only way to give rrnDB's sequences a taxonomy.
    """
    archive = cache / f"rrnDB-{RRNDB_VERSION}.tsv.zip"
    member = _zip_member(archive, ".tsv")
    mapping: dict[str, int] = {}
    with zipfile.ZipFile(archive) as z, z.open(member) as raw:
        text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
        for row in csv.DictReader(text, delimiter="\t"):
            rec = (row.get("Data source record id") or "").strip()
            tax = (row.get("NCBI tax id") or "").strip()
            if rec and tax.isdigit():
                mapping[rec] = int(tax)
    if not mapping:
        raise SystemExit(
            f"ERROR: no usable rows in {archive.name}. Expected columns "
            f"'Data source record id' and 'NCBI tax id'.")
    return mapping


def append_rrndb_fasta(cache: Path, dest: Path,
                       rec2taxid: dict[str, int]) -> dict[str, int]:
    """Append rrnDB's 16S sequences to `dest`, returning seqid -> taxid.

    rrnDB headers look like

        >Halobacterium salinarum NRC-1|GCF_000006805.1|NC_002607.1|...

    which has no accession in the first whitespace-delimited field -- the rest
    of this script keys sequences on exactly that field, and 'Halobacterium'
    is neither unique nor an identifier. Each record is therefore given a
    synthetic id, with the original header kept after it so the provenance of
    any sequence can still be traced in the finished database.
    """
    archive = cache / f"rrnDB-{RRNDB_VERSION}_16S_rRNA.fasta.zip"
    member = _zip_member(archive, ".fasta")
    mapping: dict[str, int] = {}
    n = skipped = 0
    with zipfile.ZipFile(archive) as z, z.open(member) as raw, \
            open(dest, "a", encoding="utf-8") as out:
        text = io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
        writing = False
        for line in text:
            if line.startswith(">"):
                fields = line[1:].rstrip("\n").split("|")
                rec = fields[1].strip() if len(fields) > 1 else ""
                taxid = rec2taxid.get(rec)
                writing = taxid is not None
                if not writing:
                    skipped += 1
                    continue
                n += 1
                seqid = f"rrnDB_{n:06d}"
                mapping[seqid] = taxid
                out.write(f">{seqid} {line[1:]}")
            elif writing:
                out.write(line)
    if skipped:
        log(f"  {skipped:,} rrnDB sequences had no genome-table entry and "
            f"were skipped")
    return mapping


def dedup_key(seq_lines: list[str]) -> bytes:
    """A hash of the sequence itself, case- and line-wrapping independent."""
    joined = "".join(seq_lines).replace("\n", "").replace("\r", "").upper()
    return hashlib.sha256(joined.encode("ascii", "ignore")).digest()


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
    ap.add_argument("--with-rrndb", action="store_true",
                    help="also include rrnDB's 16S set, which adds many "
                         "sequences per species. Fixes species-level "
                         "misnaming; costs roughly 4-5x the classification "
                         "time (see docs/USER_GUIDE.md)")
    args = ap.parse_args()

    if shutil.which("emu") is None:
        print("ERROR: 'emu' is not on PATH. Activate the environment first:\n"
              "         conda activate nano16s", file=sys.stderr)
        return 1

    outdir: Path = args.outdir.expanduser()
    work = outdir / "_build"
    work.mkdir(parents=True, exist_ok=True)

    cache = args.cache.expanduser()

    log("[1/6] Downloading NCBI sources")
    fetch(cache)
    if args.with_rrndb:
        log(f"       and rrnDB v{RRNDB_VERSION}")
        fetch(cache, RRNDB_SOURCES)

    log("[2/6] Merging FASTA")
    combined = work / "ncbi_16s_combined.fna"
    n_seqs = merge_fasta(cache, combined)
    log(f"  {n_seqs:,} sequences from NCBI RefSeq Targeted Loci")

    log("[3/6] Reading taxids from GenBank flat files")
    seq2taxid = parse_seq2taxid(cache)
    log(f"  {len(seq2taxid):,} accession -> taxid mappings")
    if len(seq2taxid) != n_seqs:
        log(f"  NOTE: {n_seqs - len(seq2taxid):,} sequences have no taxid and "
            f"will be dropped")

    n_rrndb = 0
    if args.with_rrndb:
        rec2taxid = parse_rrndb_taxids(cache)
        log(f"  {len(rec2taxid):,} rrnDB genome -> taxid mappings")
        rrndb_map = append_rrndb_fasta(cache, combined, rec2taxid)
        n_rrndb = len(rrndb_map)
        log(f"  {n_rrndb:,} sequences from rrnDB v{RRNDB_VERSION}")
        # rrnDB ids are synthetic and prefixed, so they cannot collide with a
        # RefSeq accession; assert it anyway rather than silently overwrite.
        clash = seq2taxid.keys() & rrndb_map.keys()
        if clash:
            print(f"ERROR: {len(clash)} ids appear in both sources, e.g. "
                  f"{sorted(clash)[:3]}", file=sys.stderr)
            return 1
        seq2taxid.update(rrndb_map)

    log("[4/6] Extracting taxonomy")
    taxonomy = work / "taxonomy"
    extract_taxonomy(cache, taxonomy)
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

    # Filter to the sequences that survived, deduplicating identical ones
    # under the same species. rrnDB lists every 16S copy of every genome, so
    # the same sequence recurs many times over -- without this the database
    # is roughly three times larger for no extra information, and minimap2
    # pays for all of it on every read. The whole record is buffered because
    # the sequence has to be hashed before it can be written.
    filtered = work / "ncbi_16s_filtered.fna"
    final: dict[str, int] = {}
    n_written = duplicates = 0
    seen: set[tuple[int, bytes]] = set()

    def flush(header, body, out):
        """Write one record if it is wanted and not a duplicate."""
        nonlocal n_written, duplicates
        if header is None:
            return
        # A header with nothing after ">" makes split()[0] an IndexError, ten
        # minutes into a build. It is not our file to validate, so skip the
        # record rather than abort.
        fields = header[1:].split()
        if not fields:
            return
        seqid = fields[0]
        taxid = kept.get(seqid)
        if taxid is None:
            return
        key = (taxid, dedup_key(body))
        if key in seen:
            duplicates += 1
            return
        seen.add(key)
        final[seqid] = taxid
        out.write(header)
        out.writelines(body)
        n_written += 1

    with open(combined, encoding="utf-8") as fin, \
         open(filtered, "w", encoding="utf-8") as out:
        header, body = None, []
        for line in fin:
            if line.startswith(">"):
                flush(header, body, out)
                header, body = line, []
            elif header is not None:
                body.append(line)
        flush(header, body, out)

    log(f"  filtered FASTA: {n_written:,} sequences")
    if duplicates:
        log(f"  {duplicates:,} identical sequences collapsed within a species")

    seq2tax_file = work / "seq2taxid.map"
    with open(seq2tax_file, "w", encoding="utf-8") as out:
        for acc, taxid in final.items():
            out.write(f"{acc}\t{taxid}\n")

    if n_written != len(final):
        print(f"ERROR: FASTA has {n_written} sequences but the map has "
              f"{len(final)} entries. Refusing to build a mismatched database.",
              file=sys.stderr)
        return 1

    log("[6/6] Running emu build-database (this shells out to minimap2)")
    dbname = "ncbi_16s"
    target = outdir / dbname

    # Build beside the existing database, not on top of it. `nano16s db build`
    # names the directory after the current month, so re-running in the same
    # month used to delete a working database and then rebuild it -- and an
    # Emu failure, a dropped connection or a Ctrl-C in the minutes that
    # followed left the user with neither. The old one is only replaced once
    # the new one is complete.
    staging = outdir / f"{dbname}.incoming"
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    proc = subprocess.run(
        ["emu", "build-database", dbname,
         "--sequences", str(filtered),
         "--seq2tax", str(seq2tax_file),
         "--ncbi-taxonomy", str(taxonomy)],
        cwd=staging,
    )
    if proc.returncode != 0:
        shutil.rmtree(staging, ignore_errors=True)
        print("ERROR: emu build-database failed. The existing database, if "
              "there was one, is untouched.", file=sys.stderr)
        return 1

    built = staging / dbname
    if not (built / "species_taxid.fasta").exists() \
            or not (built / "taxonomy.tsv").exists():
        shutil.rmtree(staging, ignore_errors=True)
        print(f"ERROR: expected outputs missing under {built}. The existing "
              f"database, if there was one, is untouched.", file=sys.stderr)
        return 1

    previous = outdir / f"{dbname}.previous"
    shutil.rmtree(previous, ignore_errors=True)
    if target.exists():
        target.rename(previous)
    built.rename(target)
    shutil.rmtree(staging, ignore_errors=True)
    shutil.rmtree(previous, ignore_errors=True)

    fasta = target / "species_taxid.fasta"
    taxtsv = target / "taxonomy.tsv"

    with open(fasta, encoding="utf-8") as fh:
        n_db = sum(1 for line in fh if line.startswith(">"))
    with open(taxtsv, encoding="utf-8") as fh:
        n_taxa = sum(1 for _ in fh) - 1

    manifest = {
        "built": date.today().isoformat(),
        "source": ("NCBI RefSeq Targeted Loci 16S (Bacteria + Archaea)"
                   + (f" + rrnDB v{RRNDB_VERSION}" if args.with_rrndb else "")),
        "source_urls": {**SOURCES,
                        **(RRNDB_SOURCES if args.with_rrndb else {})},
        "sequences": n_db,
        "taxa": n_taxa,
        "sequences_from_refseq": n_seqs,
        "sequences_from_rrndb": n_rrndb,
        "duplicates_collapsed": duplicates,
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
