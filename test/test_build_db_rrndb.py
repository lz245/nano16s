"""rrnDB as a second source for the reference database.

The point of including rrnDB is depth: RefSeq Targeted Loci carries about one
sequence per species, which is too few to separate close relatives, and rrnDB
carries many per species because it lists every 16S copy of every genome it
catalogues. That last property is also why these tests care so much about
deduplication -- the same sequence recurs constantly, and keeping every copy
would triple the database without adding information.
"""

import csv
import importlib.util
import io
import zipfile
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "tools" / "build_emu_db.py"

_spec = importlib.util.spec_from_file_location("build_emu_db", SCRIPT)
bdb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bdb)


# --- fixtures ------------------------------------------------------------

def write_genome_table(cache: Path, rows):
    """rrnDB-<v>.tsv.zip, with only the two columns the build reads."""
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(["Data source record id", "NCBI tax id", "NCBI scientific name"])
    for rec, tax, name in rows:
        w.writerow([rec, tax, name])
    path = cache / f"rrnDB-{bdb.RRNDB_VERSION}.tsv.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr(f"rrnDB-{bdb.RRNDB_VERSION}.tsv", buf.getvalue())
    return path


def write_rrndb_fasta(cache: Path, records):
    """records: (organism, genome_record_id, sequence)."""
    body = "".join(
        f">{org}|{rec}|NC_000000.1|Chromosome: X|1..2 +\n{seq}\n"
        for org, rec, seq in records
    )
    path = cache / f"rrnDB-{bdb.RRNDB_VERSION}_16S_rRNA.fasta.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr(f"rrnDB-{bdb.RRNDB_VERSION}_16S_rRNA.fasta", body)
    return path


# --- the genome table ----------------------------------------------------

def test_genome_table_maps_record_id_to_taxid(tmp_path):
    write_genome_table(tmp_path, [("GCF_1", "562", "Escherichia coli"),
                                  ("GCF_2", "1280", "Staphylococcus aureus")])
    assert bdb.parse_rrndb_taxids(tmp_path) == {"GCF_1": 562, "GCF_2": 1280}


def test_rows_without_a_numeric_taxid_are_skipped(tmp_path):
    write_genome_table(tmp_path, [("GCF_1", "562", "E. coli"),
                                  ("GCF_2", "", "no taxid"),
                                  ("GCF_3", "n/a", "not a number")])
    assert bdb.parse_rrndb_taxids(tmp_path) == {"GCF_1": 562}


def test_an_empty_genome_table_is_an_error_not_a_silent_empty_build(tmp_path):
    write_genome_table(tmp_path, [])
    with pytest.raises(SystemExit, match="no usable rows"):
        bdb.parse_rrndb_taxids(tmp_path)


# --- the FASTA -----------------------------------------------------------

def test_sequences_get_synthetic_ids_because_rrndb_headers_have_none(tmp_path):
    """The first whitespace field of an rrnDB header is a genus name, which is
    neither unique nor an identifier -- and the rest of the build keys on it."""
    write_rrndb_fasta(tmp_path, [("Escherichia coli K12", "GCF_1", "ACGT"),
                                 ("Escherichia coli O157", "GCF_2", "TTTT")])
    out = tmp_path / "combined.fna"
    out.write_text("")
    mapping = bdb.append_rrndb_fasta(tmp_path, out, {"GCF_1": 562, "GCF_2": 562})

    assert mapping == {"rrnDB_000001": 562, "rrnDB_000002": 562}
    headers = [ln for ln in out.read_text().splitlines() if ln.startswith(">")]
    assert headers[0].split()[0] == ">rrnDB_000001"
    # the original header is kept, so provenance survives into the database
    assert "GCF_1" in headers[0]


def test_sequences_whose_genome_is_not_in_the_table_are_dropped(tmp_path):
    write_rrndb_fasta(tmp_path, [("A", "GCF_1", "ACGT"), ("B", "GCF_missing", "TTTT")])
    out = tmp_path / "combined.fna"
    out.write_text("")
    mapping = bdb.append_rrndb_fasta(tmp_path, out, {"GCF_1": 562})
    assert list(mapping) == ["rrnDB_000001"]
    assert "GCF_missing" not in out.read_text()


def test_appends_rather_than_truncating_the_refseq_sequences(tmp_path):
    write_rrndb_fasta(tmp_path, [("A", "GCF_1", "ACGT")])
    out = tmp_path / "combined.fna"
    out.write_text(">NR_1.1 something\nGGGG\n")
    bdb.append_rrndb_fasta(tmp_path, out, {"GCF_1": 562})
    text = out.read_text()
    assert ">NR_1.1" in text and ">rrnDB_000001" in text


# --- deduplication -------------------------------------------------------

def test_identical_sequences_hash_the_same_whatever_the_wrapping(tmp_path):
    assert bdb.dedup_key(["ACGT\n", "ACGT\n"]) == bdb.dedup_key(["ACGTACGT\n"])


def test_case_is_ignored(tmp_path):
    assert bdb.dedup_key(["acgt\n"]) == bdb.dedup_key(["ACGT\n"])


def test_different_sequences_do_not_collide(tmp_path):
    assert bdb.dedup_key(["ACGT\n"]) != bdb.dedup_key(["ACGA\n"])


# --- the zip guard -------------------------------------------------------

def test_a_zip_with_two_candidate_members_is_refused(tmp_path):
    path = tmp_path / "two.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("a.tsv", "x")
        z.writestr("b.tsv", "y")
    with pytest.raises(SystemExit, match="exactly one"):
        bdb._zip_member(path, ".tsv")


def test_macos_resource_forks_are_not_mistaken_for_the_real_member(tmp_path):
    """A zip made on a Mac carries __MACOSX/ copies of every file."""
    path = tmp_path / "mac.zip"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("rrnDB.tsv", "x")
        z.writestr("__MACOSX/._rrnDB.tsv", "junk")
    assert bdb._zip_member(path, ".tsv") == "rrnDB.tsv"
