"""A path with a space in it is refused before the run starts.

Emu builds its minimap2 command as a single string and splits it on whitespace,
and Porechop_ABI does the same for its adapter search. A path containing a
space therefore reaches both as two arguments. Left to run, the pipeline merges
every barcode, runs NanoStat over all of them, trims them, and only then fails
with `failed to open file 'probe/out'` -- half a path, naming a program the
user never invoked.

So these check the refusal, its wording, and that nothing is written first.

Run with:
    python -m pytest test/ -v
"""

import csv
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "bin" / "nano16s"
DEMO = ROOT / "test" / "demo" / "fastq_pass"


def make_input(parent: Path, name: str = "fastq_pass"):
    """A minimal run: one barcode holding one .fastq.gz."""
    bc = parent / name / "barcode01"
    bc.mkdir(parents=True)
    (bc / "reads.fastq.gz").write_bytes(b"")
    return parent / name


def make_db(parent: Path, name: str = "db"):
    db = parent / name
    db.mkdir(parents=True)
    (db / "species_taxid.fasta").write_text(">1:x\nACGT\n")
    with (db / "taxonomy.tsv").open("w", newline="") as fh:
        csv.writer(fh, delimiter="\t").writerow(["tax_id", "species", "genus"])
    return db


def run(*args):
    return subprocess.run([str(CLI), *args], capture_output=True, text=True)


@pytest.mark.skipif(not CLI.exists(), reason="CLI not present")
class TestSingleRun:
    def test_input_with_a_space_is_refused(self, tmp_path):
        d = make_input(tmp_path / "a dir")
        r = run("-d", str(d), "-o", str(tmp_path / "out"), "-y")
        assert r.returncode != 0
        assert "contains a space" in r.stderr
        assert "-d/--input-dir" in r.stderr

    def test_output_with_a_space_is_refused(self, tmp_path):
        d = make_input(tmp_path)
        r = run("-d", str(d), "-o", str(tmp_path / "out dir"), "-y")
        assert r.returncode != 0
        assert "contains a space" in r.stderr
        assert "-o/--output-dir" in r.stderr

    def test_database_with_a_space_is_refused(self, tmp_path):
        d = make_input(tmp_path)
        db = make_db(tmp_path / "db dir")
        r = run("-d", str(d), "-o", str(tmp_path / "out"), "--db", str(db), "-y")
        assert r.returncode != 0
        assert "contains a space" in r.stderr
        assert "--db" in r.stderr

    def test_a_tab_counts_too(self, tmp_path):
        d = make_input(tmp_path / "a\tdir")
        r = run("-d", str(d), "-o", str(tmp_path / "out"), "-y")
        assert r.returncode != 0
        assert "contains a space" in r.stderr

    def test_nothing_is_written_before_the_refusal(self, tmp_path):
        """The point of checking early: the old failure came after merge,
        NanoStat and Porechop had already run."""
        d = make_input(tmp_path)
        out = tmp_path / "out dir"
        r = run("-d", str(d), "-o", str(out), "-y")
        assert r.returncode != 0
        assert not out.exists()

    def test_the_message_offers_a_symlink(self, tmp_path):
        d = make_input(tmp_path / "a dir")
        r = run("-d", str(d), "-o", str(tmp_path / "out"), "-y")
        assert "ln -s" in r.stderr


@pytest.mark.skipif(not CLI.exists(), reason="CLI not present")
class TestBatch:
    def test_input_with_a_space_is_refused(self, tmp_path):
        parent = tmp_path / "runs dir"
        make_input(parent / "run1")
        r = run("batch", "-d", str(parent), "-o", str(tmp_path / "out"))
        assert r.returncode != 0
        assert "contains a space" in r.stderr

    def test_output_with_a_space_is_refused(self, tmp_path):
        parent = tmp_path / "runs"
        make_input(parent / "run1")
        r = run("batch", "-d", str(parent), "-o", str(tmp_path / "out dir"))
        assert r.returncode != 0
        assert "contains a space" in r.stderr


# Snakemake is absent from the unit CI job, which installs nothing but pytest.
@pytest.mark.skipif(not CLI.exists() or not DEMO.exists()
                    or shutil.which("snakemake") is None,
                    reason="needs the CLI, the demo data and Snakemake")
class TestControl:
    def test_a_path_without_spaces_still_runs(self, tmp_path):
        """The check must not refuse ordinary paths."""
        r = subprocess.run(
            [str(CLI), "-d", str(DEMO), "-o", str(tmp_path / "out"),
             "--db", str(make_db(tmp_path)), "-n"],
            capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        assert "contains a space" not in r.stderr
