"""Per-read calls: one line per sequencing read, and only when asked for.

Emu does not label reads. It spreads each read over the references it matched
and estimates abundances from the whole distribution, so "what was this read?"
has no answer in the abundance tables. These cover the conversion of that
distribution into a per-read table, and the flag that turns it on -- which is
off by default because the distribution Emu writes is reads x taxa and can
reach gigabytes for a single deep barcode.
"""

import csv
import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"
CLI = ROOT / "bin" / "nano16s"
DEMO = ROOT / "test" / "demo" / "fastq_pass"


def make_db(tmp_path):
    """The two files the CLI checks for, plus a taxonomy to name taxids."""
    db = tmp_path / "db"
    db.mkdir()
    (db / "species_taxid.fasta").write_text(">1:x\nACGT\n")
    with (db / "taxonomy.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["tax_id", "species", "genus", "family", "order", "class",
                    "phylum", "clade", "superkingdom", "subspecies",
                    "species subgroup", "species group"])
        w.writerow(["562", "Escherichia coli", "Escherichia", "Enterobacteriaceae",
                    "Enterobacterales", "Gammaproteobacteria", "Pseudomonadota",
                    "", "Bacteria", "", "", ""])
        w.writerow(["1280", "Staphylococcus aureus", "Staphylococcus",
                    "Staphylococcaceae", "Bacillales", "Bacilli", "Bacillota",
                    "", "Bacteria", "", "", ""])
    return db


def make_dist(path, taxids, rows):
    """Emu's read-assignment distribution: one row per read, one column per taxid."""
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([""] + taxids)
        for read_id, probs in rows:
            w.writerow([read_id] + [repr(p) for p in probs])


def make_fastq(path, read_ids):
    with gzip.open(path, "wt") as fh:
        for r in read_ids:
            fh.write(f"@{r} runid=x\nACGT\n+\nIIII\n")


def convert(tmp_path, taxids, rows, read_ids, barcode="barcode01"):
    dist = tmp_path / "dist.tsv"
    make_dist(dist, taxids, rows)
    fastq = tmp_path / "reads.fastq.gz"
    make_fastq(fastq, read_ids)
    out = tmp_path / "per_read.tsv.gz"
    subprocess.run([sys.executable, str(SCRIPTS / "per_read_calls.py"), str(dist),
                    str(make_db(tmp_path)), str(fastq), barcode, str(out)],
                   check=True, capture_output=True, text=True)
    with gzip.open(out, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class TestConversion:
    def test_reports_the_taxon_holding_most_of_the_read(self, tmp_path):
        rows = convert(tmp_path, ["562", "1280"],
                       [("read1", [0.9, 0.1])], ["read1"])
        assert len(rows) == 1
        r = rows[0]
        assert r["read_id"] == "read1"
        assert r["barcode"] == "barcode01"
        assert r["taxid"] == "562"
        assert r["species"] == "Escherichia coli"
        assert float(r["confidence"]) == pytest.approx(0.9)
        assert r["status"] == "C"

    def test_full_lineage_is_written(self, tmp_path):
        rows = convert(tmp_path, ["562"], [("read1", [1.0])], ["read1"])
        assert rows[0]["lineage"] == ("Bacteria|Pseudomonadota|Gammaproteobacteria|"
                                      "Enterobacterales|Enterobacteriaceae|"
                                      "Escherichia|Escherichia coli")

    def test_a_split_read_is_marked_ambiguous(self, tmp_path):
        """Emu is probabilistic: a read fitting two references equally is
        normal, and rounding it to one species name without saying so is how a
        coin toss becomes a result."""
        rows = convert(tmp_path, ["562", "1280"],
                       [("read1", [0.45, 0.55])], ["read1"])
        assert rows[0]["status"] == "A"
        assert float(rows[0]["confidence"]) == pytest.approx(0.55)
        assert int(rows[0]["candidates"]) == 2

    def test_a_clear_call_is_confident(self, tmp_path):
        rows = convert(tmp_path, ["562", "1280"],
                       [("read1", [0.95, 0.05])], ["read1"])
        assert rows[0]["status"] == "C"

    def test_vanishing_probabilities_are_not_candidates(self, tmp_path):
        """Emu leaves the tail of the estimate in the row; 1e-20 is not a match."""
        rows = convert(tmp_path, ["562", "1280"],
                       [("read1", [1.0, 1e-20])], ["read1"])
        assert int(rows[0]["candidates"]) == 1

    def test_every_read_appears_even_if_it_matched_nothing(self, tmp_path):
        """Emu leaves unmatched reads out of the distribution entirely, so a
        file built from it alone would cover only part of the barcode."""
        rows = convert(tmp_path, ["562"], [("read1", [1.0])],
                       ["read1", "read2", "read3"])
        assert len(rows) == 3
        by_id = {r["read_id"]: r for r in rows}
        assert by_id["read1"]["status"] == "C"
        assert by_id["read2"]["status"] == "U"
        assert by_id["read2"]["species"] == ""
        assert by_id["read3"]["status"] == "U"

    def test_missing_distribution_says_what_to_do(self, tmp_path):
        out = tmp_path / "out.tsv.gz"
        fastq = tmp_path / "reads.fastq.gz"
        make_fastq(fastq, ["read1"])
        r = subprocess.run([sys.executable, str(SCRIPTS / "per_read_calls.py"),
                            str(tmp_path / "nope.tsv"), str(make_db(tmp_path)),
                            str(fastq), "barcode01", str(out)],
                           capture_output=True, text=True)
        assert r.returncode != 0
        assert "--per-read" in (r.stdout + r.stderr)


# These drive the real CLI, which builds a DAG and so needs Snakemake. The
# unit job installs only what the tests themselves import, so they are skipped
# there and run wherever the pipeline is actually installed -- locally, and in
# the install job.
@pytest.mark.skipif(not CLI.exists() or not DEMO.exists()
                    or shutil.which("snakemake") is None,
                    reason="needs the CLI, the demo data and Snakemake")
class TestTheFlag:
    """Off unless asked for -- including when asked for as the string "false"."""

    def plan(self, tmp_path, *extra):
        """Rule names Snakemake would run, from a dry run.

        Read from the job table, not from the whole output: Snakemake echoes
        the config, so "per_read=false" contains the rule's own name.
        """
        r = subprocess.run([str(CLI), "-d", str(DEMO), "-o", str(tmp_path / "out"),
                            "--db", str(make_db(tmp_path)), "-n", *extra],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stdout + r.stderr
        return {line.split()[0] for line in (r.stdout + r.stderr).splitlines()
                if line[:1].isalpha() and len(line.split()) == 2
                and line.split()[1].isdigit()}

    def test_absent_by_default(self, tmp_path):
        """REGRESSION: the CLI passes per_read=false to Snakemake, where it
        arrives as the *string* "false" -- and every non-empty string is true
        in Python. Read literally, the option switched itself on, and this is
        the case that caught it."""
        assert "per_read" not in self.plan(tmp_path)

    def test_present_when_asked_for(self, tmp_path):
        assert "per_read" in self.plan(tmp_path, "--per-read")
