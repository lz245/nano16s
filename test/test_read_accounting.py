"""Read accounting, and the label on Emu's unnamed row.

Two questions a user should not have to answer with a calculator: where did my
reads go, and how many taxa did each barcode produce. The counts tables answer
neither -- they hold one row per taxon plus a row Emu leaves unnamed, so a
column sum includes reads with no visible label, and dropping unnamed rows
loses them silently.
"""

import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


accounting = load("read_accounting")
labeller = load("label_unclassified")

HEADER = ["tax_id", "abundance", "species", "genus", "family", "order", "class",
          "phylum", "clade", "superkingdom", "subspecies", "species subgroup",
          "species group", "estimated counts"]


def emu_table(path, rows):
    """An Emu rel-abundance file: (tax_id, species, genus, counts) tuples."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(HEADER)
        for tax_id, species, genus, counts in rows:
            row = [tax_id, "", species, genus] + [""] * 9 + [counts]
            w.writerow(row)


def summary_csv(path, barcodes):
    """preprocessing_summary.csv: {barcode: (raw, filtered)}."""
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode", "raw_reads", "raw_bases", "raw_median_length",
                    "raw_median_quality", "filtered_reads", "filtered_bases",
                    "filtered_median_length", "filtered_median_quality"])
        for bc, (raw, filt) in barcodes.items():
            w.writerow([bc, float(raw), 0, 0, 0, float(filt), 0, 0, 0])


def run_accounting(tmp_path, barcodes, tables, want_taxa=False):
    emu_dir = tmp_path / "06_emu_output"
    for bc, rows in tables.items():
        emu_table(emu_dir / bc / f"{bc}_rel-abundance.tsv", rows)
    summary = tmp_path / "preprocessing_summary.csv"
    summary_csv(summary, barcodes)
    out = tmp_path / "read_accounting.tsv"
    taxa = tmp_path / "per_barcode_taxa.tsv"
    subprocess.run([sys.executable, str(SCRIPTS / "read_accounting.py"),
                    str(emu_dir), str(summary), str(out), str(taxa)], check=True,
                   capture_output=True, text=True)
    with out.open() as fh:
        acct = {r["barcode"]: r for r in csv.DictReader(fh, delimiter="\t")}
    if not want_taxa:
        return acct
    with taxa.open() as fh:
        return acct, list(csv.DictReader(fh, delimiter="\t"))


class TestCounting:
    def test_splits_classified_from_unclassified(self, tmp_path):
        rows = run_accounting(
            tmp_path,
            {"barcode01": (1000, 900)},
            {"barcode01": [("562", "Escherichia coli", "Escherichia", 700),
                           ("1280", "Staphylococcus aureus", "Staphylococcus", 150),
                           ("unmapped", "", "", 30),
                           ("mapped_unclassified", "", "", 20)]},
        )
        r = rows["barcode01"]
        assert int(r["reads_classified"]) == 850
        assert int(r["reads_unclassified"]) == 50
        assert int(r["reads_to_classifier"]) == 900
        assert r["check"] == "ok"

    def test_counts_species_and_genera(self, tmp_path):
        rows = run_accounting(
            tmp_path,
            {"barcode01": (100, 100)},
            {"barcode01": [("1", "Escherichia coli", "Escherichia", 40),
                           ("2", "Escherichia fergusonii", "Escherichia", 30),
                           ("3", "Bacillus subtilis", "Bacillus", 30)]},
        )
        assert int(rows["barcode01"]["species_found"]) == 3
        assert int(rows["barcode01"]["genera_found"]) == 2

    def test_chain_balances(self, tmp_path):
        """raw = removed + filtered, and filtered = classified + unclassified."""
        rows = run_accounting(
            tmp_path,
            {"barcode01": (1000, 900)},
            {"barcode01": [("562", "Escherichia coli", "Escherichia", 880),
                           ("unmapped", "", "", 20)]},
        )
        r = rows["barcode01"]
        assert int(r["raw_reads"]) == int(r["removed_by_filter"]) + int(r["filtered_reads"])
        assert int(r["filtered_reads"]) == (int(r["reads_classified"])
                                            + int(r["reads_unclassified"])
                                            + int(r["subsampled_out"]))

    def test_subsampling_is_visible(self, tmp_path):
        """With --max-reads the classifier sees fewer reads; the chain still balances."""
        rows = run_accounting(
            tmp_path,
            {"barcode01": (5000, 4000)},
            {"barcode01": [("562", "Escherichia coli", "Escherichia", 990),
                           ("unmapped", "", "", 10)]},
        )
        r = rows["barcode01"]
        assert int(r["reads_to_classifier"]) == 1000
        assert int(r["subsampled_out"]) == 3000
        assert int(r["filtered_reads"]) == (int(r["reads_classified"])
                                            + int(r["reads_unclassified"])
                                            + int(r["subsampled_out"]))

    def test_empty_barcode_is_reported_as_zero(self, tmp_path):
        """A barcode with no reads is a fact about the run, not an omission."""
        emu_dir = tmp_path / "06_emu_output"
        emu_dir.mkdir(parents=True)
        summary = tmp_path / "preprocessing_summary.csv"
        summary_csv(summary, {"barcode01": (0, 0)})
        out = tmp_path / "read_accounting.tsv"
        subprocess.run([sys.executable, str(SCRIPTS / "read_accounting.py"),
                        str(emu_dir), str(summary), str(out)], check=True,
                       capture_output=True, text=True)
        with out.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [r["barcode"] for r in rows] == ["barcode01"]
        assert int(rows[0]["reads_classified"]) == 0


def combined_table(path, rank, rows):
    """A combined counts table: (name, barcode01, barcode02) tuples."""
    # Real combined tables name each rank once: for the genus table the rank
    # column IS "genus". A duplicated column would make csv.DictReader keep
    # only the last, which is a quirk of the test, not of the data.
    extra = [c for c in ("genus", "family") if c != rank]
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([rank, *extra, "barcode01", "barcode02"])
        for name, a, b in rows:
            filled = {"genus": "Escherichia" if name else "",
                      "family": "Enterobacteriaceae" if name else ""}
            w.writerow([name, *(filled[c] for c in extra), a, b])


class TestLabelling:
    def test_unnamed_row_gets_a_name(self, tmp_path):
        p = tmp_path / "emu-combined-species-counts.tsv"
        combined_table(p, "species", [("Escherichia coli", "90", "80"),
                                      ("", "10", "20")])
        subprocess.run([sys.executable, str(SCRIPTS / "label_unclassified.py"),
                        "species", str(p)], check=True, capture_output=True)
        with p.open() as fh:
            names = [r["species"] for r in csv.DictReader(fh, delimiter="\t")]
        assert names == ["Escherichia coli", "Unclassified"]

    def test_counts_are_untouched(self, tmp_path):
        p = tmp_path / "t.tsv"
        combined_table(p, "species", [("Escherichia coli", "90", "80"),
                                      ("", "10", "20")])
        before = p.read_text().split("\t")
        subprocess.run([sys.executable, str(SCRIPTS / "label_unclassified.py"),
                        "species", str(p)], check=True, capture_output=True)
        with p.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert sum(float(r["barcode01"]) for r in rows) == 100
        assert sum(float(r["barcode02"]) for r in rows) == 100
        assert len(before) == len(p.read_text().split("\t"))

    def test_running_twice_changes_nothing(self, tmp_path):
        p = tmp_path / "t.tsv"
        combined_table(p, "species", [("Escherichia coli", "90", "80"),
                                      ("", "10", "20")])
        for _ in range(2):
            subprocess.run([sys.executable, str(SCRIPTS / "label_unclassified.py"),
                            "species", str(p)], check=True, capture_output=True)
        once = p.read_text()
        subprocess.run([sys.executable, str(SCRIPTS / "label_unclassified.py"),
                        "species", str(p)], check=True, capture_output=True)
        assert p.read_text() == once

    def test_table_without_an_unnamed_row_is_left_alone(self, tmp_path):
        p = tmp_path / "t.tsv"
        combined_table(p, "species", [("Escherichia coli", "90", "80")])
        before = p.read_text()
        subprocess.run([sys.executable, str(SCRIPTS / "label_unclassified.py"),
                        "species", str(p)], check=True, capture_output=True)
        assert p.read_text() == before

    @pytest.mark.parametrize("rank", ["species", "genus", "phylum"])
    def test_every_rank(self, tmp_path, rank):
        p = tmp_path / f"emu-combined-{rank}-counts.tsv"
        combined_table(p, rank, [("Something", "90", "80"), ("", "10", "20")])
        subprocess.run([sys.executable, str(SCRIPTS / "label_unclassified.py"),
                        rank, str(p)], check=True, capture_output=True)
        with p.open() as fh:
            assert [r[rank] for r in csv.DictReader(fh, delimiter="\t")][-1] == "Unclassified"


class TestPerBarcodeTaxa:
    """Which species, not just how many -- and the counts still add up."""

    def test_lists_each_species_with_its_reads(self, tmp_path):
        acct, taxa = run_accounting(
            tmp_path,
            {"barcode01": (100, 100)},
            {"barcode01": [("1", "Escherichia coli", "Escherichia", 60),
                           ("2", "Bacillus subtilis", "Bacillus", 30),
                           ("unmapped", "", "", 10)]},
            want_taxa=True,
        )
        rows = [r for r in taxa if r["barcode"] == "barcode01"]
        assert [r["species"] for r in rows] == ["Escherichia coli", "Bacillus subtilis"]
        assert [int(r["reads"]) for r in rows] == [60, 30]
        assert [r["genus"] for r in rows] == ["Escherichia", "Bacillus"]

    def test_unclassified_reads_are_not_listed_as_a_species(self, tmp_path):
        _, taxa = run_accounting(
            tmp_path,
            {"barcode01": (100, 100)},
            {"barcode01": [("1", "Escherichia coli", "Escherichia", 90),
                           ("mapped_unclassified", "", "", 10)]},
            want_taxa=True,
        )
        assert [r["species"] for r in taxa] == ["Escherichia coli"]

    def test_reads_sum_exactly_to_the_classified_total(self, tmp_path):
        """Emu's counts are fractional; rounding each row on its own would
        leave the column a read or two off the total this table exists to
        reconcile."""
        acct, taxa = run_accounting(
            tmp_path,
            {"barcode01": (100, 100)},
            {"barcode01": [("1", "Escherichia coli", "Escherichia", 33.3),
                           ("2", "Bacillus subtilis", "Bacillus", 33.3),
                           ("3", "Listeria monocytogenes", "Listeria", 33.4)]},
            want_taxa=True,
        )
        assert sum(int(r["reads"]) for r in taxa) == int(acct["barcode01"]["reads_classified"])

    def test_species_count_matches_the_accounting_table(self, tmp_path):
        acct, taxa = run_accounting(
            tmp_path,
            {"barcode01": (100, 100), "barcode02": (50, 50)},
            {"barcode01": [("1", "Escherichia coli", "Escherichia", 50),
                           ("2", "Bacillus subtilis", "Bacillus", 50)],
             "barcode02": [("1", "Escherichia coli", "Escherichia", 50)]},
            want_taxa=True,
        )
        for bc in ("barcode01", "barcode02"):
            listed = [r for r in taxa if r["barcode"] == bc]
            assert len(listed) == int(acct[bc]["species_found"])

