"""When a barcode loses every read, the report names the setting responsible.

The filter applies two things: a length window and a quality floor. The message
named only the window, so a run that lost everything to `--min-quality` sent the
reader to change the setting that was not the problem -- and the raw medians
needed to tell the two apart were already in the summary it reads.

Run with:
    python -m pytest test/ -v
"""

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import make_report  # noqa: E402


class Stub:
    def __init__(self, **kw):
        self.__dict__.update(kw)
        self._items = list(kw.values())

    def __getitem__(self, i):
        return self._items[i]


def render(tmp_path, *, raw_q=None, raw_len=None,
           min_quality=10, min_length=1000, max_length=2000):
    """A one-barcode run that kept nothing, with the given raw medians."""
    summary = tmp_path / "preprocessing_summary.csv"
    cols = ["barcode", "raw_reads", "filtered_reads"]
    vals = ["barcode01", 1000, 0]
    if raw_q is not None:
        cols.append("raw_median_quality")
        vals.append(raw_q)
    if raw_len is not None:
        cols.append("raw_median_length")
        vals.append(raw_len)
    with open(summary, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        w.writerow(vals)

    paths = {}
    for rank, taxon in (("species", "Escherichia coli"), ("genus", "Escherichia")):
        p = tmp_path / f"{rank}.tsv"
        p.write_text(f"{rank}\tbarcode01\n{taxon}\t1.0\n", encoding="utf-8")
        paths[rank] = str(p)

    out = tmp_path / "nano16s_report.html"
    make_report.snakemake = Stub(
        input=Stub(summary=str(summary), species=paths["species"],
                   genus=paths["genus"]),
        params=Stub(db=str(tmp_path / "db"), min_length=min_length,
                    max_length=max_length, min_quality=min_quality,
                    version="1.2.1"),
        output=[str(out)],
    )
    try:
        make_report.main()
    finally:
        del make_report.snakemake
    return out.read_text(encoding="utf-8")


def test_the_quality_floor_is_named_when_it_is_the_cause(tmp_path):
    """REGRESSION: this said "probably outside the 1000-2000 bp window" for a
    barcode whose reads were the right length and too noisy."""
    html = render(tmp_path, raw_q=15.0, raw_len=1600, min_quality=40)
    assert "median raw quality is Q15.0" in html
    assert "bp window" not in html.split("lost every read")[1][:200]


def test_the_length_window_is_named_when_it_is_the_cause(tmp_path):
    html = render(tmp_path, raw_q=15.0, raw_len=1612,
                  min_quality=10, min_length=4000, max_length=5000)
    assert "median raw read is 1,612 bp" in html
    assert "4000-5000 bp window" in html


def test_both_are_named_when_the_medians_do_not_say(tmp_path):
    """No medians in the summary, or both within range: name both settings
    rather than guessing one."""
    html = render(tmp_path)
    after = html.split("lost every read")[1][:300]
    assert "bp window" in after
    assert "below Q" in after


def test_a_barcode_with_no_reads_at_all_is_unchanged(tmp_path):
    """The filter never saw them, so neither setting is the explanation."""
    summary = tmp_path / "preprocessing_summary.csv"
    with open(summary, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["barcode", "raw_reads", "filtered_reads"])
        w.writerow(["barcode01", 0, 0])
    paths = {}
    for rank, taxon in (("species", "Escherichia coli"), ("genus", "Escherichia")):
        p = tmp_path / f"{rank}.tsv"
        p.write_text(f"{rank}\tbarcode01\n{taxon}\t1.0\n", encoding="utf-8")
        paths[rank] = str(p)
    out = tmp_path / "nano16s_report.html"
    make_report.snakemake = Stub(
        input=Stub(summary=str(summary), species=paths["species"],
                   genus=paths["genus"]),
        params=Stub(db=str(tmp_path / "db"), min_length=1000, max_length=2000,
                    min_quality=10, version="1.2.1"),
        output=[str(out)],
    )
    try:
        make_report.main()
    finally:
        del make_report.snakemake
    html = out.read_text(encoding="utf-8")
    assert "had no reads at all" in html
