"""Unit tests for the performance report.

The report exists to answer two questions: where did the time go, and does
any barcode look wrong. Both answers are computed from files the pipeline
leaves behind, so the failure modes are parsing ones — a missing column, a
job that never ran, a barcode still in flight. Those are what is covered here.

The threshold tests matter most. An earlier version of this report flagged a
barcode only when its retention fell below half the run median, which cannot
fire when a whole run is uniformly bad and stayed silent on a barcode that
kept 60% of its reads. Every threshold below is absolute for that reason.

Run with:
    python -m pytest test/ -v
"""

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from make_perf_report import (  # noqa: E402
    MIN_MEDIAN_Q, MIN_READS, MIN_RETENTION_PCT,
    build_rows, dur, find_flags, parse_nanostat, read_benchmarks,
    stage_bars_svg, timeline_svg,
)

BENCH_HEADER = [
    "s", "h:m:s", "max_rss", "max_vms", "max_uss", "max_pss",
    "io_in", "io_out", "mean_load", "cpu_time",
]

CFG = {"min_length": 1000, "max_length": 2000, "min_quality": 10}


def write_bench(root, rule, sample, seconds, cpu=None, rss=None, rows=None):
    """Write a Snakemake-shaped benchmark TSV."""
    path = root / rule / f"{sample}.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    body = rows if rows is not None else [[
        seconds, "0:00:00", rss if rss is not None else 100.0,
        0, 0, 0, 0, 0, 1.0, cpu if cpu is not None else seconds,
    ]]
    lines = ["\t".join(BENCH_HEADER)]
    lines += ["\t".join(str(c) for c in r) for r in body]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_nanostat(root, sample, reads, *, filtered=False, median_q=17.5,
                   median_len=1450):
    suffix = "_filtered_quality_summary.txt" if filtered else "_quality_summary.txt"
    path = root / sample / f"{sample}{suffix}"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        f"Number of reads:\t{reads:,.1f}\n"
        f"Total bases:\t{reads * 1500:,.1f}\n"
        f"Median read length:\t{median_len:,.1f}\n"
        f"Mean read length:\t{median_len:,.1f}\n"
        f"Read length N50:\t{median_len:,.1f}\n"
        f"Median read quality:\t{median_q}\n"
        f"Mean read quality:\t{median_q}\n"
        f">Q10:\t{int(reads * 0.98)} (98.0%) 5.7Mb\n"
        f">Q15:\t{int(reads * 0.77)} (77.0%) 4.4Mb\n"
        f">Q20:\t{int(reads * 0.10)} (10.0%) 0.6Mb\n",
        encoding="utf-8",
    )
    return path


# ---------------------------------------------------------------------------
# read_benchmarks
# ---------------------------------------------------------------------------

def test_benchmarks_are_read_with_rule_and_sample_from_the_path(tmp_path):
    write_bench(tmp_path, "porechop", "barcode01", 12.5, cpu=90.0, rss=101.0)
    jobs = read_benchmarks(str(tmp_path))
    assert len(jobs) == 1
    assert jobs[0]["rule"] == "porechop"
    assert jobs[0]["sample"] == "barcode01"
    assert jobs[0]["wall"] == 12.5
    assert jobs[0]["cpu"] == 90.0
    assert jobs[0]["max_rss"] == 101.0


def test_cpu_time_may_exceed_wall_time(tmp_path):
    """A threaded stage does more CPU-seconds than it spends elapsed.

    This is the whole reason for using benchmarks over log timestamps, so
    guard against a future "sanity" clamp quietly discarding it.
    """
    write_bench(tmp_path, "emu", "barcode01", 10.0, cpu=75.0)
    job = read_benchmarks(str(tmp_path))[0]
    assert job["cpu"] > job["wall"]


def test_start_is_derived_from_mtime_minus_duration(tmp_path):
    """The TSV holds a duration, not a clock time; the timeline needs both."""
    path = write_bench(tmp_path, "emu", "barcode01", 30.0)
    os.utime(path, (1_000_000, 1_000_000))
    job = read_benchmarks(str(tmp_path))[0]
    assert job["end"] == 1_000_000
    assert job["start"] == 1_000_000 - 30.0


def test_a_retried_job_uses_its_last_attempt(tmp_path):
    """Snakemake appends a row per attempt; the last one is what succeeded."""
    write_bench(tmp_path, "emu", "barcode01", None, rows=[
        [5.0, "0:00:05", 50.0, 0, 0, 0, 0, 0, 1.0, 5.0],
        [42.0, "0:00:42", 90.0, 0, 0, 0, 0, 0, 1.0, 300.0],
    ])
    job = read_benchmarks(str(tmp_path))[0]
    assert job["wall"] == 42.0
    assert job["cpu"] == 300.0


def test_missing_directory_and_unparseable_rows_are_survivable(tmp_path):
    assert read_benchmarks(str(tmp_path / "nope")) == []
    write_bench(tmp_path, "emu", "barcode01", None, rows=[["", "", "", "", "",
                                                           "", "", "", "", ""]])
    write_bench(tmp_path, "emu", "barcode02", 4.0, cpu=None, rows=[
        [4.0, "0:00:04", "-", 0, 0, 0, 0, 0, 1.0, "nan"],
    ])
    jobs = read_benchmarks(str(tmp_path))
    # barcode01 has no usable wall time and is dropped; barcode02 survives
    # with the unreadable fields left as None rather than zero.
    assert [j["sample"] for j in jobs] == ["barcode02"]
    assert jobs[0]["cpu"] is None
    assert jobs[0]["max_rss"] is None


# ---------------------------------------------------------------------------
# parse_nanostat
# ---------------------------------------------------------------------------

def test_quality_cutoff_block_is_recovered(tmp_path):
    """preprocessing_summary.csv drops these; the report needs them."""
    write_nanostat(tmp_path, "barcode01", 1000)
    stats = parse_nanostat(
        str(tmp_path / "barcode01" / "barcode01_quality_summary.txt")
    )
    assert stats["pct_above_Q15"] == 77.0
    assert stats["reads_above_Q10"] == 980
    assert stats["Number of reads"] == 1000.0


def test_missing_nanostat_file_is_empty_not_an_error(tmp_path):
    assert parse_nanostat(str(tmp_path / "absent.txt")) == {}


# ---------------------------------------------------------------------------
# build_rows
# ---------------------------------------------------------------------------

def test_retention_and_throughput_are_computed(tmp_path):
    raw, filt = tmp_path / "raw", tmp_path / "filt"
    write_nanostat(raw, "barcode01", 10000)
    write_nanostat(filt, "barcode01", 9500, filtered=True)
    per_stage = {"barcode01": {
        "porechop": {"wall": 100.0, "cpu": 700.0, "max_rss": 120.0},
        "emu": {"wall": 50.0, "cpu": 350.0, "max_rss": 480.0},
    }}
    row = build_rows(["barcode01"], str(raw), str(filt), per_stage)[0]
    assert row["retention_pct"] == 95.0
    assert row["total_s"] == 150.0
    assert row["total_cpu_s"] == 1050.0
    assert row["max_rss_mb"] == 480.0
    assert row["sec_per_1k_reads"] == 15.0


def test_a_barcode_still_in_flight_does_not_divide_by_none(tmp_path):
    """Raw stats exist before filtered ones; the report must survive that."""
    raw, filt = tmp_path / "raw", tmp_path / "filt"
    write_nanostat(raw, "barcode01", 5000)
    row = build_rows(["barcode01"], str(raw), str(filt), {})[0]
    assert row["retention_pct"] is None
    assert row["total_s"] is None
    assert row["raw_reads"] == 5000


def test_read_counts_are_integers_not_floats(tmp_path):
    """NanoStat prints '3,629.0'; these are counts and should read as counts."""
    raw, filt = tmp_path / "raw", tmp_path / "filt"
    write_nanostat(raw, "barcode01", 3629)
    write_nanostat(filt, "barcode01", 3606, filtered=True)
    row = build_rows(["barcode01"], str(raw), str(filt), {})[0]
    assert row["raw_reads"] == 3629
    assert isinstance(row["raw_reads"], int)


# ---------------------------------------------------------------------------
# find_flags — absolute thresholds
# ---------------------------------------------------------------------------

def row(**kw):
    base = {
        "barcode": "barcode01", "raw_reads": 10000, "retention_pct": 99.0,
        "filtered_median_len": 1450, "filtered_median_q": 17.5,
    }
    base.update(kw)
    return base


def test_uniformly_poor_retention_is_still_flagged():
    """REGRESSION: a median-relative check cannot fire when every barcode is bad.

    Three barcodes at 60% retention have a median of 60%, so nothing sits
    below half of it and the old check reported "nothing to see" on a run
    where the length window plainly did not match the amplicon.
    """
    rows = [row(barcode=f"barcode0{i}", retention_pct=60.0) for i in (1, 2, 3)]
    flags = find_flags(rows, CFG)
    assert len(flags) == 3
    assert all(level == "bad" for _, level, _ in flags)
    assert "length window" in flags[0][2]


def test_retention_just_above_the_floor_is_quiet():
    assert find_flags([row(retention_pct=MIN_RETENTION_PCT + 0.1)], CFG) == []


def test_low_depth_is_flagged_on_its_own_terms():
    flags = find_flags([row(raw_reads=MIN_READS - 1)], CFG)
    assert len(flags) == 1
    assert "raw reads" in flags[0][2]


def test_a_barcode_far_below_the_run_median_is_flagged_even_above_the_floor():
    rows = [row(barcode="barcode01", raw_reads=2000)]
    rows += [row(barcode=f"barcode0{i}", raw_reads=40000) for i in (2, 3, 4)]
    flags = find_flags(rows, CFG)
    assert [b for b, _, _ in flags] == ["barcode01"]
    assert "below the run median" in flags[0][2]


def test_median_length_outside_the_configured_window_is_flagged():
    flags = find_flags([row(filtered_median_len=650)], CFG)
    assert len(flags) == 1
    assert "outside the configured window" in flags[0][2]


def test_low_filtered_quality_is_flagged():
    flags = find_flags([row(filtered_median_q=MIN_MEDIAN_Q - 1)], CFG)
    assert any("median quality" in m for _, _, m in flags)


def test_a_clean_run_produces_no_flags():
    rows = [row(barcode=f"barcode{i:02d}") for i in range(1, 7)]
    assert find_flags(rows, CFG) == []


def test_a_barcode_with_no_data_is_not_flagged_for_every_threshold():
    """An empty barcode should not produce four separate complaints."""
    flags = find_flags([row(raw_reads=None, retention_pct=None,
                            filtered_median_len=None,
                            filtered_median_q=None)], CFG)
    assert flags == []


# ---------------------------------------------------------------------------
# formatting and charts
# ---------------------------------------------------------------------------

def test_sub_second_durations_survive():
    """Log scraping rounded these to '0s'; benchmarks resolve them."""
    assert dur(0.4) == "0.4s"
    assert dur(None) == "-"
    assert dur(90) == "1m 30s"
    assert dur(3725) == "1h 2m 5s"


def test_charts_handle_no_data():
    assert "No benchmark data" in stage_bars_svg([])
    svg, legend = timeline_svg([], [])
    assert "No timing data" in svg
    assert legend == ""


def test_timeline_places_one_segment_per_job():
    jobs = [
        {"rule": "porechop", "sample": "barcode01", "wall": 10.0,
         "start": 0.0, "end": 10.0},
        {"rule": "emu", "sample": "barcode01", "wall": 5.0,
         "start": 10.0, "end": 15.0},
        {"rule": "porechop", "sample": "barcode02", "wall": 8.0,
         "start": 0.0, "end": 8.0},
    ]
    svg, legend = timeline_svg(jobs, ["porechop", "emu"])
    assert svg.count("<rect") == 3
    assert legend.count("<i class=") == 2


def test_timeline_skips_jobs_with_no_recoverable_start():
    jobs = [{"rule": "emu", "sample": "barcode01", "wall": 5.0,
             "start": None, "end": None}]
    svg, _ = timeline_svg(jobs, ["emu"])
    assert "No timing data" in svg
