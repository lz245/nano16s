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
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import make_perf_report  # noqa: E402
from make_perf_report import (  # noqa: E402
    MIN_MEDIAN_Q, MIN_READS, MIN_RETENTION_PCT,
    build_rows, dur, ffix, find_flags, fint, parse_nanostat, read_benchmarks,
    split_sittings, stage_bars_svg, system_info, timeline_svg,
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
        "filtered_reads": 9900,
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
                            filtered_reads=None,
                            filtered_median_len=None,
                            filtered_median_q=None)], CFG)
    assert flags == []


def test_a_barcode_the_filter_emptied_is_not_given_a_quality():
    """REGRESSION: "filtered median quality is Q0.0" on a barcode with no reads.

    NanoStat writes zeros for an empty file rather than omitting the fields,
    so a barcode whose reads were all filtered out arrived here looking like a
    measurement of zero. Retention already says what happened, and names the
    length window as the likely cause; a second flag asserting a quality that
    was never measured only obscures it. This is the state a run with a length
    window that does not match the amplicon puts every barcode into.
    """
    flags = find_flags([row(retention_pct=0.0, filtered_reads=0,
                            filtered_median_len=0, filtered_median_q=0.0)], CFG)
    assert len(flags) == 1
    assert "survived filtering" in flags[0][2]
    assert not any("median quality" in m for _, _, m in flags)


# ---------------------------------------------------------------------------
# system_info
# ---------------------------------------------------------------------------

def test_system_info_reports_the_keys_the_report_renders():
    """Every key the page reads must exist, even where the value is unknown.

    The report indexes these directly; a missing key is a KeyError that takes
    the whole report down rather than leaving one field blank.
    """
    info = system_info()
    for key in ("hostname", "cpu_model", "cpu_cores_physical",
                "cpu_threads_logical", "ram_mb", "os", "platform", "kernel",
                "arch", "wsl", "python"):
        assert key in info, f"system_info() is missing {key}"


def test_system_info_values_are_sane_here():
    info = system_info()
    assert info["cpu_threads_logical"] >= 1
    assert info["ram_mb"] > 0
    assert info["python"].count(".") >= 1


def test_system_info_survives_a_machine_with_no_proc(monkeypatch):
    """Nothing about the hardware is worth failing a finished run over.

    /proc is absent on macOS and inside some containers; sysctl is absent on
    Linux. Both lookups failing must degrade to None, not raise.
    """
    def no_open(*args, **kwargs):
        raise OSError("no /proc here")

    monkeypatch.setattr("builtins.open", no_open)
    monkeypatch.setattr(make_perf_report, "_sysctl", lambda name: None)
    monkeypatch.setattr(make_perf_report.os, "sysconf",
                        lambda name: (_ for _ in ()).throw(ValueError()))

    info = system_info()
    assert info["ram_mb"] is None
    assert info["cpu_cores_physical"] is None
    # os.cpu_count() and platform do not touch /proc, so these still answer.
    assert info["cpu_threads_logical"] >= 1
    assert info["python"]


def test_a_mac_whose_version_plist_will_not_open(monkeypatch):
    """The same promise as the test above, on the branch only macOS reaches.

    That test passed on Linux and failed on every macOS runner: with `open`
    stubbed out it reached platform.mac_ver(), which reads
    /System/Library/CoreServices/SystemVersion.plist, and the OSError went
    straight through _os_name() and out of system_info(). Forcing the branch
    here means the guard is checked wherever the tests run, not only on a Mac.
    """
    def boom():
        raise OSError("plist unreadable")

    monkeypatch.setattr(make_perf_report.platform, "system", lambda: "Darwin")
    monkeypatch.setattr(make_perf_report.platform, "mac_ver", boom)

    assert make_perf_report._os_name() == "macOS"
    info = system_info()
    assert info["os"] == "macOS"
    assert info["wsl"] is False


def test_wsl_is_reported_as_both_windows_and_linux(monkeypatch):
    """"Win or Linux" has a third answer, and it is the one in use here."""
    monkeypatch.setattr(make_perf_report, "_is_wsl", lambda: True)
    monkeypatch.setattr(make_perf_report.platform, "system", lambda: "Linux")
    monkeypatch.setattr(make_perf_report, "_os_name", lambda: "Ubuntu 26.04 LTS")
    info = system_info()
    assert info["wsl"] is True
    assert "Windows" in info["platform"]
    assert "Ubuntu" in info["platform"]


def test_native_linux_is_not_labelled_windows(monkeypatch):
    monkeypatch.setattr(make_perf_report, "_is_wsl", lambda: False)
    monkeypatch.setattr(make_perf_report.platform, "system", lambda: "Linux")
    monkeypatch.setattr(make_perf_report, "_os_name", lambda: "Ubuntu 24.04 LTS")
    info = system_info()
    assert info["wsl"] is False
    assert "Windows" not in info["platform"]


# ---------------------------------------------------------------------------
# formatting and charts
# ---------------------------------------------------------------------------

def test_sub_second_durations_survive():
    """Log scraping rounded these to '0s'; benchmarks resolve them."""
    assert dur(0.4) == "0.4s"
    assert dur(None) == "-"
    assert dur(90) == "1m 30s"
    assert dur(3725) == "1h 2m 5s"


def test_a_recorded_zero_is_not_shown_as_a_measurement():
    """0.0 means "finished inside the sampling interval", not "took no time".

    Printed as "0.0s" beside a non-zero wall time it reads as a broken
    number, which is how the chopper row looked.
    """
    assert dur(0) == "<0.1s"
    assert dur(0.0) == "<0.1s"


def test_read_counts_lose_the_nanostat_decimal():
    """REGRESSION: read counts rendered as '993.0' and '1,059.0'.

    NanoStat prints counts as floats. Formatting them with one decimal place
    put a tenth of a read in every row of the table.
    """
    assert fint(993.0) == "993"
    assert fint(1059) == "1,059"
    assert fint(None) == "-"


def test_measurements_keep_a_fixed_precision():
    """REGRESSION: a quality column mixing '15', '15.1' and '18' misaligns.

    Integral floats were shortened to bare integers, so decimal points did
    not line up down the column.
    """
    assert ffix(15.0) == "15.0"
    assert ffix(15.14) == "15.1"
    assert ffix(18) == "18.0"
    assert ffix(None) == "-"


def test_chart_scale_spans_cpu_time_when_it_exceeds_wall():
    """REGRESSION: bars ran off the panel and were clipped.

    The x-axis was scaled to peak wall time while the CPU bar was drawn on
    the same scale. Any threaded stage spends more CPU-seconds than elapsed
    seconds, so porechop and emu overflowed the plot area entirely.
    """
    rollup = [
        {"stage": "porechop", "jobs": 6, "wall_total": 102.0,
         "cpu_total": 239.0, "wall_mean": 17.0, "wall_max": 22.0,
         "max_rss_mb": 100.0, "sec_per_1k_reads": 15.3},
        {"stage": "chopper", "jobs": 6, "wall_total": 0.6,
         "cpu_total": 0.0, "wall_mean": 0.1, "wall_max": 0.1,
         "max_rss_mb": 9.0, "sec_per_1k_reads": 0.09},
    ]
    svg = stage_bars_svg(rollup)
    # Widths are relative to the plot area, which is width - pad_l - pad_r.
    plot_w = 720 - 118 - 84
    widths = [float(w) for w in re.findall(r'width="([\d.]+)"', svg)]
    assert widths, "expected bars to be drawn"
    assert max(widths) <= plot_w + 0.5, (
        f"a bar of {max(widths):.1f}px overflows the {plot_w}px plot area"
    )


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


# ---------------------------------------------------------------------------
# the whole report
# ---------------------------------------------------------------------------

class FakeParams(dict):
    """Snakemake exposes params as attributes; a dict subclass is enough."""

    __getattr__ = dict.__getitem__


class FakeSnakemake:
    def __init__(self, params, output):
        self.params = FakeParams(params)
        self.output = FakeParams(output)


def run_report(tmp_path, *, benchmarks=True, bench_dir=None):
    """Generate a complete report and return its HTML.

    main() reads nothing from snakemake.input -- that exists only to order
    the DAG -- so a params/output stub drives the whole thing.
    """
    raw, filt = tmp_path / "raw", tmp_path / "filt"
    bench = Path(bench_dir) if bench_dir else tmp_path / "bench"
    for sample in ("barcode01", "barcode02"):
        write_nanostat(raw, sample, 5000)
        write_nanostat(filt, sample, 4800, filtered=True)
        if benchmarks:
            write_bench(bench, "porechop", sample, 20.0, cpu=60.0, rss=100.0)
            write_bench(bench, "emu", sample, 10.0, cpu=50.0, rss=500.0)

    out = {
        "html": str(tmp_path / "performance_report.html"),
        "csv": str(tmp_path / "performance_summary.csv"),
        "json": str(tmp_path / "performance.json"),
    }
    make_perf_report.snakemake = FakeSnakemake(
        {
            "bench_dir": str(bench), "raw_dir": str(raw),
            "filtered_dir": str(filt), "cores": 8,
            "min_length": 1000, "max_length": 2000, "min_quality": 10,
            "db": "ncbi_16s", "version": "1.1.0",
        },
        out,
    )
    try:
        make_perf_report.main()
    finally:
        del make_perf_report.snakemake
    return Path(out["html"]).read_text(encoding="utf-8")


def test_a_normal_run_reports_where_the_time_went(tmp_path):
    """Guards the test below: these sections are present when timings exist."""
    html = run_report(tmp_path)
    assert "<h2>Where the time went</h2>" in html
    assert "<h2>Machine use</h2>" in html
    assert "<h2>No timing data</h2>" not in html


def test_a_run_with_no_benchmarks_explains_itself(tmp_path):
    """REGRESSION: the report went hollow without saying why.

    Re-running into a directory that already holds finished results leaves
    Snakemake nothing to run, so no job records a benchmark. The QC half of
    the report still renders from NanoStat, so the page looked complete while
    'Machine use' and 'Where the time went' had silently vanished and every
    timing read '-'. That is indistinguishable from a broken report.
    """
    html = run_report(tmp_path, benchmarks=False)

    # The timing sections are genuinely gone -- there is nothing to put in
    # them -- so the report has to account for their absence. Match the
    # heading, not the words: "No timing data" also appears in the timeline
    # placeholder, which was already there while the report was still silent.
    assert "<h2>Where the time went</h2>" not in html
    assert "<h2>Machine use</h2>" not in html
    assert "<h2>No timing data</h2>" in html
    assert "--forceall" in html

    # And the half that does not depend on timings is still there.
    assert "barcode01" in html
    assert "Per barcode" in html


# ---------------------------------------------------------------------------
# split_sittings
# ---------------------------------------------------------------------------

def j(start, end):
    return {"start": float(start), "end": float(end)}


def test_one_uninterrupted_run_is_one_sitting():
    assert split_sittings([j(0, 30), j(10, 50), j(45, 90)]) == [(0.0, 90.0)]


def test_overlapping_jobs_do_not_split_a_sitting():
    """Jobs run concurrently, so a sitting ends only when nothing is running."""
    assert split_sittings([j(0, 100), j(20, 40), j(60, 80)]) == [(0.0, 100.0)]


def test_a_week_long_gap_is_two_sittings():
    """REGRESSION: a benchmark is only rewritten by a job that runs.

    Resuming a directory leaves the earlier stages' records in place, so the
    earliest-to-latest span measured the calendar gap between two runs rather
    than a run. A 39-minute resume of a week-old directory reported 144 hours
    elapsed and 0.5% core use, with nothing to explain it.
    """
    week = 7 * 24 * 3600
    got = split_sittings([j(0, 600), j(week, week + 120)])
    assert got == [(0.0, 600.0), (float(week), float(week + 120))]
    active = sum(e - s for s, e in got)
    assert active == 720.0                      # not a week


def test_a_short_pause_is_not_a_new_sitting():
    """Snakemake starts the next job as a core frees; seconds are normal."""
    assert split_sittings([j(0, 100), j(160, 200)]) == [(0.0, 200.0)]


def test_no_placed_jobs_is_no_sittings():
    assert split_sittings([]) == []


def test_a_resumed_run_says_so_in_the_report(tmp_path):
    """The banner must name the situation, not leave the reader to infer it."""
    bench = tmp_path / "benchmarks"
    old = write_bench(bench, "porechop", "barcode01", 100.0, cpu=200.0)
    new = write_bench(bench, "emu", "barcode01", 50.0, cpu=300.0)
    os.utime(old, (1_000_000, 1_000_000))
    os.utime(new, (1_000_000 + 7 * 24 * 3600, 1_000_000 + 7 * 24 * 3600))
    html = run_report(tmp_path, benchmarks=False, bench_dir=str(bench))
    assert "<h2>Timings from more than one run</h2>" in html
    assert "2 separate runs" in html


def test_the_timeline_does_not_stretch_across_a_gap_between_runs(tmp_path):
    """REGRESSION: the chart was drawn against the calendar, not the work.

    A resumed directory holds records from several sittings. Spanning the axis
    from the earliest to the latest left days of empty chart and squeezed every
    bar to the 1.2px floor, so all 144 of them rendered identically and none
    showed its own duration.
    """
    bench = tmp_path / "b"
    a = write_bench(bench, "porechop", "barcode01", 100.0, cpu=200.0)
    b = write_bench(bench, "emu", "barcode01", 100.0, cpu=300.0)
    os.utime(a, (1_000_000, 1_000_000))
    os.utime(b, (1_000_000 + 7 * 24 * 3600, 1_000_000 + 7 * 24 * 3600))
    jobs = read_benchmarks(str(bench))
    svg, _ = timeline_svg(jobs, ["porechop", "emu"])

    widths = [float(w) for w in re.findall(r'class="seg[^"]*"[^>]*width="([\d.]+)"', svg)]
    assert len(widths) == 2
    # Two 100s jobs over 200s of working time: each fills about half the plot,
    # not the 1.2px floor it collapsed to when the axis was a week wide.
    assert min(widths) > 100, f"bars collapsed: {widths}"
    # And the axis tops out at the working time, not the week.
    assert "7d" not in svg and "168h" not in svg
