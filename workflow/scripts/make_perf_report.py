#!/usr/bin/env python3
"""Build a self-contained HTML performance report for a nano16s run.

Companion to make_report.py, which covers the biology. This one covers the
run itself: where the time went, how well the machine was used, and which
barcodes look wrong on QC grounds.

Timings come from Snakemake's `benchmark:` directive, not from scraping the
run log. Each benchmarked job leaves a TSV holding wall seconds, CPU seconds
and peak RSS. That matters for three reasons:

  * CPU time is recorded, so "8 threads did 2 seconds of work in 1 second of
    wall clock" is visible. Log timestamps can only ever show wall clock.
  * There is nothing to parse. Log scraping depends on the C locale for
    month and day names, on the log living at a guessable path, and on
    picking the right file when two runs overlap.
  * Resumed runs keep the benchmark files of jobs that did not re-run, so
    the report describes the output directory rather than the last invocation.

Driven by Snakemake (`snakemake.input`, `snakemake.output`, `snakemake.params`).

Note: no `from __future__ import annotations` here. Snakemake's `script:`
directive prepends its own preamble to this file before running it, which
would push a __future__ import off line 1 and make it a SyntaxError.
"""

import csv
import glob
import html
import json
import os
import platform
import re
import shutil
import socket
import subprocess
from datetime import datetime
from pathlib import Path

# Stages in pipeline order, so every table reads the way the run executed.
STAGE_ORDER = [
    "merge",
    "nanostat_raw",
    "porechop",
    "chopper",
    "emu",
    "nanostat_filtered",
]

# Shared with make_report.py so the two reports look like one product.
PALETTE_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
PALETTE_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500",
                "#d55181", "#008300", "#9085e9", "#e66767"]

# --- QC thresholds ---------------------------------------------------------
# Absolute floors, not percentiles. A relative-only check cannot fire when a
# whole run is uniformly bad, which is the case most worth catching.
MIN_RETENTION_PCT = 80.0    # below this, the length window likely mismatches
MIN_READS = 1000            # below this, abundances are not worth quoting
LOW_DEPTH_FACTOR = 0.25     # also flag barcodes far below the run's median
MIN_MEDIAN_Q = 12.0         # filtered reads should clear this comfortably


# ---------------------------------------------------------------------------
# the machine
# ---------------------------------------------------------------------------
# Timings mean little without the hardware they were measured on: "8 hours"
# is a different result on 4 cores than on 19. Everything below is read from
# the standard library and /proc, with a short sysctl fallback for macOS, so
# no dependency is added and nothing here can fail the report.

def _sysctl(name):
    """macOS: one sysctl value, or None. Never raises."""
    exe = shutil.which("sysctl")
    if not exe:
        return None
    try:
        out = subprocess.run([exe, "-n", name], capture_output=True, text=True,
                             timeout=5)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def _cpu_model():
    try:
        with open("/proc/cpuinfo", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("model name", "Model name", "Hardware")):
                    return line.split(":", 1)[1].strip()
    except OSError:
        pass
    return _sysctl("machdep.cpu.brand_string") or platform.processor() or None


def _physical_cores():
    """Distinct (socket, core) pairs — hyperthreads are not extra cores."""
    try:
        pairs, cur = set(), {}
        with open("/proc/cpuinfo", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    if "physical id" in cur and "core id" in cur:
                        pairs.add((cur["physical id"], cur["core id"]))
                    cur = {}
                    continue
                if ":" in line:
                    k, _, v = line.partition(":")
                    cur[k.strip()] = v.strip()
            if "physical id" in cur and "core id" in cur:
                pairs.add((cur["physical id"], cur["core id"]))
        if pairs:
            return len(pairs)
    except OSError:
        pass
    mac = _sysctl("hw.physicalcpu")
    return int(mac) if mac and mac.isdigit() else None


def _total_ram_mb():
    try:
        with open("/proc/meminfo", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / 1024.0     # kB -> MB
    except (OSError, ValueError, IndexError):
        pass
    mac = _sysctl("hw.memsize")
    if mac and mac.isdigit():
        return int(mac) / (1024.0 * 1024.0)
    try:                                    # POSIX, present on most Unixes
        return (os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
                / (1024.0 * 1024.0))
    except (ValueError, OSError, AttributeError):
        return None


def _os_name():
    """Distribution or product name, not just the kernel."""
    if platform.system() == "Linux":
        try:
            with open("/etc/os-release", encoding="utf-8") as fh:
                fields = {}
                for line in fh:
                    if "=" in line:
                        k, _, v = line.partition("=")
                        fields[k.strip()] = v.strip().strip('"')
            name = fields.get("PRETTY_NAME") or fields.get("NAME")
            if name:
                return name
        except OSError:
            pass
        return "Linux"
    if platform.system() == "Darwin":
        return f"macOS {platform.mac_ver()[0]}".strip()
    return f"{platform.system()} {platform.release()}".strip()


def _is_wsl():
    """WSL is Linux hosted by Windows; reporting only one of the two misleads."""
    try:
        with open("/proc/version", encoding="utf-8", errors="replace") as fh:
            blob = fh.read().lower()
        return "microsoft" in blob or "wsl" in blob
    except OSError:
        return "microsoft" in platform.release().lower()


def system_info():
    """Hardware and OS the run executed on. Best effort; keys may be None."""
    logical = os.cpu_count()
    physical = _physical_cores()
    wsl = _is_wsl() if platform.system() == "Linux" else False
    osname = _os_name()
    platform_label = osname
    if wsl:
        # Dr. Li asked for "Win or Linux"; on WSL the honest answer is both,
        # and it matters — the filesystem and scheduler are not native Linux.
        platform_label = f"{osname} on Windows (WSL2)"

    return {
        "hostname": socket.gethostname() or None,
        "cpu_model": _cpu_model(),
        "cpu_cores_physical": physical,
        "cpu_threads_logical": logical,
        "ram_mb": _total_ram_mb(),
        "os": osname,
        "platform": platform_label,
        "kernel": platform.release() or None,
        "arch": platform.machine() or None,
        "wsl": wsl,
        "python": platform.python_version(),
    }


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------

def read_benchmarks(bench_dir):
    """benchmarks/{rule}/{sample}.tsv -> list of per-job dicts.

    Snakemake writes the file when the job finishes, so mtime is the end
    time and mtime - wall is the start. That is the only way to place jobs
    on a timeline; the TSV itself holds durations, not clock times.
    """
    jobs = []
    for path in sorted(glob.glob(os.path.join(bench_dir, "*", "*.tsv"))):
        rule = os.path.basename(os.path.dirname(path))
        sample = os.path.splitext(os.path.basename(path))[0]
        try:
            with open(path, newline="", encoding="utf-8") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
        except OSError:
            continue
        if not rows:
            continue

        # Snakemake writes one row per job. A rule that Snakemake retried
        # appends further rows; the last one is the attempt that succeeded.
        row = rows[-1]

        def num(key):
            raw = (row.get(key) or "").strip()
            if raw in ("", "-", "nan"):
                return None
            try:
                return float(raw)
            except ValueError:
                return None

        wall = num("s")
        if wall is None:
            continue
        try:
            end = os.path.getmtime(path)
        except OSError:
            end = None

        jobs.append({
            "rule": rule,
            "sample": sample,
            "wall": wall,
            "cpu": num("cpu_time"),
            "max_rss": num("max_rss"),
            "end": end,
            "start": (end - wall) if end is not None else None,
        })
    return jobs


def parse_nanostat(path):
    """Full NanoStat parse, including the >Q10/>Q15/>Q20 block.

    preprocessing_summary.csv keeps four fields and drops the rest, so the
    quality distribution is only available by reading the source file.
    """
    out = {}
    if not os.path.exists(path):
        return out
    with open(path, errors="replace", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip()
            # ">Q10:\t3587 (98.8%) 5.7Mb"
            q = re.match(r"^>Q(\d+):\s+(\d+)\s+\(([\d.]+)%\)", line)
            if q:
                out[f"reads_above_Q{q.group(1)}"] = int(q.group(2))
                out[f"pct_above_Q{q.group(1)}"] = float(q.group(3))
                continue
            if ":" in line:
                key, _, val = line.partition(":")
                val = val.strip().replace(",", "")
                try:
                    out[key.strip()] = float(val)
                except ValueError:
                    pass
    return out


# ---------------------------------------------------------------------------
# formatting
# ---------------------------------------------------------------------------

def esc(s):
    return html.escape(str(s), quote=True)


def dur(seconds):
    """Compact duration. Sub-second values are real here, unlike log scraping."""
    if seconds is None:
        return "-"
    # Anything that rounds to zero finished inside the sampling interval; it
    # did not take no time. Printing "0.0s" beside a non-zero wall time reads
    # as a broken measurement, so say what is actually known.
    if seconds < 0.05:
        return "<0.1s"
    if seconds < 1:
        return f"{seconds:.1f}s"
    seconds = int(round(seconds))
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h {m}m {s}s"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"


def axis(seconds):
    """Duration for an axis tick, where the origin is genuinely zero.

    dur() reports a rounded-down measurement as "<0.1s", which is right for a
    job that ran and wrong for the left edge of a chart.
    """
    return "0s" if seconds <= 0 else dur(seconds)


def fint(v):
    """Counts and lengths. NanoStat prints '993.0'; a read count has no tenths."""
    if v is None or v == "":
        return "-"
    if isinstance(v, (int, float)):
        return f"{int(round(v)):,}"
    return str(v)


def ffix(v, dp=1):
    """Measurements, always to the same precision.

    Quality scores mixing "15" and "15.1" down a column read as different
    kinds of number and break the decimal alignment that makes a column
    scannable.
    """
    if v is None or v == "":
        return "-"
    if isinstance(v, (int, float)):
        return f"{v:,.{dp}f}"
    return str(v)


def as_int(v):
    """NanoStat prints counts as '3,629.0'; they are counts, not measurements."""
    return int(v) if isinstance(v, float) else v


def pct(v, dp=1):
    return "-" if v is None else f"{v:.{dp}f}%"


# ---------------------------------------------------------------------------
# charts
# ---------------------------------------------------------------------------

def stage_bars_svg(rollup):
    """Wall time and CPU time per stage, as a pair of bars on each row.

    Not nested bars. A threaded stage spends more CPU-seconds than it spends
    elapsed seconds, so CPU is routinely the larger of the two and drawing it
    inside wall time both inverts the encoding and overflows any axis scaled
    to wall time alone. The scale below spans whichever value is larger, and
    the two bars share a baseline so the ratio between them stays readable.
    """
    if not rollup:
        return "<p class='empty'>No benchmark data available.</p>"

    bar_h, pair_gap, row_gap = 9, 3, 12
    row_h = bar_h * 2 + pair_gap
    pad_l, pad_r, pad_t = 118, 84, 28
    width = 720
    plot_w = width - pad_l - pad_r
    height = pad_t + len(rollup) * (row_h + row_gap) + 10

    # Span both series: CPU exceeds wall wherever threads are used.
    peak = max(
        max((r["wall_total"] for r in rollup), default=0),
        max((r["cpu_total"] or 0 for r in rollup), default=0),
    ) or 1

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="Wall and CPU time by pipeline stage">'
    ]
    for frac in (0.25, 0.5, 0.75, 1.0):
        x = pad_l + plot_w * frac
        parts.append(f'<line class="grid" x1="{x:.1f}" y1="{pad_t - 8}" '
                     f'x2="{x:.1f}" y2="{height - 10}"/>')
        parts.append(f'<text class="tick" x="{x:.1f}" y="{pad_t - 12}" '
                     f'text-anchor="middle">{axis(peak * frac)}</text>')

    for i, r in enumerate(rollup):
        y = pad_t + i * (row_h + row_gap)
        wall_w = plot_w * (r["wall_total"] / peak)
        cpu_w = plot_w * ((r["cpu_total"] or 0) / peak)
        parts.append(f'<text class="rowlab" x="{pad_l - 10}" y="{y + row_h - 4}" '
                     f'text-anchor="end">{esc(r["stage"])}</text>')
        parts.append(f'<rect class="raw" x="{pad_l}" y="{y}" '
                     f'width="{max(wall_w, 1):.1f}" height="{bar_h}" rx="2">'
                     f'<title>{esc(r["stage"])}: {dur(r["wall_total"])} wall '
                     f'across {r["jobs"]} jobs</title></rect>')
        parts.append(f'<rect class="kept" x="{pad_l}" y="{y + bar_h + pair_gap}" '
                     f'width="{max(cpu_w, 1):.1f}" height="{bar_h}" rx="2">'
                     f'<title>{esc(r["stage"])}: {dur(r["cpu_total"])} CPU'
                     f'</title></rect>')
        parts.append(f'<text class="rowval" x="{pad_l + plot_w + 8}" '
                     f'y="{y + row_h - 4}">{dur(r["wall_total"])}</text>')
    parts.append("</svg>")
    return "".join(parts)


def timeline_svg(jobs, stages):
    """Gantt of every job, one row per barcode, coloured by stage.

    This is where idle capacity shows up. Aggregate tables cannot distinguish
    a run that saturated the machine from one that trickled through a queue,
    because both report the same totals.
    """
    placed = [j for j in jobs if j.get("start") is not None and j.get("sample")]
    if not placed:
        return ("<p class='empty'>No timing data available for the timeline.</p>", "")

    t0 = min(j["start"] for j in placed)
    t1 = max(j["end"] for j in placed)
    span = max(t1 - t0, 1.0)

    samples = sorted({j["sample"] for j in placed})
    slot = {s: i for i, s in enumerate(stages)}

    row_h, gap, pad_l, pad_r, pad_t = 15, 3, 96, 20, 30
    width = 720
    plot_w = width - pad_l - pad_r
    height = pad_t + len(samples) * (row_h + gap) + 14

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="Job timeline by barcode">'
    ]
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = pad_l + plot_w * frac
        parts.append(f'<line class="grid" x1="{x:.1f}" y1="{pad_t - 8}" '
                     f'x2="{x:.1f}" y2="{height - 14}"/>')
        parts.append(f'<text class="tick" x="{x:.1f}" y="{pad_t - 12}" '
                     f'text-anchor="middle">{axis(span * frac)}</text>')

    for i, sample in enumerate(samples):
        y = pad_t + i * (row_h + gap)
        parts.append(f'<text class="rowlab" x="{pad_l - 10}" y="{y + 11}" '
                     f'text-anchor="end">{esc(sample)}</text>')
        for j in (x for x in placed if x["sample"] == sample):
            x0 = pad_l + plot_w * ((j["start"] - t0) / span)
            w = max(plot_w * (j["wall"] / span), 1.2)
            cls = f"s{slot.get(j['rule'], 0) % len(PALETTE_LIGHT)}"
            parts.append(f'<rect class="seg {cls}" x="{x0:.1f}" y="{y + 2}" '
                         f'width="{w:.1f}" height="{row_h - 4}" rx="2">'
                         f'<title>{esc(sample)} {esc(j["rule"])}: '
                         f'{dur(j["wall"])}</title></rect>')
    parts.append("</svg>")

    legend = ["<div class='legend'>"]
    for stage in stages:
        cls = f"s{slot[stage] % len(PALETTE_LIGHT)}"
        legend.append(f"<span class='lg'><i class='sw {cls}'></i>{esc(stage)}</span>")
    legend.append("</div>")
    return "".join(parts), "".join(legend)


# ---------------------------------------------------------------------------
# page
# ---------------------------------------------------------------------------

# Filled with str.replace, not %-formatting: the template contains literal
# percent signs (width: 100%) that %-formatting would read as conversions.
CSS = """
:root {
  --surface: #fcfcfb; --panel: #ffffff; --line: #e3e2dd;
  --ink: #0b0b0b; --ink-2: #52514e; --ink-3: #78776f;
  --raw: #dedcd4; --kept: #2a78d6;
  --flag: #fff6e0; --flag-line: #eda100; --bad: #e34948;
  @SLOTS_L@
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --surface: #161615; --panel: #1e1e1c; --line: #33322e;
    --ink: #f4f3ef; --ink-2: #c3c2b7; --ink-3: #918f85;
    --raw: #33322e; --kept: #3987e5;
    --flag: #2e2718; --flag-line: #c98500; --bad: #e66767;
    @SLOTS_D@
  }
}
:root[data-theme="dark"] {
  --surface: #161615; --panel: #1e1e1c; --line: #33322e;
  --ink: #f4f3ef; --ink-2: #c3c2b7; --ink-3: #918f85;
  --raw: #33322e; --kept: #3987e5;
  --flag: #2e2718; --flag-line: #c98500; --bad: #e66767;
  @SLOTS_D@
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--surface); color: var(--ink);
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 64rem; margin: 0 auto; padding: 3rem 1.5rem 5rem;
        display: flex; flex-direction: column; gap: 2.25rem; }
h1 { font-size: 1.9rem; margin: 0; letter-spacing: -0.02em; }
h2 { font-size: 1.15rem; margin: 0 0 .25rem; letter-spacing: -0.01em; }
.sub { color: var(--ink-2); margin: .4rem 0 0; }
.note { color: var(--ink-3); font-size: .82rem; margin: .7rem 0 0; }
.panel { background: var(--panel); border: 1px solid var(--line);
         border-radius: 6px; padding: 1.4rem 1.5rem; }
.meta { display: grid; grid-template-columns: repeat(auto-fit, minmax(11rem, 1fr));
        gap: .9rem 1.5rem; }
.meta div { display: flex; flex-direction: column; gap: .15rem; }
.meta dt { font-size: .7rem; letter-spacing: .09em; text-transform: uppercase;
           color: var(--ink-3); }
.meta dd { margin: 0; font-variant-numeric: tabular-nums; word-break: break-word; }
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(9rem, 1fr));
         gap: 1rem; }
.stat { background: var(--panel); border: 1px solid var(--line); border-radius: 6px;
        padding: .9rem 1rem; }
.stat .n { font-size: 1.6rem; font-weight: 600; font-variant-numeric: tabular-nums;
           letter-spacing: -0.02em; }
.stat .k { font-size: .72rem; letter-spacing: .08em; text-transform: uppercase;
           color: var(--ink-3); margin-top: .15rem; }
.stat .h { font-size: .78rem; color: var(--ink-3); margin-top: .3rem; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
.grid { stroke: var(--line); stroke-width: 1; }
.tick, .rowlab, .rowval { font-size: 11px; fill: var(--ink-3); }
.rowlab { font-size: 12px; fill: var(--ink-2); }
.rowval { font-size: 12px; fill: var(--ink-2); font-variant-numeric: tabular-nums; }
rect.raw { fill: var(--raw); }
rect.kept { fill: var(--kept); }
@SLOTCSS@
.legend { display: flex; flex-wrap: wrap; gap: .5rem 1.1rem; margin-top: 1rem;
          font-size: .84rem; color: var(--ink-2); }
.lg { display: inline-flex; align-items: center; gap: .4rem; }
.sw { width: .8rem; height: .8rem; border-radius: 2px; display: inline-block; flex: none; }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: .87rem; margin-top: .4rem; }
th, td { text-align: right; padding: .45rem .6rem; border-bottom: 1px solid var(--line);
         white-space: nowrap; font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
thead th { font-size: .68rem; letter-spacing: .08em; text-transform: uppercase;
           color: var(--ink-3); font-weight: 600; position: sticky; top: 0;
           background: var(--panel); cursor: pointer; user-select: none; }
thead th:hover { color: var(--ink); }
thead th::after { content: "\\2003"; font-size: .7em; color: var(--ink-3); }
thead th[data-dir="1"]::after { content: " \\25B2"; }
thead th[data-dir="-1"]::after { content: " \\25BC"; }
tbody tr:last-child td { border-bottom: none; }
tr.flag { background: var(--flag); }
.warn { border-left: 3px solid var(--flag-line); padding-left: .9rem; margin: .6rem 0; }
.warn.bad { border-left-color: var(--bad); }
.ok { color: var(--ink-2); }
.empty { color: var(--ink-3); font-style: italic; }
code { font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: .86em;
       background: var(--surface); border: 1px solid var(--line);
       border-radius: 3px; padding: .05em .35em; }
footer { color: var(--ink-3); font-size: .8rem; border-top: 1px solid var(--line);
         padding-top: 1rem; }
"""

# Sorting only. Everything else about the page works with scripting disabled.
JS = """
document.querySelectorAll('table.sortable').forEach(function (table) {
  table.querySelectorAll('thead th').forEach(function (th, col) {
    th.addEventListener('click', function () {
      var dir = th.dataset.dir === '1' ? -1 : 1;
      table.querySelectorAll('thead th').forEach(function (o) {
        delete o.dataset.dir;
      });
      th.dataset.dir = String(dir);
      var body = table.tBodies[0];
      var rows = Array.prototype.slice.call(body.rows);
      rows.sort(function (a, b) {
        // data-v carries the unformatted number; text is for humans only.
        var x = a.cells[col].dataset.v, y = b.cells[col].dataset.v;
        var nx = parseFloat(x), ny = parseFloat(y);
        if (!isNaN(nx) && !isNaN(ny)) return (nx - ny) * dir;
        return String(x || a.cells[col].textContent)
          .localeCompare(String(y || b.cells[col].textContent)) * dir;
      });
      rows.forEach(function (r) { body.appendChild(r); });
    });
  });
});
"""


def build_css():
    slots_l = " ".join(f"--s{i}: {c};" for i, c in enumerate(PALETTE_LIGHT))
    slots_d = " ".join(f"--s{i}: {c};" for i, c in enumerate(PALETTE_DARK))
    slotcss = "\n".join(
        f"rect.seg.s{i}, i.sw.s{i} {{ fill: var(--s{i}); background: var(--s{i}); }}"
        for i in range(len(PALETTE_LIGHT))
    )
    return (CSS
            .replace("@SLOTS_L@", slots_l)
            .replace("@SLOTS_D@", slots_d)
            .replace("@SLOTCSS@", slotcss))


def td(value, text=None, cls=""):
    """Cell carrying a sort key in data-v alongside its display text."""
    shown = text if text is not None else value
    sort_key = "" if value is None else value
    klass = f" class='{cls}'" if cls else ""
    return f"<td{klass} data-v='{esc(sort_key)}'>{esc(shown)}</td>"


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------

def build_rows(samples, raw_dir, filtered_dir, per_stage):
    """One record per barcode: QC from NanoStat, timings from benchmarks."""
    rows = []
    for s in samples:
        raw = parse_nanostat(os.path.join(raw_dir, s, f"{s}_quality_summary.txt"))
        flt = parse_nanostat(
            os.path.join(filtered_dir, s, f"{s}_filtered_quality_summary.txt")
        )
        rr = raw.get("Number of reads")
        fr = flt.get("Number of reads")
        # A barcode part-way through the run has raw stats but no filtered
        # ones yet, so guard both sides rather than just the divisor.
        retention = (fr / rr * 100) if (rr and fr is not None) else None

        stages = per_stage.get(s, {})
        wall = {k: v["wall"] for k, v in stages.items()}
        cpu = {k: v["cpu"] for k, v in stages.items() if v["cpu"] is not None}
        rss = [v["max_rss"] for v in stages.values() if v["max_rss"] is not None]

        total_wall = sum(wall.values()) if wall else None
        per_1k = (total_wall / (rr / 1000.0)) if (total_wall and rr) else None

        rows.append({
            "barcode": s,
            "raw_reads": as_int(rr),
            "raw_bases": as_int(raw.get("Total bases")),
            "raw_mean_len": raw.get("Mean read length"),
            "raw_median_len": raw.get("Median read length"),
            "raw_n50": raw.get("Read length N50"),
            "raw_mean_q": raw.get("Mean read quality"),
            "raw_median_q": raw.get("Median read quality"),
            "raw_pct_Q10": raw.get("pct_above_Q10"),
            "raw_pct_Q15": raw.get("pct_above_Q15"),
            "raw_pct_Q20": raw.get("pct_above_Q20"),
            "filtered_reads": as_int(fr),
            "filtered_bases": as_int(flt.get("Total bases")),
            "filtered_mean_len": flt.get("Mean read length"),
            "filtered_median_len": flt.get("Median read length"),
            "filtered_n50": flt.get("Read length N50"),
            "filtered_mean_q": flt.get("Mean read quality"),
            "filtered_median_q": flt.get("Median read quality"),
            "filtered_pct_Q10": flt.get("pct_above_Q10"),
            "filtered_pct_Q15": flt.get("pct_above_Q15"),
            "filtered_pct_Q20": flt.get("pct_above_Q20"),
            "retention_pct": retention,
            "porechop_s": wall.get("porechop"),
            "chopper_s": wall.get("chopper"),
            "emu_s": wall.get("emu"),
            "total_s": total_wall,
            "porechop_cpu_s": cpu.get("porechop"),
            "emu_cpu_s": cpu.get("emu"),
            "total_cpu_s": sum(cpu.values()) if cpu else None,
            "max_rss_mb": max(rss) if rss else None,
            "sec_per_1k_reads": per_1k,
        })
    return rows


def find_flags(rows, cfg):
    """Absolute QC checks, plus one relative depth check.

    Each flag names the barcode, the observation and the likely cause. A
    threshold with no interpretation just moves the diagnosis downstream.
    """
    depths = [r["raw_reads"] for r in rows if r["raw_reads"]]
    median_depth = sorted(depths)[len(depths) // 2] if depths else None
    flags = []

    for r in rows:
        bc = r["barcode"]
        ret, reads = r["retention_pct"], r["raw_reads"]
        med_len, med_q = r["filtered_median_len"], r["filtered_median_q"]

        if ret is not None and ret < MIN_RETENTION_PCT:
            flags.append((bc, "bad",
                          f"only {ret:.1f}% of reads survived filtering "
                          f"(floor is {MIN_RETENTION_PCT:.0f}%). The length "
                          f"window {cfg['min_length']}-{cfg['max_length']} bp "
                          f"may not match the amplicon."))
        if reads is not None and reads < MIN_READS:
            flags.append((bc, "bad",
                          f"only {reads:,} raw reads (floor is {MIN_READS:,}). "
                          f"Relative abundances from this barcode carry wide "
                          f"uncertainty."))
        elif (reads and median_depth
                and reads < median_depth * LOW_DEPTH_FACTOR):
            flags.append((bc, "warn",
                          f"{reads:,} raw reads is well below the run median "
                          f"of {median_depth:,}. Check the barcoding balance."))
        if med_q is not None and med_q < MIN_MEDIAN_Q:
            flags.append((bc, "bad",
                          f"filtered median quality is Q{med_q:.1f}, under the "
                          f"Q{MIN_MEDIAN_Q:.0f} floor."))
        if med_len is not None and med_len > 0 and not (
                cfg["min_length"] <= med_len <= cfg["max_length"]):
            flags.append((bc, "warn",
                          f"filtered median length {med_len:,.0f} bp sits "
                          f"outside the configured window "
                          f"{cfg['min_length']}-{cfg['max_length']} bp."))
    return flags


def main():
    # snakemake.input is not read here. It exists to make this rule wait for
    # every per-sample job to finish; the numbers come from the benchmark
    # directory and the NanoStat trees, both located via params.
    params = snakemake.params       # noqa: F821
    out_html = Path(snakemake.output.html)    # noqa: F821
    out_csv = Path(snakemake.output.csv)      # noqa: F821
    out_json = Path(snakemake.output.json)    # noqa: F821

    raw_dir = str(params.raw_dir)
    filtered_dir = str(params.filtered_dir)
    cores = int(params.cores or 0)
    cfg = {
        "min_length": int(params.min_length),
        "max_length": int(params.max_length),
        "min_quality": int(params.min_quality),
    }
    version = str(getattr(params, "version", "") or "").strip()
    db = str(getattr(params, "db", "") or "")

    jobs = read_benchmarks(str(params.bench_dir))

    # sample -> rule -> record, for the per-barcode table
    per_stage = {}
    for j in jobs:
        per_stage.setdefault(j["sample"], {})[j["rule"]] = j

    samples = sorted(
        d for d in os.listdir(raw_dir) if os.path.isdir(os.path.join(raw_dir, d))
    ) if os.path.isdir(raw_dir) else []

    rows = build_rows(samples, raw_dir, filtered_dir, per_stage)

    # --- stage rollup ----------------------------------------------------
    seen_stages = [s for s in STAGE_ORDER if any(j["rule"] == s for j in jobs)]
    seen_stages += sorted({j["rule"] for j in jobs} - set(STAGE_ORDER))

    total_raw_reads = sum(r["raw_reads"] or 0 for r in rows)
    rollup = []
    for stage in seen_stages:
        vals = [j for j in jobs if j["rule"] == stage]
        walls = [j["wall"] for j in vals]
        cpus = [j["cpu"] for j in vals if j["cpu"] is not None]
        rsss = [j["max_rss"] for j in vals if j["max_rss"] is not None]
        wall_total = sum(walls)
        rollup.append({
            "stage": stage,
            "jobs": len(vals),
            "wall_total": wall_total,
            "wall_mean": wall_total / len(walls),
            "wall_max": max(walls),
            "cpu_total": sum(cpus) if cpus else None,
            "max_rss_mb": max(rsss) if rsss else None,
            "sec_per_1k_reads": (wall_total / (total_raw_reads / 1000.0)
                                 if total_raw_reads else None),
        })

    total_wall = sum(j["wall"] for j in jobs)
    total_cpu = sum(j["cpu"] for j in jobs if j["cpu"] is not None)
    peak_rss = max((j["max_rss"] for j in jobs if j["max_rss"] is not None),
                   default=None)

    placed = [j for j in jobs if j["start"] is not None]
    if placed:
        span = max(j["end"] for j in placed) - min(j["start"] for j in placed)
        started = datetime.fromtimestamp(min(j["start"] for j in placed))
        finished = datetime.fromtimestamp(max(j["end"] for j in placed))
    else:
        span, started, finished = None, None, None

    # The two headline numbers. Parallelism is what the run achieved;
    # utilisation is what it achieved against what the machine offered.
    parallelism = (total_wall / span) if span else None
    utilisation = (total_cpu / (span * cores) * 100) if (span and cores) else None

    busiest = max(rollup, key=lambda r: r["wall_total"]) if rollup else None
    flags = find_flags(rows, cfg)
    flagged = {bc for bc, _, _ in flags}

    # --- CSV -------------------------------------------------------------
    if rows:
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            for r in rows:
                w.writerow({k: ("" if v is None else v) for k, v in r.items()})
    else:
        out_csv.write_text("barcode\n", encoding="utf-8")

    # --- JSON, for run-to-run comparison ---------------------------------
    # system goes in too: comparing two runs' timings without comparing the
    # machines they ran on is how a hardware difference gets read as a
    # regression.
    out_json.write_text(json.dumps({
        "nano16s_version": version or None,
        "generated": datetime.now().isoformat(timespec="seconds"),
        "system": system_info(),
        "run": {
            "started": started.isoformat(sep=" ", timespec="seconds") if started else None,
            "finished": finished.isoformat(sep=" ", timespec="seconds") if finished else None,
            "span_seconds": span,
            "cores": cores or None,
            "barcodes": len(rows),
            "total_raw_reads": total_raw_reads,
            "job_wall_seconds": total_wall,
            "cpu_seconds": total_cpu,
            "effective_parallelism": parallelism,
            "cpu_utilisation_pct": utilisation,
            "peak_rss_mb": peak_rss,
        },
        "config": cfg,
        "thresholds": {
            "min_retention_pct": MIN_RETENTION_PCT,
            "min_reads": MIN_READS,
            "min_median_q": MIN_MEDIAN_Q,
            "low_depth_factor": LOW_DEPTH_FACTOR,
        },
        "stages": rollup,
        "barcodes": rows,
        "flags": [{"barcode": b, "level": lv, "message": m} for b, lv, m in flags],
    }, indent=2, default=str) + "\n", encoding="utf-8")

    # --- HTML ------------------------------------------------------------
    h = ["<!doctype html><html lang='en'><head><meta charset='utf-8'>",
         "<meta name='viewport' content='width=device-width, initial-scale=1'>",
         "<title>nano16s performance report</title>",
         f"<style>{build_css()}</style></head><body><div class='wrap'>"]

    label = f"nano16s{f' v{esc(version)}' if version else ''}"
    h.append("<header><h1>Performance report</h1>"
             f"<p class='sub'>{label} &middot; {esc(os.path.basename(str(out_html.parent)))} "
             f"&middot; {len(rows)} barcodes</p></header>")

    # headline tiles
    h.append("<section class='stats'>")
    for n, k, hint in [
        (dur(span), "Elapsed", "first job start to last job end"),
        (dur(total_cpu), "CPU time", "summed over every job"),
        (f"{parallelism:.1f}&times;" if parallelism else "-", "Parallelism",
         "jobs in flight, on average"),
        (pct(utilisation, 0) if utilisation else "-", "Core use",
         f"of {cores} cores over the elapsed time" if cores
         else "core count unknown"),
        (f"{peak_rss:,.0f} MB" if peak_rss else "-", "Peak RSS",
         "largest single job"),
    ]:
        h.append(f"<div class='stat'><div class='n'>{n}</div>"
                 f"<div class='k'>{esc(k)}</div><div class='h'>{esc(hint)}</div></div>")
    h.append("</section>")

    # the machine, before any judgement about how well it was used
    sysinfo = system_info()
    ram_gb = (sysinfo["ram_mb"] / 1024.0) if sysinfo["ram_mb"] else None
    core_desc = "-"
    if sysinfo["cpu_cores_physical"] and sysinfo["cpu_threads_logical"]:
        core_desc = f"{sysinfo['cpu_cores_physical']} cores"
        if sysinfo["cpu_threads_logical"] != sysinfo["cpu_cores_physical"]:
            core_desc += f" / {sysinfo['cpu_threads_logical']} threads"
    elif sysinfo["cpu_threads_logical"]:
        core_desc = f"{sysinfo['cpu_threads_logical']} threads"

    h.append("<section class='panel'><h2>System</h2>"
             "<p class='sub'>The machine this run was measured on. Timings are "
             "only comparable between runs on comparable hardware.</p>"
             "<div class='meta'>")
    for k, v in [
        ("CPU", sysinfo["cpu_model"] or "-"),
        ("Cores", core_desc),
        ("Memory", f"{ram_gb:,.1f} GB" if ram_gb else "-"),
        ("Operating system", sysinfo["platform"] or "-"),
        ("Kernel", sysinfo["kernel"] or "-"),
        ("Architecture", sysinfo["arch"] or "-"),
        ("Host", sysinfo["hostname"] or "-"),
        ("Cores given to this run", cores or "-"),
    ]:
        h.append(f"<div><dt>{esc(k)}</dt><dd>{esc(v)}</dd></div>")
    h.append("</div>")
    # Asking for more cores than exist is silently capped by Snakemake, so the
    # run looks fine while quietly being narrower than intended.
    if cores and sysinfo["cpu_threads_logical"] and \
            cores > sysinfo["cpu_threads_logical"]:
        h.append(f"<p class='note'>This run was given {cores} cores but the "
                 f"machine has {sysinfo['cpu_threads_logical']}; Snakemake caps "
                 f"the difference.</p>")
    h.append("</section>")

    # verdict on machine use
    if parallelism and cores:
        # Two different things, easily conflated. Concurrency counts jobs;
        # utilisation counts cores. A single wide job can run alone and still
        # keep the machine busy, so the verdict turns on utilisation and uses
        # concurrency to explain it.
        h.append("<section class='panel'><h2>Machine use</h2>")
        h.append(
            f"<p class='sub'>{dur(total_cpu)} of CPU work finished in "
            f"{dur(span)}, using <strong>{pct(utilisation, 0)}</strong> of the "
            f"{cores} cores available. On average "
            f"<strong>{parallelism:.1f}</strong> job"
            f"{'s were' if parallelism >= 1.5 else ' was'} running at a time.</p>")
        wasted = dur(span * cores * (1 - (utilisation or 0) / 100))
        if utilisation is not None and utilisation < 50:
            h.append(
                f"<p class='sub'>That leaves <strong>{wasted}</strong> of core "
                f"time unused. A job cannot start until its full "
                f"<code>threads:</code> count is free, so large per-rule thread "
                f"counts run fewer barcodes at once and idle the remainder. "
                f"Lowering <code>resources.*.cpus</code> in "
                f"<code>config.yaml</code> trades per-job speed for concurrent "
                f"barcodes, which is normally the better trade when barcodes "
                f"outnumber cores.</p>")
        elif utilisation is not None and utilisation < 80:
            h.append(
                f"<p class='note'>Reasonable, with <strong>{wasted}</strong> of "
                f"core time still idle. Most of that is the tail of the run, "
                f"where too few barcodes remain to fill the machine; lowering "
                f"<code>resources.*.cpus</code> would recover some of it.</p>")
        else:
            h.append("<p class='note'>Little headroom left at this core "
                     "count — the stages are using what the machine has.</p>")
        h.append("</section>")

    # where the time went
    if busiest:
        share = busiest["wall_total"] / max(total_wall, 1e-9) * 100
        h.append("<section class='panel'><h2>Where the time went</h2>"
                 f"<p class='sub'><strong>{esc(busiest['stage'])}</strong> is the "
                 f"bottleneck: {dur(busiest['wall_total'])} of job wall time "
                 f"across {busiest['jobs']} jobs, {share:.0f}% of the total.</p>")
        h.append(stage_bars_svg(rollup))
        h.append("<div class='legend'>"
                 "<span class='lg'><i class='sw' style='background:var(--raw)'></i>"
                 "wall time</span>"
                 "<span class='lg'><i class='sw' style='background:var(--kept)'></i>"
                 "CPU time</span></div>")
        h.append("<div class='scroll'><table class='sortable'><thead><tr>"
                 "<th>Stage</th><th>Jobs</th><th>Wall</th><th>CPU</th>"
                 "<th>Mean/job</th><th>Slowest</th><th>Peak RSS</th>"
                 "<th>s / 1k reads</th></tr></thead><tbody>")
        for r in rollup:
            h.append(
                "<tr>"
                + td(r["stage"])
                + td(r["jobs"])
                + td(r["wall_total"], dur(r["wall_total"]))
                + td(r["cpu_total"], dur(r["cpu_total"]))
                + td(r["wall_mean"], dur(r["wall_mean"]))
                + td(r["wall_max"], dur(r["wall_max"]))
                + td(r["max_rss_mb"], f"{r['max_rss_mb']:,.0f} MB"
                    if r["max_rss_mb"] else "-")
                + td(r["sec_per_1k_reads"], ffix(r["sec_per_1k_reads"], 2))
                + "</tr>")
        h.append("</tbody></table></div>"
                 "<p class='note'>CPU time exceeds wall time for stages that use "
                 "several threads. Stage wall times sum to more than the elapsed "
                 "time because jobs overlap.</p></section>")

    # timeline
    tl, tl_legend = timeline_svg(jobs, seen_stages)
    h.append("<section class='panel'><h2>Timeline</h2>"
             "<p class='sub'>Every job, positioned by when it actually ran. "
             "Gaps are idle capacity.</p>")
    h.append(tl)
    h.append(tl_legend)
    h.append("</section>")

    # QC
    h.append("<section class='panel'>")
    if flags:
        h.append("<h2>Worth checking</h2>")
        for bc, level, msg in flags:
            cls = "warn bad" if level == "bad" else "warn"
            h.append(f"<div class='{cls}'><strong>{esc(bc)}</strong>: {esc(msg)}</div>")
    else:
        h.append("<h2>Worth checking</h2>"
                 f"<p class='ok'>Nothing flagged. Every barcode cleared "
                 f"{MIN_RETENTION_PCT:.0f}% retention, {MIN_READS:,} reads, "
                 f"Q{MIN_MEDIAN_Q:.0f} median quality, and a filtered median "
                 f"length inside the {cfg['min_length']}-{cfg['max_length']} bp "
                 f"window.</p>")
    h.append("</section>")

    # per barcode
    # Timings first. This is the performance report, and on a narrow window a
    # wide table scrolls its rightmost columns out of view — which previously
    # meant the timing columns were the ones you could not see.
    h.append("<section class='panel'><h2>Per barcode</h2>"
             "<p class='sub'>Time spent, then reads and quality before and "
             "after filtering. Click any heading to sort.</p>"
             "<div class='scroll'><table class='sortable'><thead><tr>"
             "<th>Barcode</th><th>Porechop</th><th>Emu</th><th>Total</th>"
             "<th>CPU</th><th>s / 1k</th><th>Peak RSS</th>"
             "<th>Raw reads</th><th>Filtered</th><th>Retained</th>"
             "<th>Raw med Q</th><th>Filt med Q</th><th>Raw &gt;Q15</th>"
             "<th>Filt &gt;Q15</th><th>Raw med len</th><th>Filt med len</th>"
             "</tr></thead><tbody>")
    for r in rows:
        cls = " class='flag'" if r["barcode"] in flagged else ""
        h.append(
            f"<tr{cls}>"
            + td(r["barcode"])
            + td(r["porechop_s"], dur(r["porechop_s"]))
            + td(r["emu_s"], dur(r["emu_s"]))
            + td(r["total_s"], dur(r["total_s"]))
            + td(r["total_cpu_s"], dur(r["total_cpu_s"]))
            + td(r["sec_per_1k_reads"], ffix(r["sec_per_1k_reads"], 2))
            + td(r["max_rss_mb"], f"{r['max_rss_mb']:,.0f} MB"
                if r["max_rss_mb"] else "-")
            + td(r["raw_reads"], fint(r["raw_reads"]))
            + td(r["filtered_reads"], fint(r["filtered_reads"]))
            + td(r["retention_pct"], pct(r["retention_pct"]))
            + td(r["raw_median_q"], ffix(r["raw_median_q"]))
            + td(r["filtered_median_q"], ffix(r["filtered_median_q"]))
            + td(r["raw_pct_Q15"], pct(r["raw_pct_Q15"]))
            + td(r["filtered_pct_Q15"], pct(r["filtered_pct_Q15"]))
            + td(r["raw_median_len"], fint(r["raw_median_len"]))
            + td(r["filtered_median_len"], fint(r["filtered_median_len"]))
            + "</tr>")
    h.append("</tbody></table></div></section>")

    # provenance
    h.append("<section class='panel'><h2>Run</h2><div class='meta'>")
    for k, v in [
        ("Output directory", str(out_html.parent)),
        ("Started", started.strftime("%Y-%m-%d %H:%M:%S") if started else "-"),
        ("Finished", finished.strftime("%Y-%m-%d %H:%M:%S") if finished else "-"),
        ("Cores provided", cores or "-"),
        ("Barcodes", len(rows)),
        ("Total raw reads", f"{total_raw_reads:,}"),
        ("Length window", f"{cfg['min_length']}-{cfg['max_length']} bp"),
        ("Min quality", f"Q{cfg['min_quality']}"),
        ("Emu database", db or "-"),
        ("nano16s version", version or "-"),
    ]:
        h.append(f"<div><dt>{esc(k)}</dt><dd>{esc(v)}</dd></div>")
    h.append("</div></section>")

    h.append(
        "<footer><p>Timings come from Snakemake's <code>benchmark:</code> "
        "records, one per job: wall seconds, CPU seconds and peak resident "
        "memory measured in-process. Elapsed time spans the first job start to "
        "the last job end, so on a resumed run it covers the work still on disk "
        "rather than the most recent invocation. Quality percentages are the "
        "share of reads above that Phred score, straight from NanoStat.</p>"
        "<p>This file is self-contained &mdash; no internet connection is "
        "needed to view it.</p></footer>")
    h.append(f"<script>{JS}</script>")
    h.append("</div></body></html>")

    out_html.write_text("\n".join(h), encoding="utf-8")

    print(f"barcodes      {len(rows)}")
    print(f"elapsed       {dur(span)}")
    print(f"cpu time      {dur(total_cpu)}")
    if parallelism:
        print(f"parallelism   {parallelism:.1f}x of {cores} cores")
    if busiest:
        print(f"bottleneck    {busiest['stage']} ({dur(busiest['wall_total'])})")
    print(f"flags         {len(flags)}")
    print(f"html          {out_html}")


if __name__ == "__main__":
    main()
