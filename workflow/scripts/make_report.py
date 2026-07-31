#!/usr/bin/env python3
"""Build a self-contained HTML summary of a nano16s run.

No external assets, no plotting library — charts are inline SVG generated
here. The file opens by double-click and can be emailed as-is.

Driven by Snakemake (`snakemake.input`, `snakemake.output`, `snakemake.params`).
"""

from __future__ import annotations

import csv
import html
import json
import os
import subprocess
from datetime import datetime
from pathlib import Path

# Categorical palette, validated for colorblind separation and contrast in
# both light and dark modes. Slots are assigned in fixed order and never
# cycled: beyond slot 7 everything folds into "Other".
PALETTE_LIGHT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                 "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
PALETTE_DARK = ["#3987e5", "#d95926", "#199e70", "#c98500",
                "#d55181", "#008300", "#9085e9", "#e66767"]
OTHER_LIGHT, OTHER_DARK = "#9a9a94", "#7b7b74"

TOP_N = 7  # + "Other" = 8 slots, the documented cap


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------

def read_summary(path):
    """preprocessing_summary.csv -> list of per-barcode dicts."""
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            def num(k):
                try:
                    return int(float(r.get(k) or 0))
                except (TypeError, ValueError):
                    return 0
            rows.append({
                "barcode": r.get("barcode", "?"),
                "raw": num("raw_reads"),
                "filtered": num("filtered_reads"),
                "raw_len": r.get("raw_median_length") or "-",
                "filt_len": r.get("filtered_median_length") or "-",
                "raw_q": r.get("raw_median_quality") or "-",
                "filt_q": r.get("filtered_median_quality") or "-",
            })
    return rows


def read_abundance(path, rank):
    """Emu combined table -> (sample_names, {taxon: {sample: fraction}}).

    Emu writes tax_id, the rank column, then one column per sample. Sample
    column names carry Emu's own suffixes, which are stripped for display.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return [], {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return [], {}
        if header and header[0].startswith("#"):
            return [], {}

        # The label column is the one named for the rank; fall back to the
        # first non-tax_id column.
        try:
            label_i = header.index(rank)
        except ValueError:
            label_i = 1 if len(header) > 1 else 0

        sample_idx = [i for i, h in enumerate(header)
                      if i != label_i and h.lower() not in ("tax_id", "taxid", "")]
        samples = [_clean_sample(header[i]) for i in sample_idx]

        table: dict[str, dict[str, float]] = {}
        for row in reader:
            if not row or len(row) <= label_i:
                continue
            taxon = (row[label_i] or "").strip() or "Unassigned"
            per = table.setdefault(taxon, {s: 0.0 for s in samples})
            for s, i in zip(samples, sample_idx):
                if i < len(row):
                    try:
                        per[s] += float(row[i])
                    except (TypeError, ValueError):
                        pass
    return samples, table


def _clean_sample(name: str) -> str:
    for suffix in ("_rel-abundance", "_filtered", ".tsv", "_counts"):
        name = name.replace(suffix, "")
    return name.strip()


def tool_versions():
    """Best-effort version capture, for the Methods paragraph."""
    wanted = {
        "emu": ["emu", "--version"],
        "minimap2": ["minimap2", "--version"],
        "chopper": ["chopper", "--version"],
        "porechop_abi": ["porechop_abi", "--version"],
        "NanoStat": ["NanoStat", "--version"],
        "snakemake": ["snakemake", "--version"],
    }
    out = {}
    for name, cmd in wanted.items():
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
            line = (r.stdout or r.stderr).strip().splitlines()
            out[name] = line[0].split()[-1] if line else "unknown"
        except Exception:
            out[name] = "unknown"
    return out


# ---------------------------------------------------------------------------
# charts
# ---------------------------------------------------------------------------

def esc(s) -> str:
    return html.escape(str(s), quote=True)


def funnel_svg(rows) -> str:
    """Horizontal bars: reads retained after filtering, per barcode.

    Bars are read for magnitude and compared across barcodes, so a shared
    x-scale is used rather than per-row normalisation.
    """
    if not rows:
        return "<p class='empty'>No preprocessing summary available.</p>"

    row_h, gap, pad_l, pad_r, pad_t = 26, 6, 96, 76, 26
    width = 720
    plot_w = width - pad_l - pad_r
    height = pad_t + len(rows) * (row_h + gap) + 8
    peak = max((r["raw"] for r in rows), default=1) or 1

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="Reads per barcode before and after filtering">'
    ]
    # x gridlines
    for frac in (0.25, 0.5, 0.75, 1.0):
        x = pad_l + plot_w * frac
        parts.append(f'<line class="grid" x1="{x:.1f}" y1="{pad_t - 8}" '
                     f'x2="{x:.1f}" y2="{height - 8}"/>')
        parts.append(f'<text class="tick" x="{x:.1f}" y="{pad_t - 12}" '
                     f'text-anchor="middle">{int(peak * frac):,}</text>')

    for i, r in enumerate(rows):
        y = pad_t + i * (row_h + gap)
        raw_w = plot_w * (r["raw"] / peak) if peak else 0
        filt_w = plot_w * (r["filtered"] / peak) if peak else 0
        pct = (100.0 * r["filtered"] / r["raw"]) if r["raw"] else 0.0

        parts.append(f'<text class="rowlab" x="{pad_l - 10}" y="{y + 17}" '
                     f'text-anchor="end">{esc(r["barcode"])}</text>')
        # raw = recessive track, filtered = emphasised fill on the same baseline
        parts.append(f'<rect class="raw" x="{pad_l}" y="{y + 4}" '
                     f'width="{max(raw_w, 1):.1f}" height="{row_h - 8}" rx="4"/>')
        parts.append(f'<rect class="kept" x="{pad_l}" y="{y + 4}" '
                     f'width="{max(filt_w, 1):.1f}" height="{row_h - 8}" rx="4"/>')
        parts.append(f'<title>{esc(r["barcode"])}: {r["raw"]:,} raw, '
                     f'{r["filtered"]:,} kept ({pct:.0f}%)</title>')
        parts.append(f'<text class="rowval" x="{pad_l + plot_w + 8}" y="{y + 17}">'
                     f'{pct:.0f}%</text>')
    parts.append("</svg>")
    return "".join(parts)


def composition_svg(samples, table, rank) -> tuple[str, str, list[tuple[str, int]]]:
    """Stacked horizontal bars of relative abundance, one row per barcode.

    Returns (svg, legend_html, ordered_taxa). Taxa are ranked by mean
    abundance across barcodes so a colour always means the same organism in
    every row.
    """
    if not samples or not table:
        return ("<p class='empty'>No classification results available.</p>", "", [])

    totals = {t: sum(v.values()) for t, v in table.items()}
    ranked = sorted(totals, key=totals.get, reverse=True)
    top = ranked[:TOP_N]
    slots = [(t, i) for i, t in enumerate(top)]
    has_other = len(ranked) > TOP_N

    row_h, gap, pad_l, pad_r, pad_t = 28, 6, 96, 56, 24
    width = 720
    plot_w = width - pad_l - pad_r
    height = pad_t + len(samples) * (row_h + gap) + 8

    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" '
        f'aria-label="Relative abundance by {esc(rank)} for each barcode">'
    ]
    for i, s in enumerate(samples):
        y = pad_t + i * (row_h + gap)
        col_total = sum(table[t].get(s, 0.0) for t in table) or 1.0
        parts.append(f'<text class="rowlab" x="{pad_l - 10}" y="{y + 19}" '
                     f'text-anchor="end">{esc(s)}</text>')

        x = pad_l
        for taxon, slot in slots:
            frac = table[taxon].get(s, 0.0) / col_total
            w = plot_w * frac
            if w <= 0.4:
                continue
            parts.append(
                f'<rect class="seg s{slot}" x="{x:.2f}" y="{y + 4}" '
                f'width="{max(w - 2, 0.6):.2f}" height="{row_h - 8}" rx="2">'
                f'<title>{esc(s)} — {esc(taxon)}: {frac * 100:.1f}%</title></rect>'
            )
            # direct label only where the segment can hold it
            if w > 54:
                parts.append(f'<text class="seglab" x="{x + w / 2 - 1:.1f}" '
                             f'y="{y + 19}" text-anchor="middle">'
                             f'{frac * 100:.0f}%</text>')
            x += w
        if has_other:
            other = plot_w + pad_l - x
            if other > 0.4:
                oth_frac = other / plot_w
                parts.append(
                    f'<rect class="seg other" x="{x:.2f}" y="{y + 4}" '
                    f'width="{max(other - 2, 0.6):.2f}" height="{row_h - 8}" rx="2">'
                    f'<title>{esc(s)} — Other taxa: {oth_frac * 100:.1f}%</title></rect>'
                )
    parts.append("</svg>")

    legend = ['<div class="legend">']
    for taxon, slot in slots:
        legend.append(f'<span class="lg"><i class="sw s{slot}"></i>'
                      f'<em>{esc(taxon)}</em></span>')
    if has_other:
        legend.append(f'<span class="lg"><i class="sw other"></i>'
                      f'Other ({len(ranked) - TOP_N} taxa)</span>')
    legend.append("</div>")
    return "".join(parts), "".join(legend), [(t, 0) for t in ranked]


# ---------------------------------------------------------------------------
# page
# ---------------------------------------------------------------------------

# Filled with str.replace, not %-formatting: the template contains literal
# percent signs (width: 100%) that %-formatting would read as conversion
# specifiers.
CSS = """
:root {
  --surface: #fcfcfb; --panel: #ffffff; --line: #e3e2dd;
  --ink: #0b0b0b; --ink-2: #52514e; --ink-3: #78776f;
  --raw: #dedcd4; --kept: #2a78d6; --other: @OTHER_L@;
  @SLOTS_L@
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --surface: #161615; --panel: #1e1e1c; --line: #33322e;
    --ink: #f4f3ef; --ink-2: #c3c2b7; --ink-3: #918f85;
    --raw: #33322e; --kept: #3987e5; --other: @OTHER_D@;
    @SLOTS_D@
  }
}
:root[data-theme="dark"] {
  --surface: #161615; --panel: #1e1e1c; --line: #33322e;
  --ink: #f4f3ef; --ink-2: #c3c2b7; --ink-3: #918f85;
  --raw: #33322e; --kept: #3987e5; --other: @OTHER_D@;
  @SLOTS_D@
}
* { box-sizing: border-box; }
body {
  margin: 0; background: var(--surface); color: var(--ink);
  font: 15px/1.6 -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  -webkit-font-smoothing: antialiased;
}
.wrap { max-width: 60rem; margin: 0 auto; padding: 3rem 1.5rem 5rem;
        display: flex; flex-direction: column; gap: 2.25rem; }
h1 { font-size: 1.9rem; margin: 0; letter-spacing: -0.02em; }
h2 { font-size: 1.15rem; margin: 0 0 .25rem; letter-spacing: -0.01em; }
.sub { color: var(--ink-2); margin: .4rem 0 0; }
.panel { background: var(--panel); border: 1px solid var(--line);
         border-radius: 6px; padding: 1.4rem 1.5rem; }
.meta { display: grid; grid-template-columns: repeat(auto-fit, minmax(11rem, 1fr));
        gap: .9rem 1.5rem; }
.meta div { display: flex; flex-direction: column; gap: .15rem; }
.meta dt { font-size: .7rem; letter-spacing: .09em; text-transform: uppercase;
           color: var(--ink-3); }
.meta dd { margin: 0; font-variant-numeric: tabular-nums; word-break: break-word; }
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(8rem, 1fr));
         gap: 1rem; }
.stat { background: var(--panel); border: 1px solid var(--line); border-radius: 6px;
        padding: .9rem 1rem; }
.stat .n { font-size: 1.6rem; font-weight: 600; font-variant-numeric: tabular-nums;
           letter-spacing: -0.02em; }
.stat .k { font-size: .72rem; letter-spacing: .08em; text-transform: uppercase;
           color: var(--ink-3); margin-top: .15rem; }
svg { display: block; width: 100%; height: auto; overflow: visible; }
.grid { stroke: var(--line); stroke-width: 1; }
.tick, .rowlab, .rowval { font-size: 11px; fill: var(--ink-3); }
.rowlab { font-size: 12px; fill: var(--ink-2); }
.rowval { font-size: 12px; fill: var(--ink-2); font-variant-numeric: tabular-nums; }
.seglab { font-size: 10.5px; fill: #ffffff; font-variant-numeric: tabular-nums;
          paint-order: stroke; stroke: rgba(0,0,0,.28); stroke-width: 2.5px; }
rect.raw { fill: var(--raw); }
rect.kept { fill: var(--kept); }
rect.other, i.sw.other { background: var(--other); fill: var(--other); }
@SLOTCSS@
.legend { display: flex; flex-wrap: wrap; gap: .5rem 1.1rem; margin-top: 1rem;
          font-size: .84rem; color: var(--ink-2); }
.lg { display: inline-flex; align-items: center; gap: .4rem; }
.sw { width: .8rem; height: .8rem; border-radius: 2px; display: inline-block; flex: none; }
.legend em { font-style: italic; }
.scroll { overflow-x: auto; }
table { border-collapse: collapse; width: 100%; font-size: .87rem; margin-top: .4rem; }
th, td { text-align: right; padding: .45rem .6rem; border-bottom: 1px solid var(--line);
         white-space: nowrap; font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
thead th { font-size: .68rem; letter-spacing: .08em; text-transform: uppercase;
           color: var(--ink-3); font-weight: 600; }
tbody tr:last-child td { border-bottom: none; }
.warn { border-left: 3px solid #eda100; padding-left: .9rem; margin: .5rem 0; }
.warn.bad { border-left-color: #e34948; }
.ok { color: var(--ink-2); }
.empty { color: var(--ink-3); font-style: italic; }
details { margin-top: 1rem; }
summary { cursor: pointer; color: var(--ink-2); font-size: .88rem; }
code { font-family: ui-monospace, "SF Mono", Menlo, monospace; font-size: .86em;
       background: var(--surface); border: 1px solid var(--line);
       border-radius: 3px; padding: .05em .35em; }
footer { color: var(--ink-3); font-size: .8rem; border-top: 1px solid var(--line);
         padding-top: 1rem; }
"""


def build_css() -> str:
    slots_l = " ".join(f"--s{i}: {c};" for i, c in enumerate(PALETTE_LIGHT))
    slots_d = " ".join(f"--s{i}: {c};" for i, c in enumerate(PALETTE_DARK))
    slotcss = "\n".join(
        f"rect.seg.s{i}, i.sw.s{i} {{ fill: var(--s{i}); background: var(--s{i}); }}"
        for i in range(len(PALETTE_LIGHT))
    )
    return (CSS
            .replace("@SLOTS_L@", slots_l)
            .replace("@SLOTS_D@", slots_d)
            .replace("@SLOTCSS@", slotcss)
            .replace("@OTHER_L@", OTHER_LIGHT)
            .replace("@OTHER_D@", OTHER_DARK))


def main():
    inp = snakemake.input           # noqa: F821
    params = snakemake.params       # noqa: F821
    out = Path(snakemake.output[0])  # noqa: F821

    rows = read_summary(inp.summary)
    samples, species = read_abundance(inp.species, "species")
    _, genus = read_abundance(inp.genus, "genus")
    versions = tool_versions()

    total_raw = sum(r["raw"] for r in rows)
    total_kept = sum(r["filtered"] for r in rows)
    retention = (100.0 * total_kept / total_raw) if total_raw else 0.0

    db_path = Path(str(params.db))
    manifest = db_path / "MANIFEST.json"
    db_desc = db_path.name
    if manifest.exists():
        try:
            m = json.loads(manifest.read_text())
            db_desc = (f"{db_path.parent.name} — {m.get('sequences', '?'):,} sequences, "
                       f"{m.get('taxa', '?'):,} taxa, built {m.get('built', '?')}")
        except Exception:
            pass

    # --- warnings -----------------------------------------------------------
    warnings = []
    for r in rows:
        if r["raw"] == 0:
            warnings.append(("bad", f"<strong>{esc(r['barcode'])}</strong> had no reads "
                                    f"at all. Check the barcode was used in this run."))
        elif r["filtered"] == 0:
            warnings.append(("bad", f"<strong>{esc(r['barcode'])}</strong> lost every read "
                                    f"during filtering ({r['raw']:,} raw). Its reads are "
                                    f"probably outside the "
                                    f"{params.min_length}-{params.max_length} bp window."))
        else:
            pct = 100.0 * r["filtered"] / r["raw"]
            if pct < 25:
                warnings.append(("warn", f"<strong>{esc(r['barcode'])}</strong> retained only "
                                         f"{pct:.0f}% of its reads. Low retention usually means "
                                         f"short fragments or low quality, not a pipeline fault."))
        if 0 < r["filtered"] < 1000:
            warnings.append(("warn", f"<strong>{esc(r['barcode'])}</strong> has only "
                                     f"{r['filtered']:,} reads after filtering. Diversity "
                                     f"estimates from this few reads are unreliable."))

    # --- assemble -----------------------------------------------------------
    comp_svg, legend, ranked = composition_svg(samples, species, "species")
    gcomp_svg, glegend, _ = composition_svg(samples, genus, "genus")

    top_rows = []
    for taxon in [t for t, _ in ranked][:15]:
        vals = species.get(taxon, {})
        tot = sum(vals.values())
        mean = (100.0 * tot / len(samples)) if samples else 0.0
        top_rows.append(
            "<tr><td><em>%s</em></td><td>%.2f</td>%s</tr>" % (
                esc(taxon), mean,
                "".join(f"<td>{100.0 * vals.get(s, 0.0):.1f}</td>" for s in samples))
        )

    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width, initial-scale=1'>",
        "<title>nano16s report</title>",
        f"<style>{build_css()}</style></head><body><div class='wrap'>",

        "<header>",
        "<h1>nano16s report</h1>",
        f"<p class='sub'>{len(rows)} barcodes &middot; "
        f"{total_raw:,} reads in &middot; {total_kept:,} kept "
        f"({retention:.0f}%) &middot; "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M')}</p>",
        "</header>",

        "<section class='stats'>",
        f"<div class='stat'><div class='n'>{len(rows)}</div><div class='k'>Barcodes</div></div>",
        f"<div class='stat'><div class='n'>{total_kept:,}</div><div class='k'>Reads classified</div></div>",
        f"<div class='stat'><div class='n'>{retention:.0f}%</div><div class='k'>Retained</div></div>",
        f"<div class='stat'><div class='n'>{len(species)}</div><div class='k'>Species detected</div></div>",
        "</section>",
    ]

    if warnings:
        parts.append("<section class='panel'><h2>Worth checking</h2>")
        for kind, msg in warnings:
            parts.append(f"<p class='warn {'bad' if kind == 'bad' else ''}'>{msg}</p>")
        parts.append("</section>")
    else:
        parts.append("<section class='panel'><h2>Worth checking</h2>"
                     "<p class='ok'>Nothing unusual — every barcode produced reads and "
                     "retained a reasonable fraction of them.</p></section>")

    parts += [
        "<section class='panel'><h2>Reads through filtering</h2>",
        "<p class='sub'>Pale bar is what came off the sequencer; solid bar is what "
        "survived adapter trimming and the length/quality filter. Percentage is retention.</p>",
        funnel_svg(rows),
        "</section>",

        "<section class='panel'><h2>Composition by species</h2>",
        f"<p class='sub'>Top {TOP_N} species by mean relative abundance; everything else "
        f"is grouped as Other. Hover any segment for its exact value.</p>",
        comp_svg, legend,
        "</section>",

        "<section class='panel'><h2>Composition by genus</h2>",
        "<p class='sub'>The same reads summarised one rank up. Genus-level calls are "
        "more reliable than species-level ones for nanopore 16S.</p>",
        gcomp_svg, glegend,
        "</section>",

        "<section class='panel'><h2>Top species (% relative abundance)</h2>",
        "<div class='scroll'><table><thead><tr><th>Species</th><th>Mean</th>",
        "".join(f"<th>{esc(s)}</th>" for s in samples),
        "</tr></thead><tbody>", "".join(top_rows), "</tbody></table></div>",
        "<details><summary>Where are the full tables?</summary>"
        "<p class='sub'>Complete tables for every rank are in "
        "<code>07_emu_combined/</code> — both relative abundance and read counts, "
        "as tab-separated files you can open in R, Python, or Excel.</p></details>",
        "</section>",

        "<section class='panel'><h2>Methods</h2>",
        "<div class='meta'>",
        f"<div><dt>Database</dt><dd>{esc(db_desc)}</dd></div>",
        f"<div><dt>Read filter</dt><dd>{params.min_length}–{params.max_length} bp, "
        f"Q &ge; {params.min_quality}</dd></div>",
        "".join(f"<div><dt>{esc(k)}</dt><dd>{esc(v)}</dd></div>"
                for k, v in versions.items()),
        "</div>",
        "<details><summary>Paragraph for a manuscript</summary>",
        f"<p class='sub'>Reads were merged per barcode, adapter-trimmed with "
        f"Porechop_ABI v{esc(versions.get('porechop_abi'))}, and filtered to "
        f"{params.min_length}–{params.max_length} bp at Q&ge;{params.min_quality} with "
        f"Chopper v{esc(versions.get('chopper'))}. Read statistics before and after "
        f"filtering were computed with NanoStat v{esc(versions.get('NanoStat'))}. "
        f"Taxonomic abundances were estimated with Emu "
        f"v{esc(versions.get('emu'))} against {esc(db_desc)}. "
        f"The workflow was executed with Snakemake "
        f"v{esc(versions.get('snakemake'))}.</p></details>",
        "</section>",

        "<footer>Generated by nano16s. This file is self-contained — "
        "no internet connection is needed to view it.</footer>",
        "</div></body></html>",
    ]

    out.write_text("".join(parts), encoding="utf-8")
    print(f"  wrote {out}")


main()
