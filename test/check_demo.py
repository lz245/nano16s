#!/usr/bin/env python3
"""Check that a demo run produced sensible output.

This verifies the installation, not the science. It asserts structure and
plausibility — every barcode processed, reads survived, taxa were assigned,
abundances sum to one — rather than specific species, so a database update
cannot turn it red.

Usage:
    python test/check_demo.py <results_dir>

Exits 0 if everything looks right, 1 otherwise.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

EXPECTED_BARCODES = 6
RANKS = ("species", "genus", "phylum")

GREEN, RED, YELLOW, RESET = "\033[32m", "\033[31m", "\033[33m", "\033[0m"
if not sys.stdout.isatty():
    GREEN = RED = YELLOW = RESET = ""

problems: list[str] = []
notes: list[str] = []


def ok(msg: str) -> None:
    print(f"  {GREEN}ok{RESET}    {msg}")


def fail(msg: str) -> None:
    print(f"  {RED}FAIL{RESET}  {msg}")
    problems.append(msg)


def warn(msg: str) -> None:
    print(f"  {YELLOW}note{RESET}  {msg}")
    notes.append(msg)


def check_files(root: Path) -> None:
    expected = [root / "preprocessing_summary.csv", root / "nano16s_report.html"]
    for rank in RANKS:
        expected.append(root / "07_emu_combined" / f"emu-combined-{rank}.tsv")
        expected.append(root / "07_emu_combined" / f"emu-combined-{rank}-counts.tsv")

    for path in expected:
        rel = path.relative_to(root)
        if not path.exists():
            fail(f"missing output: {rel}")
        elif path.stat().st_size == 0:
            fail(f"empty output: {rel}")
        else:
            ok(f"{rel} ({path.stat().st_size:,} bytes)")


def check_preprocessing(root: Path) -> None:
    path = root / "preprocessing_summary.csv"
    if not path.exists():
        return
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))

    if len(rows) != EXPECTED_BARCODES:
        fail(f"expected {EXPECTED_BARCODES} barcodes in the summary, found {len(rows)}")
    else:
        ok(f"all {EXPECTED_BARCODES} barcodes present")

    total_raw = total_kept = 0
    for r in rows:
        bc = r.get("barcode", "?")
        try:
            raw = int(float(r.get("raw_reads") or 0))
            kept = int(float(r.get("filtered_reads") or 0))
        except ValueError:
            fail(f"{bc}: unparseable read counts")
            continue
        total_raw += raw
        total_kept += kept
        if raw == 0:
            fail(f"{bc}: no raw reads — the merge step produced nothing")
        elif kept == 0:
            fail(f"{bc}: every read was filtered out")
        elif kept / raw < 0.10:
            warn(f"{bc}: retained only {100 * kept / raw:.0f}% of reads")

    if total_kept:
        ok(f"{total_kept:,} of {total_raw:,} reads survived filtering "
           f"({100 * total_kept / total_raw:.0f}%)")


def check_abundance(root: Path) -> None:
    path = root / "07_emu_combined" / "emu-combined-species.tsv"
    if not path.exists() or path.stat().st_size == 0:
        return

    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header or header[0].startswith("#"):
            fail("species table has no data — no barcode produced classifications")
            return
        try:
            label_i = header.index("species")
        except ValueError:
            label_i = 1

        # Emu writes the full lineage before the sample columns, so anything
        # named for a taxonomic rank is taxonomy rather than a barcode.
        lineage = {"tax_id", "taxid", "species", "genus", "family", "order",
                   "class", "phylum", "superkingdom", "domain", "kingdom",
                   "subspecies", "clade", "lineage", "estimated counts", ""}
        rows = [r for r in reader if len(r) > label_i]
        sample_idx = []
        for i, h in enumerate(header):
            if i == label_i or h.strip().lower() in lineage:
                continue
            for r in rows:
                if i < len(r) and r[i].strip():
                    try:
                        float(r[i])
                        sample_idx.append(i)
                    except ValueError:
                        pass
                    break

        sums = {i: 0.0 for i in sample_idx}
        n_taxa = 0
        for row in rows:
            n_taxa += 1
            for i in sample_idx:
                if i < len(row) and row[i].strip():
                    try:
                        sums[i] += float(row[i])
                    except ValueError:
                        pass

    if n_taxa == 0:
        fail("species table has a header but no rows")
        return
    ok(f"{n_taxa:,} species assigned across {len(sample_idx)} barcodes")

    # Emu reports relative abundance, so each barcode's column should total ~1.
    for i, total in sums.items():
        name = header[i]
        if total == 0:
            fail(f"{name}: all abundances are zero")
        elif not (0.95 <= total <= 1.05):
            warn(f"{name}: abundances sum to {total:.3f}, expected ~1.0")
    if all(0.95 <= v <= 1.05 for v in sums.values() if v):
        ok("relative abundances sum to 1.0 in every barcode")


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1])
    if not root.is_dir():
        print(f"ERROR: not a directory: {root}", file=sys.stderr)
        return 1

    print(f"Checking {root}\n")
    check_files(root)
    print()
    check_preprocessing(root)
    print()
    check_abundance(root)
    print()

    if problems:
        print(f"{RED}{len(problems)} problem(s) found.{RESET}")
        print("The install is not working correctly. The messages above say what failed.")
        return 1

    if notes:
        print(f"{GREEN}Install verified.{RESET} "
              f"{len(notes)} note(s) above are about the demo data, not the software.")
    else:
        print(f"{GREEN}Install verified — everything works.{RESET}")
    print(f"\nOpen the report:  {root / 'nano16s_report.html'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
