#!/usr/bin/env python3
"""Subsample a full nanopore run into the small demo dataset shipped with nano16s.

Reads are drawn evenly across the source's many per-barcode FASTQ files rather
than taken from the first few, so the sample spans the whole sequencing run
instead of just its opening minutes.

Output keeps the two-file-per-barcode shape of real MinKNOW output, so the
merge step is genuinely exercised by `nano16s test` rather than bypassed.

Usage:
    python tools/make_demo.py --source /path/to/fastq_pass --dest test/demo/fastq_pass
"""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

READS_PER_BARCODE = 2000
FILES_PER_BARCODE = 2


def subsample_barcode(src: Path, dest: Path, n_reads: int, n_files: int) -> int:
    """Write n_reads from src, spread across its files, into n_files outputs."""
    sources = sorted(src.glob("*.fastq.gz"))
    if not sources:
        return 0

    # Even stride across the source files, so we sample the whole run.
    step = max(1, len(sources) // 12)
    picked = sources[::step] or sources

    records: list[str] = []
    per_file = max(1, n_reads // len(picked) + 1)
    for path in picked:
        taken = 0
        with gzip.open(path, "rt") as fh:
            while taken < per_file:
                block = [fh.readline() for _ in range(4)]
                if not block[0]:
                    break
                records.append("".join(block))
                taken += 1
        if len(records) >= n_reads:
            break

    records = records[:n_reads]
    if not records:
        return 0

    dest.mkdir(parents=True, exist_ok=True)
    chunk = (len(records) + n_files - 1) // n_files
    for i in range(n_files):
        part = records[i * chunk:(i + 1) * chunk]
        if not part:
            continue
        out = dest / f"{dest.name}_demo_{i}.fastq.gz"
        with gzip.open(out, "wt", compresslevel=9) as fh:
            fh.writelines(part)
    return len(records)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source", type=Path, required=True,
                    help="fastq_pass directory of the full run")
    ap.add_argument("--dest", type=Path, required=True,
                    help="where to write the demo dataset")
    ap.add_argument("--reads", type=int, default=READS_PER_BARCODE)
    ap.add_argument("--files", type=int, default=FILES_PER_BARCODE)
    args = ap.parse_args()

    if args.dest.exists():
        shutil.rmtree(args.dest)
    args.dest.mkdir(parents=True)

    total = 0
    for bc in sorted(args.source.glob("barcode*")):
        if not bc.is_dir():
            continue
        n = subsample_barcode(bc, args.dest / bc.name, args.reads, args.files)
        total += n
        print(f"  {bc.name:<12} {n:>6,} reads")

    size = sum(f.stat().st_size for f in args.dest.rglob("*.fastq.gz"))
    print(f"\n  total: {total:,} reads, {size / 1e6:.1f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
