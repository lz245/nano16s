#!/usr/bin/env python3
"""Give Emu's unnamed row a name, in a combined table.

    label_unclassified.py <rank> <table.tsv> [<table.tsv> ...]

`emu combine-outputs` writes one row whose taxonomy columns are all empty. It
carries the reads Emu could not place, and at a glance it looks like a blank
line. Summing a column therefore includes reads with no visible label, and
anyone who filters out unnamed rows -- a reasonable thing to do -- drops them
without noticing.

The row keeps its counts; it gains the label `Unclassified` in the rank column,
so a reader sees what it is and a script can keep or drop it deliberately.
"""
import csv
import os
import sys

LABEL = "Unclassified"
# Columns that hold taxonomy rather than a barcode's counts.
TAXONOMY = {"tax_id", "taxid", "species", "genus", "family", "order", "class",
            "phylum", "superkingdom", "domain", "kingdom", "subspecies",
            "clade", "species subgroup", "species group"}


def label(path, rank):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return 0
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows or len(rows) < 2:
        return 0
    header = rows[0]
    tax_cols = [i for i, c in enumerate(header) if c.strip().lower() in TAXONOMY]
    if rank in header:
        rank_col = header.index(rank)
    elif tax_cols:
        rank_col = tax_cols[0]
    else:
        return 0

    changed = 0
    for row in rows[1:]:
        if len(row) <= rank_col:
            continue
        if all(not (row[i].strip() if i < len(row) else "") for i in tax_cols):
            row[rank_col] = LABEL
            changed += 1
    if changed:
        with open(path, "w", newline="") as fh:
            csv.writer(fh, delimiter="\t").writerows(rows)
    return changed


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    rank, paths = sys.argv[1], sys.argv[2:]
    for p in paths:
        n = label(p, rank)
        if n:
            print(f"labelled {n} unnamed row(s) as {LABEL} in {os.path.basename(p)}")


if __name__ == "__main__":
    main()
