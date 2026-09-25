#!/usr/bin/env python3
"""Per-barcode read accounting: where every read went, and what it found.

    read_accounting.py <emu_output_dir> <preprocessing_summary.csv> \
        <read_accounting.tsv> [<per_barcode_taxa.tsv>]

The tables the pipeline already writes answer "what is in this sample?" but not
"does this add up?". The counts table holds one row per taxon, plus a row Emu
leaves unnamed that carries the reads it could not place, so a reader who sums
a column gets a number with no obvious relation to the reads that went in --
and a reader who drops the unnamed row silently loses those reads.

This writes one row per barcode:

    raw -> filtered -> given to the classifier -> classified + unclassified

with the number of species and genera found. `check` is `ok` when
classified + unclassified equals the reads the classifier was given, which is
the arithmetic a reader should not have to do themselves.

With --max-reads, the classifier is given fewer reads than passed the filter;
`subsampled_out` is that difference, so the chain still balances.

The second table answers the other half of the question: not how many species a
barcode found, but which. One row per barcode and species, with its read count
and its share of that barcode. The combined tables hold the same numbers as a
grid of taxa against barcodes, which is what a heatmap wants; this is the shape
a person or a spreadsheet filter wants -- "show me barcode07" or "show me
everything above 1%".
"""
import csv
import os
import sys

# Emu's own labels for reads it did not place. They appear in the tax_id
# column with every taxonomy field empty.
UNPLACED_IDS = {"unmapped", "mapped_filtered", "mapped_unclassified"}


def whole_reads(values, total):
    """Round fractional per-species counts so they still sum to `total`.

    Emu estimates counts, so they are fractional. Rounding each one on its own
    leaves the column a read or two away from the barcode's total, which is
    exactly the reconciliation this table exists to provide. Largest remainder:
    floor everything, then give the remaining reads to the largest fractions.
    """
    floors = [int(v) for v in values]
    short = total - sum(floors)
    if short <= 0:
        return floors
    order = sorted(range(len(values)), key=lambda i: -(values[i] - floors[i]))
    for i in order[:short]:
        floors[i] += 1
    return floors


def count_file(path):
    """One barcode's Emu table -> (classified, unclassified, taxa).

    taxa is [(species, genus, reads)], one entry per species with reads.
    """
    classified = unclassified = 0.0
    taxa = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                n = float(row.get("estimated counts") or 0)
            except ValueError:
                continue
            if n <= 0:
                continue
            name = (row.get("species") or "").strip()
            tax_id = (row.get("tax_id") or "").strip()
            if tax_id in UNPLACED_IDS or not name:
                unclassified += n
                continue
            classified += n
            genus = (row.get("genus") or "").strip()
            key = (name, genus)
            taxa[key] = taxa.get(key, 0.0) + n
    rows = [(sp, gen, n) for (sp, gen), n in taxa.items()]
    rows.sort(key=lambda r: -r[2])
    return classified, unclassified, rows


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    emu_dir, summary_csv, out_tsv = sys.argv[1:4]
    taxa_tsv = sys.argv[4] if len(sys.argv) == 5 else None

    pre = {}
    if os.path.exists(summary_csv):
        with open(summary_csv) as fh:
            for r in csv.DictReader(fh):
                pre[r["barcode"]] = r

    rows = []
    taxa_rows = []
    for barcode in sorted(pre) or sorted(os.listdir(emu_dir)):
        f = os.path.join(emu_dir, barcode, f"{barcode}_rel-abundance.tsv")
        if not os.path.exists(f) or os.path.getsize(f) == 0:
            # A barcode with no reads is reported with zeros rather than
            # dropped: a missing row reads as an oversight, a zero row as a
            # fact about the run.
            classified = unclassified = 0.0
            found = []
        else:
            classified, unclassified, found = count_file(f)
        n_species = len(found)
        n_genera = len({g for _, g, _ in found if g})
        whole = whole_reads([r[2] for r in found], int(round(classified)))
        for (species, genus, reads), n in zip(found, whole):
            taxa_rows.append({
                "barcode": barcode,
                "species": species,
                "genus": genus,
                "reads": n,
                "pct_of_classified": f"{reads / classified:.4%}" if classified else "-",
            })

        def num(key):
            try:
                return int(float(pre.get(barcode, {}).get(key, 0) or 0))
            except ValueError:
                return 0

        raw, filtered = num("raw_reads"), num("filtered_reads")
        to_classifier = int(round(classified + unclassified))
        rows.append({
            "barcode": barcode,
            "raw_reads": raw,
            "removed_by_filter": raw - filtered,
            "filtered_reads": filtered,
            "subsampled_out": max(0, filtered - to_classifier),
            "reads_to_classifier": to_classifier,
            "reads_classified": int(round(classified)),
            "reads_unclassified": int(round(unclassified)),
            "unclassified_pct": (f"{unclassified / to_classifier:.2%}"
                                 if to_classifier else "-"),
            "species_found": n_species,
            "genera_found": n_genera,
            "check": "ok" if to_classifier == int(round(classified)) + int(round(unclassified))
                     else "MISMATCH",
        })

    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["barcode"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if taxa_tsv:
        with open(taxa_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
                "barcode", "species", "genus", "reads", "pct_of_classified"])
            w.writeheader()
            w.writerows(taxa_rows)

    total = sum(r["reads_to_classifier"] for r in rows)
    placed = sum(r["reads_classified"] for r in rows)
    print(f"read accounting: {len(rows)} barcodes, {total:,} reads classified or not, "
          f"{placed:,} placed ({placed / total:.1%})" if total else
          f"read accounting: {len(rows)} barcodes, no classified reads")


if __name__ == "__main__":
    main()
