#!/usr/bin/env python3
"""One line per sequencing read: what that read was called, and how sure Emu is.

    per_read_calls.py <read-assignment-distributions.tsv> <db_dir> \\
        <classifier_input.fastq.gz> <barcode> <out.tsv.gz>

The abundance tables say what is in a barcode. They cannot say which read
supported which call, so a read cannot be traced back, pulled out for a second
opinion, or counted by hand.

Emu is probabilistic: rather than labelling a read, it spreads that read across
the references it matched and estimates abundances from the whole distribution.
Its `--keep-read-assignments` file is that distribution, one row per read and
one column per taxid. This turns it into the per-read view people expect:

    read_id  barcode  status  taxid  species  confidence  candidates  lineage

`confidence` is the probability mass on the reported taxon, so it says how much
of the call rests on that one reference: 1.0000 is unambiguous, 0.4 means the
read fits several references nearly equally and the species name is a coin
toss. `status` is C at 0.9 or above, A below it, U for a read that matched
nothing. `candidates` is how many taxa the read matched at all. Reads that reached
the classifier and matched nothing are written as status U, so the file has a
line for every read the classifier saw and the count reconciles with
read_accounting.tsv.

The lineage is pipe-joined, domain first, the way Kraken2 and wf-16s write it.
"""
import csv
import gzip
import os
import sys

# Emu's EM output is a full distribution; a read fitting several references
# equally is normal, not an error. A bare majority is still a coin toss, so
# "C" is reserved for reads whose call really does rest on one reference.
# Below it the read is reported as A: the taxon is still named, with the
# probability beside it, so the reader can decide.
CONFIDENT = 0.9


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def lineage_by_taxid(db_dir):
    """taxid -> (species, genus, pipe-joined lineage) from the database."""
    ranks = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    out = {}
    path = os.path.join(db_dir, "taxonomy.tsv")
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            taxid = (row.get("tax_id") or "").strip()
            if not taxid:
                continue
            names = [(row.get(r) or "").strip() for r in ranks]
            out[taxid] = ((row.get("species") or "").strip(),
                          (row.get("genus") or "").strip(),
                          "|".join(n for n in names if n))
    return out


def read_ids(fastq):
    """Every read id the classifier was given, in order."""
    ids = []
    with open_maybe_gzip(fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                ids.append(line[1:].split()[0].strip())
    return ids


def main():
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    dist_tsv, db_dir, fastq, barcode, out_tsv = sys.argv[1:6]

    if not os.path.exists(dist_tsv):
        sys.exit(f"no read assignments at {dist_tsv}\n"
                 "  Emu writes them only with --keep-read-assignments, which "
                 "nano16s passes when --per-read is set. Re-run with --per-read.")

    tax = lineage_by_taxid(db_dir)
    os.makedirs(os.path.dirname(out_tsv) or ".", exist_ok=True)

    called = set()
    n_conf = n_amb = 0
    with open(dist_tsv) as fh, gzip.open(out_tsv, "wt", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["read_id", "barcode", "status", "taxid", "species",
                    "confidence", "candidates", "lineage"])
        rd = csv.reader(fh, delimiter="\t")
        taxids = next(rd)[1:]
        for row in rd:
            if not row:
                continue
            read_id = row[0]
            best_i = best_p = -1.0
            candidates = 0
            for i, v in enumerate(row[1:]):
                if not v:
                    continue
                try:
                    p = float(v)
                except ValueError:
                    continue
                # Emu leaves vanishing probabilities in the row; they are not
                # matches, they are the tail of the estimate.
                if p < 1e-6:
                    continue
                candidates += 1
                if p > best_p:
                    best_p, best_i = p, i
            called.add(read_id)
            if best_i < 0:
                w.writerow([read_id, barcode, "U", "", "", "", 0, ""])
                continue
            taxid = taxids[best_i]
            species, _genus, lineage = tax.get(taxid, ("", "", ""))
            status = "C" if best_p >= CONFIDENT else "A"
            n_conf += status == "C"
            n_amb += status == "A"
            w.writerow([read_id, barcode, status, taxid, species,
                        f"{best_p:.4f}", candidates, lineage])

        # Reads the classifier was given that matched nothing: Emu leaves them
        # out of the distribution entirely, so a file built from it alone would
        # quietly cover only part of the barcode.
        n_unmapped = 0
        for read_id in read_ids(fastq):
            if read_id not in called:
                w.writerow([read_id, barcode, "U", "", "", "", 0, ""])
                n_unmapped += 1

    total = len(called) + n_unmapped
    print(f"{barcode}: {total:,} reads -> {n_conf:,} confident, {n_amb:,} ambiguous, "
          f"{n_unmapped:,} unclassified")


if __name__ == "__main__":
    main()
