# nano16s

Taxonomic profiling of Oxford Nanopore full-length 16S rRNA amplicon data, on
an ordinary computer.

Point it at the `fastq_pass` folder your sequencer produced and it returns
abundance tables and a summary report. No cluster, no container, no
bioinformatics setup beyond one install command.

```bash
nano16s -d /path/to/fastq_pass -o my_results
```

---

## What it does

```
fastq_pass/
  barcode01/  ─┐
  barcode02/   │
  ...          │
               ▼
        merge per barcode
               ▼
        read QC  (NanoStat)
               ▼
        adapter trimming  (Porechop_ABI)
               ▼
        length + quality filter  (Chopper)
               ▼
        read QC again  (NanoStat)
               ▼
        classification  (Emu, against NCBI 16S RefSeq)
               ▼
        combined abundance tables + HTML report
```

Emu estimates abundance with an expectation-maximization algorithm rather than
assigning each read to a single best hit. On error-prone long reads that
recovers low-abundance organisms a best-hit approach tends to lose.

---

## Requirements

- macOS (Intel or Apple Silicon) or Linux. Windows works through WSL2.
- [Miniforge](https://github.com/conda-forge/miniforge#install) or any conda.
- About 3× your input size in free disk, and 8 GB RAM.
- Internet once, to build the reference database.

Every tool installs natively — no Rosetta, no emulation, no Docker.

CI exercises Linux and Apple Silicon on every change. Intel macOS is
supported and all dependencies publish `osx-64` builds, but it is not
covered by CI — Intel runners are too scarce to be useful. If you are on
an Intel Mac, `nano16s test` after installing is worth the two minutes.

## Install

```bash
git clone https://github.com/lz245/nano16s-local.git
cd nano16s-local
bash install.sh
conda activate nano16s
```

Then build the reference database. This downloads about 100 MB from NCBI and
takes roughly ten minutes. You only do it once.

```bash
nano16s db build
```

Check the install on the bundled demo data — six barcodes, about two minutes:

```bash
nano16s test
```

If that prints **Install verified**, you are ready.

## Use

```bash
nano16s -d /path/to/fastq_pass -o my_results
```

`-d` must point at the folder that *contains* the `barcode01/`, `barcode02/`
… directories, which is usually called `fastq_pass`.

Common adjustments:

```bash
# stricter quality, tighter length window
nano16s -d ./fastq_pass --min-quality 12 --min-length 1300 --max-length 1800

# see what would run, without running it
nano16s -d ./fastq_pass -n

# limit CPU use
nano16s -d ./fastq_pass -c 4
```

| Option | Default | What it does |
|---|---|---|
| `-d, --input-dir` | *required* | folder containing `barcode*` directories |
| `-o, --output-dir` | `results` | where output goes |
| `--db` | newest installed | Emu database directory |
| `--min-length` | `1000` | shortest read to keep, bp |
| `--max-length` | `2000` | longest read to keep, bp |
| `--min-quality` | `10` | minimum mean Phred quality |
| `-c, --cores` | all but one | CPU cores to use |
| `-n, --dry-run` | | list the steps and stop |
| `-y, --yes` | | skip confirmation prompts |

Full list: `nano16s --help`.

### Choosing the length window

Full-length 16S is about 1,500 bp, so the default 1,000–2,000 bp window keeps
near-full-length reads and drops fragments and concatemers. If you amplified a
different region, set the window to match — leaving it at the default will
silently discard most of your data.

## Output

```
my_results/
├── nano16s_report.html          ← open this first
├── preprocessing_summary.csv    reads surviving each stage, per barcode
├── 01_merged/                   intermediates; safe to delete when done
├── 02_nanostat_raw/
├── 03_trimmed/
├── 04_filtered/
├── 05_nanostat_filtered/
├── 06_emu_output/               per-barcode Emu results
└── 07_emu_combined/
    ├── emu-combined-species.tsv         relative abundance
    ├── emu-combined-species-counts.tsv  read counts
    ├── emu-combined-genus.tsv
    ├── emu-combined-genus-counts.tsv
    ├── emu-combined-phylum.tsv
    └── emu-combined-phylum-counts.tsv
```

The combined tables are plain tab-separated files: rows are taxa, columns are
barcodes. They load directly into R (`read.delim`), Python (`pandas.read_csv`),
phyloseq, or Excel.

The HTML report is self-contained — one file, no internet needed, safe to
email. It carries the tool versions, parameters, and database version used, in
a form you can paste into a Methods section.

**A note on interpretation.** Species-level assignment from 16S is genuinely
hard: many genera contain species whose 16S genes are near-identical. Genus
level is the more defensible resolution for most nanopore 16S work. Both are
reported so you can decide.

## The database

`nano16s db build` builds an Emu database from the current NCBI 16S RefSeq
Targeted Loci collection — about 28,000 curated sequences across 21,000 taxa.

Emu's own default database was built from NCBI in September 2020 and has not
been refreshed since, so species described in the last several years are simply
absent from it. Building your own fixes that.

```bash
nano16s db build                    # newest NCBI release, stamped with today's month
nano16s db build --version 2026.07  # name it explicitly
nano16s db list                     # what's installed
```

Databases install to `~/.nano16s/db/<version>/` and never overwrite each other,
so re-running an old analysis against its original database stays possible.
Override the location with `--db` or the `NANO16S_DB` environment variable.

## Troubleshooting

**`nano16s: command not found`** — activate the environment: `conda activate nano16s`.

**`no barcode* directories inside ...`** — `-d` is pointing one level too high or
too low. It wants the folder that directly contains `barcode01/`.

**`no Emu database found`** — run `nano16s db build` first.

**A barcode retained almost no reads** — its reads are probably outside the
length window. Check the report's filtering section, then widen
`--min-length`/`--max-length` if that's expected for your amplicon.

**The run stopped partway** — just run the same command again. Completed steps
are skipped and it picks up where it left off.

**Porechop is slow** — it is the slowest stage by a wide margin, because it
infers adapter sequences from your data rather than assuming them. Budget a few
minutes per barcode.

## Development

Unit tests cover the parsing functions — the ones that read Emu's tables and
NCBI's taxonomy. They need only `pytest`, no bioinformatics tools and no
database, and run in under a second:

```bash
python -m pip install pytest
python -m pytest test/test_parsers.py -v
```

`nano16s test` is the other half: an end-to-end run on the bundled demo that
verifies an actual installation. Between them, the unit tests catch parsing
regressions in seconds and the demo run catches everything else.

CI runs the unit tests on Linux and Apple Silicon on every push, and also
installs from scratch on both to check that the documented install path
still works. A full pipeline run — NCBI database build plus the demo —
runs weekly; its job is to catch a bioconda dependency release breaking
the environment before a user hits it.

## Citing

If you use nano16s, please cite the underlying tools, which do the actual work:

- **Emu** — Curry et al. (2022) *Nature Methods* 19:845–853
- **Porechop_ABI** — Bonenfant et al. (2023) *Bioinformatics Advances* 3:vbac085
- **Chopper / NanoStat** — De Coster & Rademakers (2023) *Bioinformatics* 39:btad311
- **minimap2** — Li (2018) *Bioinformatics* 34:3094–3100
- **Snakemake** — Mölder et al. (2021) *F1000Research* 10:33

See [CITATION.cff](CITATION.cff) for nano16s itself.

## License

MIT. See [LICENSE](LICENSE).
