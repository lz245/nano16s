# nano16s user guide

A start-to-finish walkthrough: installing the software, preparing your own
sequencing run, choosing settings that match your amplicon, running the
pipeline, and reading what comes out.

The [README](../README.md) is the quick reference. This guide is the longer
version, written for someone doing this for the first time on their own data.

**Nothing here is specific to a particular dataset.** Every path is written as
a placeholder — substitute your own run directory wherever you see
`/path/to/your_run`.

---

## Contents

1. [What nano16s does](#1-what-nano16s-does)
2. [What you need](#2-what-you-need)
3. [Set up your computer](#3-set-up-your-computer)
    - [A word on the terminal](#a-word-on-the-terminal)
    - [Windows](#windows)
    - [macOS](#macos)
    - [Linux](#linux)
    - [Install Miniforge (all platforms)](#install-miniforge-all-platforms)
    - [Check you are ready](#check-you-are-ready)
4. [Install nano16s](#4-install-nano16s)
5. [Build the reference database](#5-build-the-reference-database)
6. [Verify the install](#6-verify-the-install)
7. [Prepare your data](#7-prepare-your-data)
8. [Choose settings for your amplicon](#8-choose-settings-for-your-amplicon)
9. [Plan the run](#9-plan-the-run-time-disk-memory)
10. [Run it](#10-run-it)
11. [What you get](#11-what-you-get)
12. [Read the results report](#12-read-the-results-report)
13. [Read the performance report](#13-read-the-performance-report)
14. [Quality control: what gets flagged](#14-quality-control-what-gets-flagged)
15. [Use the tables downstream](#15-use-the-tables-downstream)
    - [R](#r)
    - [Python](#python)
    - [phyloseq](#phyloseq)
    - [Excel](#excel)
16. [Interpreting 16S results](#16-interpreting-16s-results)
17. [Re-running and changing settings](#17-re-running-and-changing-settings)
18. [Processing several runs](#18-processing-several-runs)
19. [Troubleshooting](#19-troubleshooting)
20. [Reference](#20-reference)
    - [Command line](#command-line)
    - [Paths](#paths)
    - [Tuning per-rule resources](#tuning-per-rule-resources)
    - [Citing](#citing)

**Never used a terminal before?** Start at section 3 and follow it in order.
It assumes nothing is installed.

**Something went wrong?** Section 19 lists every error message this pipeline
produces, with its cause and its fix.

---

## 1. What nano16s does

You give it the `fastq_pass` folder from an Oxford Nanopore run of full-length
16S rRNA amplicons. It gives you back tables of which organisms are present in
each sample and in what proportion, plus two HTML reports.

```
your fastq_pass/
  barcode01/  ─┐
  barcode02/   │  one directory per sample
  ...          │
               ▼
        merge the FASTQ files in each barcode directory
               ▼
        measure read quality                    (NanoStat)
               ▼
        trim sequencing adapters                (Porechop_ABI)
               ▼
        filter by length and quality            (Chopper)
               ▼
        measure quality again                   (NanoStat)
               ▼
        identify organisms                      (Emu, vs NCBI 16S RefSeq)
               ▼
        abundance tables + reports
```

Each stage runs on every barcode independently, so the work parallelises
across your samples.

**Why Emu.** It estimates abundances with an expectation–maximisation
algorithm rather than assigning each read to its single best database hit. On
error-prone long reads, best-hit assignment loses low-abundance organisms;
Emu recovers them.

### The one thing to understand up front

**One `barcode*` directory is one physical sample.** The pipeline treats each
as a separate sample from beginning to end, and every output column is named
after the barcode directory it came from — `barcode01`, `barcode02`, and so on.

Barcode numbers are *not* sample numbers. If you loaded samples on barcodes 5,
6, 9 and 20, your results have columns `barcode05`, `barcode06`, `barcode09`
and `barcode20` — not 1 through 4. Section 7 covers keeping track of which is
which.

---

## 2. What you need

Section 3 installs the software. This section is what you need to have or
know before that is worth doing.

### A computer

| | |
|---|---|
| **Operating system** | Windows 10/11, macOS (Intel or Apple Silicon), or Linux |
| **RAM** | 8 GB minimum, 16 GB or more comfortable |
| **CPU** | any; more cores means proportionally faster runs |
| **Disk** | about 2–3× your data, plus ~2 GB for the software and ~150 MB for the database and its download cache |
| **Internet** | needed for setup and once to build the database; runs are offline |

On Windows everything runs inside WSL2, which is a real Linux environment
provided by Windows. Section 3 sets it up; you do not need to install Linux
separately or dual-boot.

### Your sequencing data

An Oxford Nanopore run of full-length 16S amplicons, **basecalled and
demultiplexed**, with one directory per barcode. This is what MinKNOW and
Dorado produce by default — the `fastq_pass` folder.

If your reads are in one undivided folder, the run was not demultiplexed, and
that has to happen first. nano16s does not do it.

You do not need your data to start: sections 3 to 6 install and verify
everything against bundled demo data.

### Time

| | |
|---|---|
| Setting up the computer (section 3) | 15–40 minutes, once |
| Installing nano16s (section 4) | 5–15 minutes, once |
| Building the database (section 5) | ~10 minutes, once |
| Verifying it works (section 6) | ~5 minutes |
| A real run | minutes to hours, depending on data size (section 9) |

### What you do not need

No prior bioinformatics experience, no programming, no cluster account, no
Docker, and no administrator rights beyond the initial WSL2 or Miniforge
install. Every tool the pipeline uses is installed for you in section 4.

---

## 3. Set up your computer

Skip to section 4 if you already have conda working in a terminal.

Everything here happens **once per machine**. Follow the part for your
operating system, then the Miniforge step, which is the same for everyone.

### A word on the terminal

The rest of this guide is typed commands. A terminal is a window where you
type a line and press Enter, and the computer replies with text. Commands are
shown in boxes like this:

```bash
echo hello
```

Type or paste the contents, press Enter, and read what comes back. Nothing
here will damage your machine.

---

### Windows

Windows runs nano16s through **WSL2** — Windows Subsystem for Linux — which
gives you a genuine Ubuntu system inside Windows. You do not lose Windows, and
your files stay accessible from both sides.

**1. Open PowerShell as Administrator.** Press Start, type `PowerShell`,
right-click *Windows PowerShell*, choose *Run as administrator*.

**2. Install WSL2 with Ubuntu:**

```powershell
wsl --install
```

> This single command needs Windows 11, or Windows 10 version 2004 or newer.
> On an older Windows 10, `wsl --install` will not be recognised; update
> Windows first, or follow Microsoft's manual WSL2 install steps.

**3. Restart your computer** when it asks.

**4. Finish setting up Ubuntu.** After restarting, an Ubuntu window opens and
asks for a username and password. These are for Linux and are separate from
your Windows login. The password will not appear as you type it — that is
normal.

If no Ubuntu window opens, press Start and run *Ubuntu*.

**5. From here on, use the Ubuntu terminal, not PowerShell.** Every command in
the rest of this guide is typed there. Open it any time from the Start menu.

**6. Install the two tools Ubuntu does not always ship with.** In the Ubuntu
terminal:

```bash
sudo apt update && sudo apt install -y curl git
```

It will ask for the Linux password you chose in step 4. This is harmless if
they are already installed.

> **Where to keep your data.** Work inside the Linux home directory — the
> place the Ubuntu terminal starts in. You can reach Windows drives under
> `/mnt/c/`, but reading data across that boundary is several times slower, and
> it is a common reason for a run taking far longer than it should. Copy your
> sequencing data into the Linux side first.

WSL has two failure modes worth knowing about before a long run — a clock that
drifts from Windows and stops a run near the end, and Windows line endings
breaking scripts. Both are in section 19 with their fixes.

Continue at *Install Miniforge* below.

---

### macOS

**1. Open the Terminal.** Press Cmd-Space, type `Terminal`, press Enter.

**2. Install Apple's command line tools:**

```bash
xcode-select --install
```

A dialog appears; accept it. If it says the tools are already installed, that
is fine — carry on.

Continue at *Install Miniforge* below.

---

### Linux

Open a terminal and make sure `curl` and `git` are present:

```bash
# Debian / Ubuntu
sudo apt update && sudo apt install -y curl git

# Fedora / RHEL
sudo dnf install -y curl git
```

Continue below.

---

### Install Miniforge (all platforms)

Miniforge provides `conda`, which installs and isolates the bioinformatics
tools the pipeline needs. Without it nothing else in this guide will work.

**1. Download and run the installer:**

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

Press Enter to page through the licence, type `yes` to accept, press Enter to
accept the default location, and answer `yes` when it offers to initialise
conda in your shell.

**2. Close the terminal and open a new one.** The change only applies to
terminals started afterwards.

**3. Check it worked:**

```bash
conda --version
```

You should see something like `conda 24.x.x`.

> **`conda: command not found`**
> The shell has not picked up the install. Close the terminal and open a new
> one. If it still fails, run `source ~/.bashrc` (or `source ~/.zshrc` on
> macOS) and try again.

### Check you are ready

```bash
conda --version     # a version number
git --version       # a version number
echo $HOME          # your home directory
```

Two version numbers and a path means the machine is ready — continue to
section 4.

If one of them is missing:

| Missing | Fix |
|---|---|
| `conda` | open a new terminal; if it persists, `source ~/.bashrc` |
| `git` | `sudo apt install -y git` (Linux/WSL), or `xcode-select --install` (macOS) |

---

## 4. Install nano16s

```bash
git clone https://github.com/lz245/nano16s.git
cd nano16s
bash install.sh
conda activate nano16s
```

Check it worked:

```bash
nano16s --version
```

> **`nano16s: command not found`**
> You need `conda activate nano16s` in **every new terminal**. This is the
> single most common problem people hit. If a command that worked yesterday
> fails today, this is almost always why.

Consider adding it to your shell profile so it happens automatically:

```bash
echo 'conda activate nano16s' >> ~/.bashrc      # or ~/.zshrc
```

---

## 5. Build the reference database

Organisms can only be identified against a reference. Build one now — it takes
about ten minutes and downloads roughly 100 MB from NCBI. **You do this once**,
not per run.

```bash
nano16s db build
```

This builds from the current NCBI 16S RefSeq Targeted Loci collection — about
28,000 curated sequences across 21,000 taxa.

**Why not use Emu's own database?** Emu ships one built from NCBI in September
2020 and never refreshed. Species described since then are simply absent from
it, and a read from an absent species gets assigned to whatever relative *is*
present. Building your own avoids that.

### Managing databases

```bash
nano16s db list                     # what you have
nano16s db build --version 2026.07  # name a build explicitly
```

Databases install to `~/.nano16s/db/<version>/` and never overwrite each other.
That matters for reproducibility: a result from six months ago can still be
re-run against the database that produced it. The database version is recorded
in every report.

Runs use the newest installed database unless you say otherwise:

```bash
nano16s -d /path/to/your_run/fastq_pass --db ~/.nano16s/db/2026.07/ncbi_16s
```

Set `NANO16S_DB` to move the whole database root somewhere else, such as a
shared drive.

---

## 6. Verify the install

Before touching your own data, run the bundled six-barcode demo. It takes about
five minutes.

```bash
nano16s test
```

If it prints **Install verified**, everything works: the environment, every
tool, the database, and the reports.

Do this after installing, and again after any change to your conda environment.
It is much easier to debug a broken install on demo data than three hours into a
real run.

### Trying it on a full sequencing run

The bundled demo is deliberately tiny — six barcodes, enough to prove the
install works. If you want to see what nano16s does with a real run before you
commit your own data, the five runs it was developed against are published:

**<https://doi.org/10.5281/zenodo.21998286>**

| Dataset | Barcodes | Download | Unpacked |
|---|---:|---:|---:|
| `Flongle_Demo01` | 24 | 416 MB | 401 MB |
| `Flongle_Demo02` | 16 | 1.0 GB | 973 MB |
| `Flongle_Demo03` | 24 | 708 MB | 682 MB |
| `MinION_Demo01` | 24 | 4.6 GB | 4.3 GB |
| `PromethION_Demo01` | 24 | 5.2 GB | 4.9 GB |

Each archive unpacks to a single directory holding a `fastq_pass/` folder with
one `barcode*` subdirectory per sample, exactly as MinKNOW writes it, plus that
run's original MinKNOW report — so it doubles as a worked example of the layout
section 7 describes. `fastq_fail` and `unclassified` reads are not included.

Start with `Flongle_Demo01`: it is the smallest, and takes about 45 minutes on
a 20-core machine — see the table in section 9 for the others.

```bash
mkdir -p ~/data && cd ~/data
curl -L -O https://zenodo.org/records/21998287/files/Flongle_Demo01.zip
unzip Flongle_Demo01.zip

nano16s -d ~/data/Flongle_Demo01/fastq_pass -o ~/nano16s_out/Flongle_Demo01
```

Leave room for the download, the unpacked reads, and roughly 3× the unpacked
size for the run itself — so about 2 GB all told for `Flongle_Demo01`, and
around 25 GB for `PromethION_Demo01`. The whole set is about 12 GB downloaded.

`CHECKSUMS.txt` in the record lets you confirm a download arrived intact:

```bash
curl -L -O https://zenodo.org/records/21998287/files/CHECKSUMS.txt
sha256sum -c CHECKSUMS.txt --ignore-missing
```

The DOI always resolves to the current version. The data is licensed CC BY 4.0
— cite it if you use it in your own work (see section 20).

---

## 7. Prepare your data

This is where most first runs go wrong, and it is worth five minutes of
checking.

### The layout nano16s expects

```
/path/to/your_run/
└── fastq_pass/              ← point -d at THIS directory
    ├── barcode01/
    │   ├── something.fastq.gz
    │   └── another.fastq.gz     (several files per barcode is normal)
    ├── barcode02/
    │   └── ...
    └── barcode20/
```

Rules:

- `-d` points at the directory **containing** the `barcode*` directories, not
  at a barcode directory and not at the run directory above it.
- Barcode directories must be named `barcode` followed by digits. This is what
  MinKNOW and Dorado produce.
- Barcode numbers need not be contiguous. `barcode05`, `barcode06`, `barcode20`
  is perfectly fine.
- Each barcode directory needs at least one `.fastq.gz` file. Multiple files
  are merged automatically.

### Check before you run

```bash
# 1. How many barcodes will be processed?
ls -d /path/to/your_run/fastq_pass/barcode* | wc -l

# 2. Which ones? (confirms numbering matches what you loaded)
ls -d /path/to/your_run/fastq_pass/barcode* | xargs -n1 basename

# 3. Does every barcode actually contain data?
for d in /path/to/your_run/fastq_pass/barcode*/; do
    n=$(ls "$d"*.fastq.gz 2>/dev/null | wc -l)
    echo "$(basename "$d"): $n files"
done
```

A barcode showing `0 files` will stop the run. Either the copy was incomplete,
or that barcode genuinely produced nothing — in which case remove the empty
directory and continue.

### Situations you may need to handle

**Your files are `.fastq`, not `.fastq.gz`.** Compress them:

```bash
gzip /path/to/your_run/fastq_pass/barcode*/*.fastq
```

**Everything is in one folder with no barcode directories.** The run was not
demultiplexed. Demultiplex it first (with Dorado or the MinKNOW re-basecalling
options); nano16s does not do this step.

**`fastq_pass` is nested deeper**, for instance inside a date-and-flowcell
folder. That is fine — give the full path. Find it with:

```bash
find /path/to/your_run -type d -name fastq_pass
```

**You also have a `fastq_fail` folder.** Ignore it. Those reads failed the
basecaller's quality filter.

**Your data is on a network drive or an external disk.** Copy it to a local
disk first. The pipeline reads every file several times and network latency
dominates the runtime. On WSL2 specifically, keep data under the Linux home
directory and *not* under `/mnt/c/` — crossing the Windows filesystem boundary
is several times slower.

### Record which barcode is which sample

The pipeline cannot know your sample names, so every output column is a barcode
identifier. Before you forget, write the mapping down next to your results:

```bash
cat > /path/to/your_run/barcode_map.csv <<'CSV'
barcode,sample
barcode01,Field_plot_A_rep1
barcode02,Field_plot_A_rep2
barcode05,Control_soil
CSV
```

Section 15 shows how to apply it when loading the tables. Doing this at the
start rather than at analysis time saves real confusion — barcode numbering and
sample numbering rarely line up.

---

## 8. Choose settings for your amplicon

### Length window — the setting that matters most

Reads shorter than `--min-length` or longer than `--max-length` are discarded.
The defaults suit the full-length 16S gene:

| | Default |
|---|---|
| `--min-length` | 1000 |
| `--max-length` | 2000 |

Full-length 16S is about 1,500 bp, so this keeps near-full-length reads and
drops fragments and concatemers.

**If you amplified something else, change this.** A different region left at
the default window silently discards most of your data, and the run will look
like it worked.

| What you amplified | Approximate product | Suggested window |
|---|---|---|
| Full-length 16S (27F–1492R) | ~1,500 bp | 1000–2000 (default) |
| 16S + 23S rRNA operon | ~4,500 bp | 3500–5500 |
| V3–V4 | ~460 bp | 300–700 |
| V1–V9 with long primers | ~1,600 bp | 1200–2000 |

Not sure what you have? Run the pipeline on a couple of barcodes with a wide
window, then look at the median read length in the report and narrow it:

```bash
nano16s -d /path/to/your_run/fastq_pass -o /tmp/length_check \
        --min-length 200 --max-length 10000
```

### Quality

`--min-quality` (default `10`) drops reads whose mean Phred quality is below
the threshold. Q10 means roughly 90% base accuracy — a reasonable floor for
modern nanopore chemistry.

Raise it to 12 or 15 if you have reads to spare and want cleaner
classification. Watch the retention figures in the report: if you are throwing
away more than about 20% of reads, you are being too strict for your data.

### Cores

`-c` defaults to every core but one. Lower it if you need the machine for
something else:

```bash
nano16s -d /path/to/your_run/fastq_pass -c 4
```

---

## 9. Plan the run (time, disk, memory)

### How long

Runtime scales with the number of reads, not the number of barcodes. As
measured on a 20-core workstation:

| Reads in the run | Barcodes | Elapsed | Dataset (section 6) |
|---|---|---|---|
| ~200,000 | 24 | ~45 min | `Flongle_Demo01` |
| ~340,000 | 24 | ~75 min | `Flongle_Demo03` |
| ~510,000 | 16 | ~85 min | `Flongle_Demo02` |
| ~2,300,000 | 24 | ~4h 20m | `MinION_Demo01` |
| ~3,200,000 | 24 | ~8h 25m | `PromethION_Demo01` |

Treat these as a rough guide — a machine with a quarter of the cores takes
substantially longer. Two stages dominate: **Porechop**, which infers adapter
sequences from your data rather than assuming them, and **Emu**, which does the
classification.

For a long run, start it in a way that survives losing your terminal:

```bash
nohup nano16s -d /path/to/your_run/fastq_pass -o my_results -y \
    > my_run.log 2>&1 &

tail -f my_run.log        # watch progress; Ctrl-C stops watching, not the run
```

### How much disk

Between 2× and 3× your input — 2.0× to 2.3× across the five demo runs. nano16s
budgets the higher figure, checks free space before starting, and warns you if
it looks tight.

Afterwards you can reclaim most of it:

```bash
rm -rf my_results/01_merged my_results/03_trimmed
```

Keep `04_filtered/` if you might re-run the classification step against a newer
database; it is the input to Emu.

The database build also leaves its NCBI downloads in `~/.nano16s/cache`, about
100 MB. It is only needed if you rebuild, so it is safe to delete.

### How much memory

Peak memory is dominated by Emu and scales with database size, not with your
number of reads — roughly 500 MB to 2.5 GB per concurrent job. The performance
report records the actual peak for every run.

---

## 10. Run it

Preview first. This lists the steps without executing them, and catches a bad
path or a missing database in seconds rather than minutes:

```bash
nano16s -d /path/to/your_run/fastq_pass -o my_results -n
```

Then run it:

```bash
nano16s -d /path/to/your_run/fastq_pass -o my_results
```

Before starting, nano16s prints a summary — input, barcode count, database,
filter settings, cores — and warns if free disk looks insufficient. Add `-y` to
skip the confirmation when running unattended.

A full example with non-default settings:

```bash
nano16s \
    -d /path/to/your_run/fastq_pass \
    -o /path/to/results/my_experiment \
    --min-length 1300 \
    --max-length 1800 \
    --min-quality 12 \
    -c 8
```

**If it stops partway** — a crash, a power cut, a closed laptop — run exactly
the same command again. Completed work is detected and skipped, and the run
picks up where it stopped.

---

## 11. What you get

```
my_results/
├── nano16s_report.html          ← open this first
├── performance_report.html      the run itself: timings, machine, QC flags
├── performance_summary.csv      per-barcode numbers behind that report
├── performance.json             machine-readable, for comparing runs
├── preprocessing_summary.csv    reads surviving each stage, per barcode
├── benchmarks/                  per-job wall time, CPU time, peak memory
├── 01_merged/                   intermediates — safe to delete when done
├── 02_nanostat_raw/
├── 03_trimmed/
├── 04_filtered/
├── 05_nanostat_filtered/
├── 06_emu_output/               per-barcode classification
└── 07_emu_combined/             ← the results you will analyse
    ├── emu-combined-species.tsv         relative abundance
    ├── emu-combined-species-counts.tsv  estimated read counts
    ├── emu-combined-genus.tsv
    ├── emu-combined-genus-counts.tsv
    ├── emu-combined-phylum.tsv
    └── emu-combined-phylum-counts.tsv
```

**Relative abundance vs counts.** Abundance tables give each taxon's proportion
of the sample, summing to 1 per column. Counts tables give Emu's estimated
number of reads. Use abundances to compare composition between samples; use
counts for methods that expect count data, such as differential-abundance
testing.

Both reports are self-contained single files — no internet needed to view them,
safe to email or attach to a manuscript.

---

## 12. Read the results report

### Opening it

The report is one self-contained HTML file — no internet connection, no
external files, safe to email. Open it the way you would open any web page:

| Where you are | Command |
|---|---|
| Linux desktop | `xdg-open nano16s_report.html` |
| macOS | `open nano16s_report.html` |
| Windows, via WSL | `explorer.exe nano16s_report.html` |
| No desktop (SSH, server) | copy it to your own machine first — see below |

On Ubuntu, `open` is an alias for `xdg-open`, so both work there.

You can also double-click the file in your file manager. From WSL, the reports
appear under `\\wsl.localhost\Ubuntu\home\<you>\...` in Windows Explorer.

Over SSH there is no desktop for the server to open a window on, so bring the
file to you:

```bash
scp you@server:/path/to/my_results/nano16s_report.html .
```

**If the report is under `/tmp`, copy it to your home directory first.**

```bash
cp /tmp/nano16s_test_*/nano16s_report.html ~/
xdg-open ~/nano16s_report.html
```

Two reasons. `/tmp` is cleared on reboot. More immediately, on Ubuntu the
default Firefox is a **snap**, and a snap gets its own private `/tmp` — so the
browser genuinely cannot see a file that `ls` shows you in the terminal, and
reports *File not found* for a path that exists. Anything under your home
directory is readable.

Older versions of `nano16s test` wrote there; real runs go wherever you point
`-o`, so they were never affected.

### What is in it

**Read funnel.** Reads per barcode before and after filtering. What you want is
consistency: barcodes retaining broadly similar proportions. One barcode
retaining far less than the rest points at a problem with that sample or that
library.

**Composition charts.** Relative abundance per barcode at species and genus
level. The top taxa get consistent colours across all barcodes, so a colour
means the same organism in every row.

**Per-barcode table.** Read counts, median length, and median quality before and
after filtering.

**Methods paragraph.** Tool versions, parameters, and database version in prose
you can paste into a manuscript. Everything needed to make the run reproducible
is recorded here.

---

## 13. Read the performance report

Open `performance_report.html`. This one is about the run rather than the
biology — reach for it when something took longer than expected, or a barcode
looks wrong.

**Headline figures.** Elapsed time, total CPU time, average number of jobs
running at once, percentage of your cores used, and peak memory.

**System.** The machine the run was measured on: CPU model, cores, memory,
operating system, kernel, architecture. Timings only compare meaningfully
between runs on comparable hardware, so this travels with the report.

**Machine use.** Whether the run kept your machine busy, and if not, what to do
about it. The usual cause of low utilisation is that each stage reserves more
threads than it can use: a job cannot start until its full thread count is
free, so large per-rule thread counts run fewer barcodes at once. Lowering
`resources.*.cpus` in `config/config.yaml` trades per-job speed for more
barcodes in flight, which is normally the better trade when barcodes outnumber
cores.

**Where the time went.** Wall time and CPU time per stage. CPU time larger than
wall time means the stage used several threads; the ratio between them tells
you how well it used them.

**Timeline.** Every job placed by when it actually ran, one row per barcode.
Gaps are idle capacity.

**Worth checking.** QC flags — see the next section.

**Per barcode.** Timings first, then read counts and quality. Click any column
heading to sort.

---

## 14. Quality control: what gets flagged

The performance report checks every barcode against fixed thresholds and
explains anything that fails. These are absolute, not relative to the rest of
your run — a check that only compares barcodes to each other cannot fire when
every barcode is equally bad, which is the case most worth catching.

| Flag | Threshold | What it usually means |
|---|---|---|
| **Low retention** | under 80% of reads survive filtering | Your length window does not match the amplicon. Check the median read length and adjust `--min-length` / `--max-length`. |
| **Low depth** | fewer than 1,000 raw reads | Too few reads for reliable proportions. Treat that barcode's abundances as indicative only. |
| **Below run median** | under a quarter of the median depth | Barcoding imbalance — that library was under-represented in the pool. Usually a loading issue, not a data problem. |
| **Low quality** | filtered median below Q12 | Unusual after filtering. Points at a basecalling or chemistry problem. |
| **Length outside window** | filtered median outside your configured range | Your amplicon is not the length you configured for. |

**No flags is the expected outcome for a healthy run.** If everything is quiet,
the reads look the way full-length 16S data should.

---

## 15. Use the tables downstream

The combined tables are plain tab-separated text. Rows are taxa, columns are
barcodes.

### R

```r
ab <- read.delim("my_results/07_emu_combined/emu-combined-genus.tsv",
                 check.names = FALSE)

# apply your barcode -> sample mapping
map <- read.csv("/path/to/your_run/barcode_map.csv")
idx <- match(colnames(ab), map$barcode)
colnames(ab)[!is.na(idx)] <- map$sample[idx[!is.na(idx)]]
```

### Python

```python
import pandas as pd

ab = pd.read_csv("my_results/07_emu_combined/emu-combined-genus.tsv", sep="\t")

mapping = pd.read_csv("/path/to/your_run/barcode_map.csv")
ab = ab.rename(columns=dict(zip(mapping["barcode"], mapping["sample"])))
```

### phyloseq

Use the counts table as the OTU table and the lineage columns as taxonomy.
Emu writes the full lineage — species through superkingdom — before the sample
columns, so split the frame at the first barcode column.

### Excel

Open the `.tsv` directly, or import as tab-delimited. Watch out for Excel
converting taxon names that look like dates.

---

## 16. Interpreting 16S results

**Genus is the defensible resolution.** Species-level assignment from 16S is
genuinely hard — many genera contain species whose 16S genes are near-identical,
and nanopore error rates make it harder. Species tables are provided because
they are sometimes informative, but conclusions are safer at genus level.

**Relative abundance is compositional.** Proportions sum to 1, so one taxon
rising means others fall by arithmetic, not biology. Use methods designed for
compositional data when testing for differences.

**16S copy number varies between organisms**, from one to over fifteen. A
species with many copies is over-represented relative to its true cell
abundance. nano16s does not correct for this; no tool does it reliably.

**Absence of evidence.** A taxon missing from your results may be absent from
the sample, or absent from the reference database, or below detection at your
sequencing depth.

**Your database version is part of your result.** It is recorded in every
report. Cite it alongside the tool versions.

---

## 17. Re-running and changing settings

nano16s tracks what has been done. Re-running the same command finishes in
seconds; changing a setting re-runs only what that setting affects.

```bash
# resume an interrupted run — same command, nothing lost
nano16s -d /path/to/your_run/fastq_pass -o my_results

# different filtering: re-runs from the filter step onward,
# reusing the merged and trimmed reads
nano16s -d /path/to/your_run/fastq_pass -o my_results --min-length 1200

# newer database: re-runs classification only
nano16s -d /path/to/your_run/fastq_pass -o my_results \
        --db ~/.nano16s/db/2026.09/ncbi_16s
```

**Compare settings side by side** by writing to separate output directories:

```bash
nano16s -d /path/to/your_run/fastq_pass -o results_q10 --min-quality 10
nano16s -d /path/to/your_run/fastq_pass -o results_q15 --min-quality 15
```

Trimming is repeated for each, so this costs a full run rather than a partial
one.

**To start completely fresh**, delete the output directory. Re-running does not
discard existing work by design, which is what makes resuming possible.

---

## 18. Processing several runs

Give each run its own output directory and loop:

```bash
for run in /path/to/runs/*/; do
    name=$(basename "$run")
    [ -d "$run/fastq_pass" ] || continue
    echo "=== $name ==="
    nano16s -d "$run/fastq_pass" -o "/path/to/results/$name" -y \
        || echo "$name FAILED — continuing"
done
```

Points worth noting:

- **Run them one at a time.** Concurrent runs compete for the same cores and
  finish no sooner, and a failure is harder to attribute.
- `|| echo ... ` keeps one bad run from ending the loop.
- Because completed work is skipped, re-running the loop after a failure costs
  only the runs that did not finish.
- **Keep settings identical across runs you intend to compare.** A different
  length window or database makes the results incomparable.

---

## 19. Troubleshooting

**`nano16s: command not found`**
Run `conda activate nano16s`. Needed in every new terminal.

**`no barcode* directories inside ...`**
`-d` is one level too high or too low. It wants the directory that directly
contains `barcode01/`. Find it with
`find /path/to/your_run -type d -name fastq_pass`.

**`no Emu database found`**
Run `nano16s db build` (about ten minutes, once).

**`no .fastq.gz files in ...`**
That barcode directory is empty, or holds uncompressed `.fastq`. Compress them,
or remove the directory if the barcode genuinely produced nothing.

**A barcode kept almost no reads**
Its reads fall outside your length window. Check the median read length in the
report and widen `--min-length` / `--max-length` if that length is expected for
your amplicon.

**The run stopped partway**
Run the same command again. Completed steps are skipped.

**Porechop is slow**
Expected — it is the slowest stage by a wide margin because it infers adapter
sequences from your data instead of assuming them. Budget a few minutes per
barcode.

**The run is using less of my CPU than expected**
See the Machine use section of the performance report, and section 13.

**The browser says *File not found* for a report that `ls` shows is there**

Almost always a report under `/tmp` opened with Ubuntu's snap Firefox. A snap
runs with its own private `/tmp`, so the file you can see in the terminal is
genuinely not in the `/tmp` the browser sees. Confirm with:

```bash
snap list firefox                      # a version listed means it is a snap
readlink -f "$(command -v firefox)"    # /snap/bin/firefox means the same
```

Copy the report somewhere under your home directory and open it there:

```bash
cp /tmp/nano16s_test_*/nano16s_report.html ~/
xdg-open ~/nano16s_report.html
```

The same applies to `performance_report.html`, and to Chromium installed as a
snap.

**`Gtk-Message: Not loading module "atk-bridge"` when opening a report**

Not an error, and nothing to install. `Gtk-Message:` is informational, and the
browser prints it *after* it has already started — the report has almost
certainly opened, possibly behind the terminal window or on another workspace.

There is no `atk-bridge` package, so `sudo apt install atk-bridge` cannot work.
(`install` on its own is a file-copying command, which is why
`sudo install atk-bridge` answers `missing destination file operand`.)

If no window appeared at all, the file is fine — the question is which
application your desktop chose for it:

```bash
xdg-mime query default text/html     # what is registered to open it
xdg-open nano16s_report.html         # try again
firefox nano16s_report.html          # or name a browser directly
```

If the first command prints nothing, no application is registered for HTML —
naming a browser directly, as in the third line, is the fix.

Over SSH there is no desktop to open anything on; copy the file to your own
machine instead (section 12).

### On WSL2

**`has older modification time`, or `the counts table … is empty`**
WSL's clock drifts from the Windows host and can resynchronise mid-run, leaving
a file with a timestamp behind its own input. Snakemake reads that as a
corrupted build, deletes the output and stops — often near the end of a long
run.

From PowerShell:

```powershell
wsl --shutdown
```

Then re-run the same command. Completed work is kept, so it finishes quickly.
Check for drift by comparing `date` in Linux with `Get-Date` in PowerShell;
more than a second or two apart is the cause. Keeping the machine awake during
a run avoids it.

**Everything is slow**
Check your data is not under `/mnt/c/`. Crossing the Windows filesystem
boundary is several times slower than the Linux filesystem. Copy the run into
your Linux home directory first.

**`syntax error near unexpected token $'{\r'`**
A file has Windows line endings, usually from editing the repository through
Windows. Re-clone, or run `dos2unix` on the affected file.

### Still stuck

Open an issue at
<https://github.com/lz245/nano16s/issues>. Include the output of
`nano16s --version`, your operating system, the exact command, and the error.
The `performance.json` from a failed run is small and records the machine and
settings, which usually answers the first three questions at once.

---

## 20. Reference

### Command line

| Option | Default | What it does |
|---|---|---|
| `-d, --input-dir` | *required* | directory containing `barcode*` directories |
| `-o, --output-dir` | `results` | where output goes |
| `--db` | newest installed | Emu database directory |
| `--min-length` | `1000` | shortest read to keep, bp |
| `--max-length` | `2000` | longest read to keep, bp |
| `--min-quality` | `10` | minimum mean Phred quality |
| `-c, --cores` | all but one | CPU cores to use |
| `-n, --dry-run` | | list the steps and stop |
| `-y, --yes` | | skip confirmation prompts |
| `-h, --help` | | full help |

Subcommands: `nano16s db build`, `nano16s db list`, `nano16s test`,
`nano16s --version`.

### Paths

| Path | What |
|---|---|
| `~/.nano16s/db/<version>/` | reference databases |
| `~/.nano16s/cache/` | NCBI downloads kept for rebuilding, ~100 MB, safe to delete |
| `config/config.yaml` | defaults, including per-rule CPU and memory |
| `<output>/benchmarks/` | per-job timing and memory records |

`NANO16S_DB` overrides the database root.

### Tuning per-rule resources

`config/config.yaml` sets threads and memory per stage. The thread counts
determine how many barcodes run concurrently: a job cannot start until its full
thread count is free, so a large value means fewer barcodes in flight. If the
performance report shows low core utilisation, lowering these is the lever.

```yaml
resources:
  porechop:
    cpus: 8            # the slowest stage
    mem_mb: 8000
  emu:
    cpus: 8
    mem_mb: 8000
```

Change settings there rather than on the command line when you want them to
apply to every run.

### Citing

Cite the underlying tools, which do the actual work:

- **Emu** — Curry et al. (2022) *Nature Methods* 19:845–853
- **Porechop_ABI** — Bonenfant et al. (2023) *Bioinformatics Advances* 3:vbac085
- **Chopper / NanoStat** — De Coster & Rademakers (2023) *Bioinformatics* 39:btad311
- **minimap2** — Li (2018) *Bioinformatics* 34:3094–3100
- **Snakemake** — Mölder et al. (2021) *F1000Research* 10:33

See [CITATION.cff](../CITATION.cff) for nano16s itself, and record the database
version from your report.

If you use the demo datasets from section 6, cite them as well:

- Zhang, L. & Adapa, P. D. *nano16s demo datasets: Oxford Nanopore full-length
  16S rRNA amplicon sequencing runs (Flongle, MinION, PromethION)*. Zenodo (2026).
  <https://doi.org/10.5281/zenodo.21998286>
