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
7. [Try it on real data](#7-try-it-on-real-data)
8. [Prepare your data](#8-prepare-your-data)
    - [The layout nano16s expects](#the-layout-nano16s-expects)
    - [Check before you run](#check-before-you-run)
    - [Situations you may need to handle](#situations-you-may-need-to-handle)
    - [Record which barcode is which sample](#record-which-barcode-is-which-sample)
9. [Choose settings for your amplicon](#9-choose-settings-for-your-amplicon)
10. [Plan the run (time, disk, memory)](#10-plan-the-run-time-disk-memory)
11. [Run it](#11-run-it)
12. [What you get](#12-what-you-get)
13. [Read the results report](#13-read-the-results-report)
14. [Read the performance report](#14-read-the-performance-report)
15. [Quality control: what gets flagged](#15-quality-control-what-gets-flagged)
16. [Use the tables downstream](#16-use-the-tables-downstream)
    - [R](#r)
    - [Python](#python)
    - [phyloseq](#phyloseq)
    - [Excel](#excel)
17. [Interpreting 16S results](#17-interpreting-16s-results)
18. [Re-running and changing settings](#18-re-running-and-changing-settings)
19. [Processing several runs](#19-processing-several-runs)
    - [Leave it running overnight](#leave-it-running-overnight)
    - [Comparing runs](#comparing-runs)
20. [Troubleshooting](#20-troubleshooting)
    - [Setup](#setup)
    - [Your data](#your-data)
    - [During the run](#during-the-run)
    - [Opening the reports](#opening-the-reports)
    - [On WSL2](#on-wsl2)
    - [Still stuck](#still-stuck)
21. [Reference](#21-reference)
    - [Command line](#command-line)
    - [Paths](#paths)
    - [Tuning per-rule resources](#tuning-per-rule-resources)
    - [Citing](#citing)

**Reading this for the first time?** Sections 1 to 11 are a walkthrough: they
start with nothing installed and end with a finished run. Follow them in order.

**Already have results?** Sections 12 to 17 explain what the files and reports
contain and how to analyse them.

**Something went wrong?** Section 20 lists the error messages this pipeline
produces, each with its cause and its fix.

**Looking one thing up?** Section 21 is the reference — every command-line
option, every path, and how to change the defaults.

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
and `barcode20` — not 1 through 4. Section 8 covers keeping track of which is
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

> **On a managed work or university machine**, the WSL2 step needs
> Administrator rights, and the Microsoft Store and external DNS are often
> blocked by policy. Section 3 and section 20 cover both.

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
| A real run | minutes to hours, depending on data size (section 10) |

### What you do not need

No prior bioinformatics experience, no programming, no cluster account, and no
Docker. The only step needing Administrator rights is the WSL2 install on
Windows; Miniforge and everything after it install into your own home
directory. Every tool the pipeline uses is installed for you in section 4.

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
right-click *Windows PowerShell*, choose *Run as administrator*. This is the
only step that needs Administrator rights.

**2. Install WSL2 with Ubuntu:**

```powershell
wsl --install
```

> If this sits at **0%** for more than a few minutes, it is blocked from the
> Microsoft Store rather than installing slowly — common on work machines. Close
> it and install in two steps instead, restarting the computer between them:
>
> ```powershell
> wsl --install --no-distribution --web-download
> wsl --install -d Ubuntu --web-download
> ```
>
> On Windows 10 older than version 2004, `wsl --install` is not recognised at
> all; update Windows first.

**3. Restart your computer.** Required even if the install reported success —
until you restart, WSL commands answer
`Wsl/WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED`, which looks like a failure and is
not one.

**4. Set your Linux username and password.** Ubuntu asks the first time it
starts. These are for Linux only, unrelated to your Windows login, and **the
password stays invisible as you type** — not even dots. Depending on how Ubuntu
started, it asks either in a separate Ubuntu window or inside the PowerShell
window you are already in. Both are normal; answer wherever you are asked.

If nothing prompts you, start Ubuntu yourself: press Start and run *Ubuntu*, or
type `wsl` in PowerShell.

**5. From here on, use Ubuntu, not PowerShell.**

| Terminal | Its prompt | Used for |
|---|---|---|
| Ubuntu | `you@MACHINE:~$` | every command in this guide |
| PowerShell | `PS C:\Users\you>` | the `wsl` commands above, and section 20's WSL fixes |

Go by the prompt, not the window. Typing `wsl` in PowerShell starts Ubuntu
**inside that same window**, so the title bar still says PowerShell while
everything you type now goes to Linux. Type `exit` to come back — you will
need to, because the WSL2 fixes in section 20 are PowerShell commands and
cannot run at an Ubuntu prompt.

**6. Install the two tools Ubuntu does not always ship with:**

```bash
sudo apt update && sudo apt install -y curl git
```

It asks for the Linux password from step 4, and is harmless if they are already
installed.

If this reports `Could not resolve host`, Ubuntu has no working DNS — see
section 20, *Ubuntu cannot download anything*. Everything from here on
downloads something, so fix it before continuing.

> **Where to keep your data.** Work inside the Linux home directory — where the
> Ubuntu terminal starts. Windows drives are reachable under `/mnt/c/`, but
> reading across that boundary is several times slower and is a common reason a
> run takes far longer than it should. Copy data to the Linux side first;
> section 12 shows how.

WSL has two failure modes worth knowing about before a long run — a clock that
drifts from Windows and stops a run near the end, and Windows line endings
breaking scripts. Both are in section 20 with their fixes.

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
terminals started afterwards. On Windows that means the **Ubuntu** window, not
PowerShell — conda was installed inside Ubuntu and exists nowhere else.

**3. Check it worked:**

```bash
conda --version
```

You should see something like `conda 24.x.x`.

> **`conda: command not found`**
> The shell has not picked up the install. Close the terminal and open a new
> one. If it still fails, run `source ~/.bashrc` (or `source ~/.zshrc` on
> macOS) and try again.
>
> On Windows, check the prompt first. Both of those are Ubuntu commands, and
> in PowerShell they fail too — so if the prompt reads `PS C:\Users\you>` you
> are simply in the wrong window, and nothing is wrong with the install.

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
| `conda` | open a new terminal — the Ubuntu one on Windows; if it persists, `source ~/.bashrc` |
| `git` | `sudo apt install -y git` (Linux, or Ubuntu on Windows), or `xcode-select --install` (macOS) |

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

---

## 7. Try it on real data

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
section 8 describes. `fastq_fail` and `unclassified` reads are not included.

Start with `Flongle_Demo01`: it is the smallest, and takes about 45 minutes on
a 20-core machine — see the table in section 10 for the others.

```bash
mkdir -p ~/data && cd ~/data
curl -L -O https://zenodo.org/records/21998287/files/Flongle_Demo01.zip
unzip Flongle_Demo01.zip

nano16s -d ~/data/Flongle_Demo01/fastq_pass -o ~/nano16s_out/Flongle_Demo01
```

Leave room for the download, the unpacked reads, and 2–3× the unpacked size
for the run itself — so about 2 GB all told for `Flongle_Demo01`, and around
25 GB for `PromethION_Demo01`. The whole set is about 12 GB downloaded.

`CHECKSUMS.txt` in the record lets you confirm a download arrived intact:

```bash
curl -L -O https://zenodo.org/records/21998287/files/CHECKSUMS.txt
sha256sum -c CHECKSUMS.txt --ignore-missing
```

The DOI always resolves to the current version. The data is licensed CC BY 4.0
— cite it if you use it in your own work (see section 21).

---

## 8. Prepare your data

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
is several times slower. Section 12, *Opening your results*, shows how to move
data between Windows and Ubuntu in both directions.

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

Section 16 shows how to apply it when loading the tables. Doing this at the
start rather than at analysis time saves real confusion — barcode numbering and
sample numbering rarely line up.

**Barcode numbers repeat between runs.** A barcode number identifies a sample
only within one sequencing run. Kits reuse the same barcodes, so `barcode07`
from one run and `barcode07` from the next are different samples — in one
dataset, eight separate positive controls were all `barcode07`, each on its own
flow cell. nano16s never reads the barcode sequence; MinKNOW has already sorted
the reads by then, and a sample is simply whatever its directory is called. So:

- Keep each run's `fastq_pass/` separate, or use `nano16s batch`, which gives
  every run its own output directory. Never copy two runs' barcode directories
  into one folder: the second `barcode07` lands on top of the first.
- To analyse samples from several runs together, give each directory a unique
  name first. Any name that starts with `barcode` and uses only letters,
  digits, dot, dash and underscore works, and the name is carried into every
  table and report:

```bash
mkdir combined
cp -r run_A/fastq_pass/barcode07 combined/barcodeCatfish_pool5
cp -r run_B/fastq_pass/barcode07 combined/barcodeCatfish_pool6
```

- If you are unsure where a file came from, its reads say. Every MinKNOW read
  header records the barcode, flow cell, sample ID and start time:

```bash
zcat barcode07/*.fastq.gz | head -1 | tr ' ' '\n' \
    | grep -E '^(barcode|flow_cell_id|sample_id|start_time)='
```

---

## 9. Choose settings for your amplicon

Three settings are worth a decision before a real run: the length window, the
quality threshold, and how many cores to give it. The defaults suit
full-length 16S; the first is the one that will ruin a run if it is wrong.

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
nano16s -d /path/to/your_run/fastq_pass -o ~/length_check \
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

### Read depth — the setting that shortens a long run

`--max-reads` classifies at most that many reads per barcode. It is off by
default, and every read that passes the filter is classified.

```bash
nano16s -d /path/to/your_run/fastq_pass --max-reads 25000
```

It exists because classification is where a large run spends its time, and
that cost scales with the number of reads. Trimming settings and core counts
do not change it. On a run of several million reads, this is the only setting
that makes a real difference.

**It changes your results, so treat it as a decision about your samples rather
than a speed knob.** Community profiles are usually stable well below full
depth, but "usually" is not "yours". Before using it across a run, check it on
one barcode:

```bash
# the same barcode, both ways, into two output directories
nano16s -d /path/to/one_barcode_only -o full_depth
nano16s -d /path/to/one_barcode_only -o reduced --max-reads 25000
```

Compare `07_emu_combined/emu-combined-species.tsv` between the two. If the taxa
you care about hold their abundances, the reduced depth is enough for your
samples.

What that comparison looked like on one barcode of `MinION_Demo01`, at
142,331 reads against 25,000:

| | full depth | 25,000 reads |
|---|---|---|
| classification | 32 min | 5.6 min |
| top five species | — | the same five, same order |
| taxa above 1% | 15 | all 15 still found |
| largest change among them | — | 0.31 percentage points |

The most abundant organism moved from 35.01% to 34.70%. The count of taxa
detected at all fell from 247 to 151, which is the rare tail: organisms seen a
handful of times in 142,000 reads, where a single read is the difference
between present and absent.

Depth is what decides this, not the fraction you keep. The bundled demo data
holds only about a thousand reads a barcode, and cutting that to 400 reordered
the top five species in five of six barcodes — there was too little signal to
begin with. At the depths that make a run slow there is far more room.

The same reads are chosen every time for a given input, so a result can be
reproduced. `subsample_seed` in `config/config.yaml` changes the draw, which is
worth doing once to check a result is not an artefact of one particular sample.

Reads are drawn from across the whole barcode, not from the start of it —
nanopore writes reads in the order the pores produced them and pore quality
drifts over a run, so the first N reads are the run's beginning rather than a
sample of it.

---

## 10. Plan the run (time, disk, memory)

### How long

Runtime scales with the number of reads, not the number of barcodes. As
measured on a 20-core workstation:

| Reads in the run | Barcodes | Elapsed | Dataset (section 7) |
|---|---|---|---|
| ~200,000 | 24 | ~30 min | `Flongle_Demo01` |
| ~340,000 | 24 | ~55 min | `Flongle_Demo03` |
| ~510,000 | 16 | ~65 min | `Flongle_Demo02` |
| ~2,300,000 | 24 | ~3h 35m | `MinION_Demo01` |
| ~3,200,000 | 24 | ~7h | `PromethION_Demo01` |

Measured in one batch on a 20-core workstation, so the five are comparable
with each other.

Treat these as a rough guide — a machine with a quarter of the cores takes
substantially longer. Two stages dominate: **Porechop**, which infers adapter
sequences from your data rather than assuming them, and **Emu**, which does the
classification.

Their balance shifts with the size of the run, which is why the large datasets
above are not simply scaled-up versions of the small ones. Porechop's share
falls as the read count rises and Emu's grows, so a run of a few hundred
thousand reads spends most of its time trimming, while one of several million
spends most of it classifying.

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

## 11. Run it

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
filter settings, cores — and begins. The one time it stops to ask is when free
disk looks insufficient for the run; `-y` answers that in advance, which is
what you want when running unattended.

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

**More than one run to process?** `nano16s batch` takes a directory of them and
does the lot in one command — section 19.

---

## 12. What you get

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
    ├── emu-combined-phylum-counts.tsv
    ├── read_accounting.tsv              where every read went, per barcode
    └── per_barcode_taxa.tsv             which species, per barcode

with --per-read, one more directory:
└── 08_per_read/                  one line per sequencing read
    └── <barcode>_per_read.tsv.gz
```

**Relative abundance vs counts.** Abundance tables give each taxon's proportion
of the sample, summing to 1 per column. Counts tables give Emu's estimated
number of reads. Use abundances to compare composition between samples; use
counts for methods that expect count data, such as differential-abundance
testing.

**Reads that were not classified are in the tables too**, on a row labelled
`Unclassified`. Emu writes that row without a name; nano16s labels it, so a
column still sums to the reads the classifier was given and you can drop it
deliberately rather than by accident.

**`read_accounting.tsv`** is the same arithmetic laid out per barcode:

| column | meaning |
|---|---|
| `raw_reads` | what came off the sequencer |
| `removed_by_filter` | removed by the length and quality filter |
| `filtered_reads` | what survived it |
| `subsampled_out` | left out by `--max-reads`, 0 when it is off |
| `reads_to_classifier` | what Emu was actually given |
| `reads_classified` | placed on a species |
| `reads_unclassified` | Emu could not place |
| `species_found`, `genera_found` | how many distinct taxa that barcode produced |
| `check` | `ok` when classified + unclassified equals what the classifier was given |

So `raw_reads = removed_by_filter + filtered_reads`, and `filtered_reads =
subsampled_out + reads_classified + reads_unclassified`. The report shows the
same table under **Every read accounted for**. If a species-level total ever
looks smaller than you expect, this is the first place to look: the reads are
either filtered out, not classified, or on a taxon you filtered away.

**`per_barcode_taxa.tsv`** answers the other half: not how many species a
barcode found, but which. One row per barcode and species, with the reads on it
and its share of that barcode's classified reads:

```
barcode     species              genus        reads   pct_of_classified
barcode01   Aeromonas veronii    Aeromonas    607     61.6434%
barcode01   Hafnia paralvei      Hafnia       243     24.6958%
barcode01   Ewingella americana  Ewingella    28      2.8729%
```

The combined tables hold the same numbers as a grid of taxa against barcodes,
which is the shape a heatmap wants. This is the shape a person wants: filter to
one barcode, or sort by reads, or keep everything above 1%. The reads column
sums exactly to that barcode's `reads_classified`, so the two tables agree to
the read.

**Already have results from an earlier version?** You do not need to re-run the
analysis. Point nano16s at the same output directory with the same command and
it writes the two new tables and refreshes the report — the reads are not
re-trimmed and not re-classified:

```bash
nano16s -d /path/to/fastq_pass -o my_results -y
```

To relabel the combined tables as well, so the unclassified row is named,
delete them first and run the same command. They are rebuilt from the
per-barcode results already on disk:

```bash
rm my_results/07_emu_combined/emu-combined-*.tsv
nano16s -d /path/to/fastq_pass -o my_results -y
```

On a six-barcode run that took six seconds, and the counts came out identical
to the original run. For `nano16s batch`, use the same command with `batch` and
the directory of runs.

### One line per read: `--per-read`

The tables above summarise a barcode. They cannot say *which* read supported a
call, so a read cannot be traced back, pulled out for a second opinion, or
counted by hand. Run with `--per-read` and each barcode also gets
`08_per_read/<barcode>_per_read.tsv.gz`:

```
read_id                               barcode    status  taxid   species            confidence  candidates  lineage
a2f7a15c-8ee9-4800-9dbc-5268366ec261  barcode01  C       654     Aeromonas veronii  1.0000      1           Bacteria|Pseudomonadota|...|Aeromonas veronii
b7098d5e-5f8a-40cc-a1ac-60a9ec97d0d3  barcode01  C       324617  Aeromonas tecta    0.9972      4           Bacteria|Pseudomonadota|...|Aeromonas tecta
```

| column | meaning |
|---|---|
| `status` | `C` confident, `A` ambiguous, `U` matched nothing |
| `confidence` | how much of that read's probability sits on the reported taxon |
| `candidates` | how many taxa the read matched at all |
| `lineage` | full lineage, domain first, pipe-separated |

**Why a confidence, and not just a name.** Emu does not label reads. It spreads
each read across the references it matched and estimates abundances from the
whole distribution. This file reports the taxon holding most of a read's
probability, with that probability beside it — so `1.0000` is one clear match,
while `0.55` means the read fits two references nearly equally and the species
name is close to a coin toss. Reads below 0.9 are marked `A` rather than
presented as decided. Every read the classifier saw has a line, including ones
that matched nothing, so the line count equals `reads_to_classifier` in
`read_accounting.tsv`.

Pull out one organism's reads:

```bash
zcat 08_per_read/barcode01_per_read.tsv.gz | awk -F'\t' '$5 == "Aeromonas veronii"'
```

**Cost.** No measurable extra time: Emu does the same work and writes one more
file. That file is reads × taxa and can reach gigabytes for one deep barcode,
so nano16s converts it to the compact table above and deletes it. The result is
roughly 30 bytes per read. It is off by default for one reason: Emu has to be
asked for the distribution while it classifies, so a run that already finished
without `--per-read` has to classify again to produce it.

### Opening your results

Both reports are single self-contained HTML files — no internet needed to view
them, safe to email or attach to a manuscript. Open one the way you would open
any web page, or open the whole output directory in your file manager:

| Where you are | Open the report | Open the folder |
|---|---|---|
| Linux desktop | `xdg-open nano16s_report.html` | `xdg-open .` |
| macOS | `open nano16s_report.html` | `open .` |
| Windows (WSL) | `explorer.exe nano16s_report.html` | `explorer.exe .` |
| No desktop (SSH, server) | copy it to your own machine — see below | |

The `.` means *here*, so run it from inside your output directory; `pwd` prints
the full path if you need it. On Ubuntu, `open` is an alias for `xdg-open`. On
Linux, if nothing happens, name the file manager directly — `nautilus .`,
`dolphin .` or `thunar .`.

**Do not put results under `/tmp`.** On Ubuntu the default Firefox is a
**snap**, and a snap gets its own private `/tmp` — so the browser genuinely
cannot see a file that `ls` shows you in the terminal, and reports *File not
found* for a path that exists. `/tmp` is also cleared on reboot.

`nano16s test` writes into your home directory for this reason, and `-o`
should point somewhere under it too. If you already have a report stuck under
`/tmp`, copy it out:

```bash
cp /tmp/nano16s_test_*/nano16s_report.html ~/
xdg-open ~/nano16s_report.html
```

#### On Windows

Ubuntu's files are not on your `C:` drive. Windows reaches them over a
network-style path, which is normal — they are still on your own machine:

```
\\wsl.localhost\Ubuntu\home\<your-linux-username>
```

Paste that into the Explorer address bar, or use the **Linux** entry at the
bottom of the Explorer sidebar. Right-click your results folder there and
choose *Pin to Quick access* to keep it one click away. `whoami` in Ubuntu
prints your Linux username if you have forgotten it.

Your Windows drives work the other way round, which is how you copy sequencing
data in:

```bash
cp -r /mnt/c/Users/YourWindowsName/Desktop/fastq_pass ~/
```

Copy rather than running from `/mnt/c/` directly — reading across that boundary
is several times slower.

> Browse and copy through Explorer freely, but **edit files inside Ubuntu**.
> Saving into `\\wsl.localhost\...` from a Windows editor can change line
> endings and strip the executable bit, which breaks scripts — the `$'{\r'`
> error in section 20.

#### On a server with no desktop

There is no window to open, so copy the results to your own machine. Run this
on your own computer, not the server:

```bash
scp you@server:/home/you/my_results/*.html ~/Desktop/
scp -r you@server:/home/you/my_results/07_emu_combined ~/Desktop/
```

That is the reports and the tables — a few megabytes. The rest of the output
directory is intermediate FASTQ and rarely worth moving.

---

## 13. Read the results report

Open `nano16s_report.html`. This is the biology — what was found in each
sample, and how much of it.

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

## 14. Read the performance report

Open `performance_report.html`. This one is about the run rather than the
biology — reach for it when something took longer than expected, or a barcode
looks wrong.

**Headline figures.** Elapsed time, total CPU time, average number of jobs
running at once, percentage of your cores used, and peak memory.

CPU time, core use and peak memory read `-` on macOS. Snakemake collects them
by sampling each job's process tree, which macOS does not permit, so it
records none of them there; the report says so beneath the stage table.
Elapsed time and the per-stage wall times are measured directly and are
correct on every platform.

Elapsed counts the time the machine was **working**, which for an
uninterrupted run is simply start to finish. If you resumed a run, only the
stages that actually re-ran wrote new timings, so the report would otherwise
measure the calendar gap between your two sittings rather than either run. It
excludes that idle time instead, and says so above the tables — naming the
dates and how much was excluded. Core use and parallelism are computed on the
same basis.

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
Gaps within a run are idle capacity — cores that were free while something
else finished. A dashed vertical line marks a point where the run stopped and
was started again later; the wait between sittings is cut out rather than
drawn, so the chart stays readable.

Stages that take under a second — merging, and both NanoStat passes — appear
as thin ticks. That is honest scale against an axis measured in minutes, not
a rendering fault.

**No timing data.** If the whole report is blank of timings, the run had
nothing to do: every output was already current, so no job ran and no job
recorded a benchmark. The report says so rather than leaving you to guess.
Point `-o` at a new directory, or add `--forceall`, to get timings.

**Worth checking.** QC flags — see the next section.

**Per barcode.** Timings first, then read counts and quality. Click any column
heading to sort.

---

## 15. Quality control: what gets flagged

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

## 16. Use the tables downstream

The combined tables are plain tab-separated text. Rows are taxa. The first
columns are the taxonomic lineage — how many depends on the rank — and the rest
are one column per barcode.

**Match columns by name, not by position.** Barcode columns are not in sorted
order, and the lineage columns before them differ by rank — seven in the
species tables, six in genus, two in phylum. Every example below matches on
name, which is why none of them breaks when you switch rank.

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
Emu writes the lineage before the sample columns, so split the frame at the
first column whose name begins `barcode` — position 8 in the species tables,
7 in genus, 3 in phylum. Find it rather than hard-coding it, so the same script
works at every rank.

### Excel

Open the `.tsv` directly, or import as tab-delimited. Watch out for Excel
converting taxon names that look like dates.

---

## 17. Interpreting 16S results

**Genus is the defensible resolution.** Species-level assignment from 16S is
genuinely hard — many genera contain species whose 16S genes are near-identical,
and nanopore error rates make it harder. Species tables are provided because
they are sometimes informative, but conclusions are safer at genus level.

**Running the same data twice can move species-level numbers.** This is worth
knowing before you put a species table in a paper. Two from-scratch runs of the
same MinION dataset, same database, same settings, gave:

| Rank | Largest difference between the two runs |
|---|---|
| Species | **17.5 percentage points** |
| Genus | 0.43 pp |
| Phylum | 0.04 pp |

Raw read counts were identical in all 24 barcodes. Four barcodes ended with
slightly different *filtered* counts — 31 to 66 reads out of 80,000 to 120,000,
under 0.1% — because Porechop infers adapter sequences from the reads
themselves and does not always infer exactly the same ones. Every large shift
was in a barcode affected by that.

The shifts were entirely within one genus: *Lactobacillus gasseri*,
*L. johnsonii* and *L. taiwanensis* traded abundance with each other while the
*Lactobacillus* total moved by 0.02 pp. Those species have near-identical 16S
genes, so apportioning reads between them is an ill-conditioned problem — a
0.05% change in the input moves the answer a lot, and neither answer is more
correct than the other.

The other four demo datasets reproduced exactly, so this is occasional rather
than routine. Report genus-level abundances; if a species-level claim matters,
state that it sits inside a closely related complex and check it holds across
repeat runs.

**Some species are named wrongly, and consistently.** On the ZymoBIOMICS mock
community, whose species are known, nano16s with its default database got
every genus right but named four species as a close relative: *Escherichia
coli* as *E. fergusonii*, *Staphylococcus aureus* as *S. roterodami*, *Listeria
monocytogenes* as *L. cossartiae*, and the mock's *Bacillus spizizenii*
(sold as *B. subtilis*) as *B. rugosus*. Unlike the shifts above, these calls
were the same in every sample and every run, so repeating a run will not reveal
them. Each pair has near-identical 16S genes, and the default database holds
one to ten sequences per species, so a strain can sit closer to a relative's
single reference than to its own. The misnaming follows the database: the same
reads named against a database with hundreds of sequences per species came out
right. Read a species name inside such a group as "a member of this group", not
as an identification.

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

## 18. Re-running and changing settings

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

## 19. Processing several runs

Point `nano16s batch` at a directory holding several runs and it processes
every one of them, into its own subdirectory of the output:

```bash
nano16s batch -d ~/data -o ~/results
```

If `~/data` holds `Flongle_Demo01/`, `Flongle_Demo02/` and `MinION_Demo01/`,
that is the whole job — three runs, three sets of tables and reports, one
command. It prints what it found before starting:

```
nano16s 1.2.0 — batch of 3 run(s)
  input     /home/you/data
  output    /home/you/results
  runs      Flongle_Demo01 Flongle_Demo02 MinION_Demo01

=== [1/3] Flongle_Demo01 ===
    done
=== [2/3] Flongle_Demo02 ===
    done
=== [3/3] MinION_Demo01 ===
    done

Batch finished: 3 of 3 succeeded.
  results   /home/you/results/<run>/
  reports   /home/you/results/reports  (6 files)
```

**What counts as a run.** Any subdirectory of `-d` holding either a
`fastq_pass/` directory — what MinKNOW writes — or `barcode*` directories
directly. Anything else in there is ignored, so a stray `notes.txt` or an
old analysis folder does no harm.

**Where things land.**

```
~/results/
├── Flongle_Demo01/              a complete run directory, exactly as a single run
├── Flongle_Demo01.log           everything that run printed
├── Flongle_Demo02/
├── Flongle_Demo02.log
├── MinION_Demo01/
├── MinION_Demo01.log
└── reports/                     every report, named by run
    ├── Flongle_Demo01_report.html
    ├── Flongle_Demo01_performance_report.html
    └── ...
```

The `reports/` folder is the one to share. The reports are self-contained
single files, so it can be zipped and emailed as it is — five runs come to
well under a megabyte.

**Settings apply to the whole batch.** Every option a single run takes is
passed through unchanged:

```bash
nano16s batch -d ~/data -o ~/results --min-quality 12 --min-length 1300 -c 8
```

This is the point of the batch command for comparison work: one set of
settings, applied identically, with no chance of a run drifting. Preview the
whole thing first with `-n`, which lists the steps for every run and stops.

**One failure does not stop the rest.** A run that fails is reported, its log
named, and the batch carries on:

```
=== [2/3] Flongle_Demo02 ===
    FAILED — see /home/you/results/Flongle_Demo02.log

Batch finished: 2 of 3 succeeded.
  failed    Flongle_Demo02
  logs      /home/you/results/<run>.log
```

The command exits non-zero if anything failed, so it can be used in a script.
Because completed work is skipped, running the same batch again costs only the
runs that did not finish — fix the cause and repeat the command.

### Leave it running overnight

A batch of large runs takes hours. `nohup` keeps it going after the terminal
closes:

```bash
nohup nano16s batch -d ~/data -o ~/results > ~/batch.log 2>&1 &

tail -f ~/batch.log     # Ctrl-C stops watching, not the batch
```

The batch never prompts, so nothing can stall waiting for an answer.

### Comparing runs

**Keep settings identical across runs you intend to compare.** Running the
whole set through one `nano16s batch` command is the simplest way to guarantee
it. If you compare runs done separately, check the Methods paragraph of each
report — it records the length window, quality threshold and database version.

The combined tables stay per-run: nano16s does not merge samples from different
sequencing runs into one table, because barcode names collide — every run has a
`barcode01`. To analyse across runs, load each run's table separately and apply
that run's barcode map, as in section 16. That is the point at which barcode
names become sample names and the collision disappears.

---

## 20. Troubleshooting

Grouped by where the problem appears. Each entry is the message you will see,
so searching this page for a phrase from your error is the quickest way in.

### Setup

**`nano16s: command not found`**

Run `conda activate nano16s`. Needed in every new terminal.

**`no Emu database found`**

Run `nano16s db build` (about ten minutes, once).

### Your data

**`no barcode* directories inside ...`**

`-d` is one level too high or too low. It wants the directory that directly
contains `barcode01/`. Find it with
`find /path/to/your_run -type d -name fastq_pass`.

**`no .fastq.gz files in ...`**

That barcode directory is empty, or holds uncompressed `.fastq`. Compress them,
or remove the directory if the barcode genuinely produced nothing.

### During the run

**The run stopped partway**

Run the same command again. Completed steps are skipped.

**Porechop is slow**

Expected — it is the slowest stage by a wide margin because it infers adapter
sequences from your data instead of assuming them. Budget a few minutes per
barcode.

**The run is using less of my CPU than expected**

See the Machine use section of the performance report, and section 14.

**A barcode kept almost no reads**

Its reads fall outside your length window. Check the median read length in the
report and widen `--min-length` / `--max-length` if that length is expected for
your amplicon.

### Opening the reports

**The browser says *File not found* for a report that `ls` shows is there**

Almost always a report under `/tmp` opened with Ubuntu's snap Firefox. A snap
runs with its own private `/tmp`, so the file you can see in the terminal is
genuinely not in the `/tmp` the browser sees. `nano16s test` no longer writes
there, so this means either an older version or an `-o` you chose. Confirm
the browser with:

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

This group uses both windows, so check which one each command wants: a
`powershell` block is PowerShell, a `bash` block is Ubuntu. If a command is
not recognised, being in the wrong window is the likeliest reason before
anything else. Type `exit` to leave Ubuntu, or open PowerShell from the Start
menu.

**`wsl --install` sits at 0%**

Blocked from the Microsoft Store, not installing slowly. Use the two-step
`--web-download` commands in section 3.

**`Wsl/WSL_E_WSL_OPTIONAL_COMPONENT_REQUIRED`**

Expected between installing WSL and restarting, even when the install reported
success. Restart Windows and run the command again.

**Ubuntu cannot download anything**

```
curl: (6) Could not resolve host: github.com
```

Ubuntu is running but cannot turn a name into an address, so `apt`, `git clone`,
Miniforge and the database build all fail. Common on managed work networks.

First find out whether it is only name lookup that is broken, or the network
itself. These two need different fixes and give the same error:

```bash
getent hosts github.com                                    # can it look up a name?
timeout 5 bash -c 'exec 3<>/dev/tcp/1.1.1.1/80' && echo OK  # can it connect at all?
```

If the first prints nothing but the second prints `OK`, only DNS is broken and
the fixes below apply. If neither works, the network itself is blocked — an IT
request, not something WSL can be configured around.

**Fix — let WSL use the Windows network stack.** In **PowerShell**, as one
command:

```powershell
Set-Content -Path "$env:USERPROFILE\.wslconfig" -Encoding ascii -Value '[wsl2]','networkingMode=mirrored','dnsTunneling=true'
```

Then `wsl --shutdown`, reopen Ubuntu and retry. Needs Windows build 22621 or
higher — `wsl --version` reports the build on its last line. Windows 10 is
below it, and says so on the next launch:

```
wsl: Mirrored networking mode is not supported: Windows version 19045.6456
does not have the required features.
Falling back to NAT networking.
```

If you see that, this route is closed on your machine — go straight to the
next fix, which does work there.

> Use that command as written rather than typing the file by hand. `.wslconfig`
> needs three separate lines, and if it ends up on one, WSL reports
> `Expected ' ' or '\n' in ...\.wslconfig:1` at every launch and **ignores the
> file**, so the fix appears not to work when it was never applied. If you see
> that message, delete the file with
> `Remove-Item "$env:USERPROFILE\.wslconfig"` and run the command above.

**If that is unavailable or does not help**, point Ubuntu at the DNS servers
Windows itself uses. Public resolvers such as `8.8.8.8` are often blocked on
managed networks, which is why these are the ones to copy. In **PowerShell**:

```powershell
Get-DnsClientServerAddress -AddressFamily IPv4 | Where-Object ServerAddresses
```

Then in **Ubuntu**, using the addresses for the adapter you are connected
through:

```bash
grep -q generateResolvConf /etc/wsl.conf 2>/dev/null || sudo tee -a /etc/wsl.conf >/dev/null <<'END'

[network]
generateResolvConf = false
END
sudo rm -f /etc/resolv.conf
printf 'nameserver 10.0.0.1\nnameserver 10.0.0.2\n' | sudo tee /etc/resolv.conf
```

Replace those two addresses with the ones PowerShell printed.

> The `grep` guard matters if you try this more than once. The command appends,
> so without it a second attempt adds a second `[network]` section and WSL then
> reports at every launch:
>
> ```
> wsl: Duplicated config key 'network.generateResolvConf' in /etc/wsl.conf:11
> ```
>
> That is a warning, not a failure — WSL uses the first one and DNS still
> works. To clear it, open the file with `sudo nano /etc/wsl.conf` and delete
> the repeated `[network]` sections until only one remains.

Run `wsl --shutdown` in PowerShell, reopen Ubuntu and retry.

A VPN can also take DNS over in a way WSL2 does not follow — disconnect and
retry before assuming anything else.

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

---

## 21. Reference

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
| `--per-read` | off | also write one line per read: what it was called, how sure, full lineage |
| `-y, --yes` | | skip confirmation prompts |
| `-h, --help` | | full help |

Subcommands: `nano16s db build`, `nano16s db list`, `nano16s test`,
`nano16s --version`.

`nano16s batch -d <dir-of-runs> -o <dir>` processes every run under one
directory, into `<dir>/<run-name>/`, passing all the options above through
to each. See section 19.

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
performance report shows low core utilisation, lowering these is the lever. On
macOS core use is not measured, so tune against elapsed time instead.

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

If you use the demo datasets from section 7, cite them as well:

- Zhang, L. & Adapa, P. D. *nano16s demo datasets: Oxford Nanopore full-length
  16S rRNA amplicon sequencing runs (Flongle, MinION, PromethION)*. Zenodo (2026).
  <https://doi.org/10.5281/zenodo.21998286>
