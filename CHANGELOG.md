# Changelog

Notable changes to nano16s. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Reference database versions are separate from the tool version. A database is
identified by its build month (for example `2026.08`) and recorded in the run
report, so a result can always be traced to the database that produced it.

## [Unreleased]

### Fixed
- The performance report now says when it has no timings, instead of quietly
  dropping the sections that would have held them. A run into an output
  directory that already holds finished results gives Snakemake nothing to
  run, so no job records a benchmark; the report still rendered its read
  counts and quality — which come from NanoStat — while the timing sections
  vanished and every timing column read `-`, with nothing to explain why.
  Directories created before 1.1.0 hit this the first time they are re-used.
- A barcode whose reads were all filtered out is no longer reported as having
  a median quality of Q0.0. NanoStat writes zeros for an empty file rather
  than omitting the fields, so an emptied barcode read as a measurement of
  zero and drew a second flag beside the retention one that already gave the
  real reason. A run whose length window does not match the amplicon puts
  every barcode into this state, so it doubled the flags exactly when the
  report most needed to be readable.
- On macOS, a system version file that will not open no longer takes the whole
  performance report down. Reading the OS name goes through
  `platform.mac_ver()`, which parses a plist under `/System/Library`; every
  other hardware lookup in the report already tolerated its source being
  unreadable, and this one did not.

## [1.1.0] — 2026-08-11

### Added
- A performance report, written on every run alongside the existing one:
  `performance_report.html`, `performance_summary.csv` and `performance.json`.
  It reports where the run spent its time, how much of the machine it used,
  and which barcodes fall below the QC floors.
- `benchmark:` records for every per-sample rule, so timings carry CPU time and
  peak memory rather than wall clock alone. These land in `benchmarks/`.
- The performance report records the machine it ran on — CPU model, physical
  cores and threads, memory, operating system, kernel and architecture — since
  a runtime means little without the hardware behind it. WSL is reported as
  both Linux and Windows, because it is both and the distinction affects I/O.
  Detection uses only the standard library, so no dependency is added and a
  machine that will not answer leaves the field blank rather than failing the
  run.
- QC checks with absolute floors — retention, read depth, median quality, and
  whether the filtered median read length falls inside the configured length
  window. A purely relative check cannot fire when a whole run is uniformly
  bad, which is the case most worth catching.
- `CONTRIBUTING.md`, this changelog, `ruff.toml` and pre-commit hooks.
- Issue templates that ask for the version and environment details needed to
  reproduce a problem.
- `.gitattributes` pinning text files to LF, so a Windows checkout cannot break
  the shell scripts and the heredocs inside the Snakemake rules.

### Fixed
- The report's manuscript paragraph now cites nano16s. It previously credited
  Porechop_ABI, Chopper, NanoStat, Emu and Snakemake but not the pipeline
  itself, so a user pasting it into a Methods section cited five other tools and
  not this one.
- Tool versions no longer render as `Emu vv3.6.2`. Emu reports its version with
  a leading `v` where the other tools do not; that prefix is now stripped at
  capture, for every tool rather than as a special case.

## [1.0.0]

First release.

### Added
- Snakemake workflow: merge per barcode, QC with NanoStat, adapter trimming with
  Porechop_ABI, length and quality filtering with Chopper, QC again, then
  taxonomic classification with Emu.
- Combined abundance tables at species, genus and phylum level, as relative
  abundance and as read counts.
- Self-contained HTML report — one file, no external assets, safe to email.
  Records tool versions, parameters and database version.
- `nano16s db build`, which builds an Emu database from the current NCBI 16S
  RefSeq Targeted Loci collection. Emu's own database has not been refreshed
  since September 2020. Databases install side by side and never overwrite each
  other, so re-running an old analysis against its original database stays
  possible.
- `nano16s test`, an end-to-end run on bundled demo data that verifies an actual
  installation.
- Unit tests for the parsing functions, and CI covering Linux and Apple Silicon.
- One-command install via `install.sh`.
