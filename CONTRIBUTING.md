# Contributing to nano16s

Thanks for your interest. This is a small project, so the process is light.

## Reporting a problem

Open an issue. The most useful report includes:

- what you ran (the exact `nano16s` command)
- what happened, and what you expected
- your operating system, and whether you are using WSL2
- the output of:

```bash
nano16s --version
conda list -n nano16s | grep -E "emu|minimap2|chopper|porechop|snakemake"
```

If a run failed, the error text matters more than a description of it — please
paste it rather than summarising.

## Setting up for development

```bash
git clone https://github.com/lz245/nano16s.git
cd nano16s
bash install.sh
conda activate nano16s
```

Then build a reference database once, so you can run the pipeline end to end:

```bash
nano16s db build
```

## Testing

There are two levels, and they answer different questions.

**Unit tests** cover the parsing functions — reading Emu's tables, reading the
preprocessing summary, walking NCBI's taxonomy — plus the performance report,
the batch CLI, and a static check that no rule's shell body uses an array
expansion bash 3.2 rejects. They need only `pytest`, no bioinformatics tools
and no database, and finish in a few seconds:

```bash
python -m pip install pytest
python -m pytest test/ -q
```

**The demo run** exercises a real installation end to end on the bundled
six-barcode dataset. It takes a few minutes, most of it Porechop:

```bash
nano16s test
```

It must print `Install verified — everything works.`

Run the unit tests while developing; run `nano16s test` before opening a pull
request. Between them, the first catches parsing regressions in seconds and the
second catches everything else.

**What CI runs.** Both of the above, on every pull request: the unit tests on
Linux and macOS across the supported Python range, and the demo run on Linux,
database build included — seven to ten minutes, nearly all of it the
pipeline.
The scheduled weekly job runs the demo again on both platforms — that one is
about the dependencies rather than the code, since bioconda moves underneath us
and a tool release can break the environment without anything here changing.

So a pull request that breaks a rule fails before it merges, and a dependency
that breaks one fails within a week.

## Linting

```bash
python -m pip install ruff
ruff check .
```

The rule set in `ruff.toml` is deliberately narrow — pyflakes plus the
pycodestyle error subset. It is meant to catch defects, not to impose a style.

Optionally, run the checks automatically before each commit:

```bash
python -m pip install pre-commit
pre-commit install
```

## Submitting a change

1. Fork the repository and create a branch from `main`.
2. Make the change. Keep one topic per branch — separate pull requests are
   easier to review, and one can be accepted without waiting on another.
3. Run the unit tests, and `nano16s test` if you touched the pipeline or the
   report.
4. Open a pull request describing what changes and why. If it alters results,
   say so explicitly.

### A note on Snakemake rules

Editing the shell body of a rule makes Snakemake re-run every job that used it.
Changes to `emu`, `porechop` or `minimap2_align` are therefore expensive to
verify on a real dataset — a single deep barcode can take an hour. Downstream
rules such as `emu_combine` and `report` are cheap. Worth knowing before you
start a run to check a comment change.

### Line endings

`.gitattributes` pins text files to LF. The pipeline runs under bash, and a CRLF
checkout breaks every shell script and the heredocs inside the Snakemake rules
at once, with an error that gives no clue as to the cause. Please do not remove
it, and do not commit files with CRLF endings.

## Scope

nano16s is intentionally small: full-length 16S from Oxford Nanopore, on one
ordinary computer, with one install command. Changes that keep it easy to
install and easy to explain are more welcome than ones that add breadth.

If you are planning something substantial, open an issue first so we can talk
about it before you write the code.

## Releasing

The version appears in three places and they must agree, because each is read
by something different: `bin/nano16s` reports it to the user and stamps it into
both reports, `CITATION.cff` is what a citation resolves to, and
`conda-recipe/meta.yaml` decides which tarball a package build downloads.

```bash
grep -n 'VERSION=' bin/nano16s | head -1     # bin/nano16s
grep -n '^version' CITATION.cff              # CITATION.cff
grep -n 'set version' conda-recipe/meta.yaml # conda-recipe/meta.yaml
```

Bump all three, rename the changelog's `[Unreleased]` heading to the new
version with today's date and leave a fresh empty `[Unreleased]` above it, then
open the pull request. Once it merges, tag the merge commit and publish the
release from it.

Tag promptly after a fix that matters to users. A tag that sits behind `main`
is not merely stale: `nano16s db build` and the install instructions send
people to the release, and a conda package is built from the tagged tarball, so
everyone who does not clone `main` keeps getting the bug until the tag moves.

## Licence

By contributing you agree that your contributions are licensed under the MIT
Licence, the same as the rest of the project.
