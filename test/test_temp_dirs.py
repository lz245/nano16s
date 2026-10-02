"""No rule may build a temporary directory from a name two runs could share.

Porechop_ABI writes its k-mer counts to ./tmp relative to the current
directory, so each job is given a working directory of its own. That directory
was named `porechop_tmp_<barcode>` under TMPDIR -- which is the same path for
every job with that barcode anywhere on the machine.

Two nano16s runs at once is two flow cells, or a batch beside a single sample,
and both have a barcode01. The first job to finish ran its cleanup trap and
deleted the directory the second was still working in. Reproduced on demo data:
two simultaneous runs failed with six and eight failed porechop jobs, and
porechop reported "COULD NOT OPEN FILE ./tmp/temp_approx_kmer_count_sup_2.start"
-- a file in a directory that had been removed underneath it.

A static scan, so it holds for every rule and on every pull request.

Run with:
    python -m pytest test/ -v
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RULES = sorted((ROOT / "workflow" / "rules").glob("*.smk"))

SHELL_BODY = re.compile(r'shell:\s*\n\s*"""(.*?)"""', re.S)
# WORK="$TMPDIR/something" -- an assignment whose value is a literal path under
# a temporary root, rather than a directory mktemp has just created.
LITERAL_TMP = re.compile(r'^\s*(\w+)="\$(?:TMPDIR|\{\{?TMPDIR)[^"]*"', re.M)


def bodies():
    for f in RULES:
        for body in SHELL_BODY.findall(f.read_text()):
            yield f.name, body


def test_rules_are_present():
    """A scan that finds nothing must not pass quietly."""
    assert list(bodies()), "no rule shell bodies found to scan"


def test_no_temporary_directory_is_named_from_a_wildcard():
    """REGRESSION: `WORK="$TMPDIR/porechop_tmp_{wildcards.sample}"` is shared
    by every run on the machine that has that barcode."""
    offenders = []
    for name, body in bodies():
        for m in LITERAL_TMP.finditer(body):
            if m.group(1) == "TMPDIR":
                continue        # setting the root to fall back on, not a directory
            line = m.group(0).strip()
            if "mktemp" not in line:
                offenders.append(f"{name}: {line}")
    assert not offenders, (
        "a temporary directory is built from a literal name; use mktemp -d so "
        "two runs cannot collide:\n  " + "\n  ".join(offenders))


def test_every_working_directory_comes_from_mktemp():
    """The positive form: wherever a rule makes itself a working directory
    under TMPDIR, it is one mktemp created."""
    for name, body in bodies():
        if "TMPDIR" not in body:
            continue
        assert "mktemp -d" in body, (
            f"{name} uses TMPDIR without mktemp -d; two concurrent runs would "
            f"share the path")


def test_the_working_directory_is_still_cleaned_up():
    """mktemp leaves it behind on its own, so the trap has to stay."""
    for name, body in bodies():
        if "mktemp -d" in body:
            assert "trap" in body and "rm -rf" in body, (
                f"{name} creates a temporary directory and never removes it")
