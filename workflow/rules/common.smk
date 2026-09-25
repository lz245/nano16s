# =============================================================================
# common.smk — configuration access and sample discovery
# =============================================================================
# Included by every top-level Snakefile. Defines no rules.
# =============================================================================

import glob
import os
import re

INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]

# A sample name reaches the shell inside every rule and becomes part of every
# output path. Anything outside this set -- a space, a bracket, a quote --
# breaks the command that uses it. `barcode02 (copy)` is the case that found
# this: duplicating a folder in Finder or Explorer produces exactly that name,
# it was picked up as a sample, and the run died in `merge` on the unquoted
# path with an error that quoted the rule's own comment text rather than
# naming the directory.
SAMPLE_RE = re.compile(r"^[A-Za-z0-9._-]+$")

# Also applied to the wildcard itself, so a path Snakemake parses back into a
# sample cannot smuggle in a name the check above rejected.
wildcard_constraints:
    sample = r"[A-Za-z0-9._-]+",
    rank   = r"[A-Za-z]+",


def _check_names(names):
    """Reject names that cannot survive a shell command, naming each one."""
    bad = [n for n in names if not SAMPLE_RE.match(n)]
    if bad:
        listed = "\n".join(f"    {n}" for n in sorted(bad))
        raise ValueError(
            "These barcode directory names contain characters that cannot be "
            "used in a sample name:\n"
            f"{listed}\n"
            "  Letters, digits, dot, dash and underscore only.\n"
            "  Rename them, or move them out of the input directory if they "
            "are not samples,\n"
            "  then run the same command again."
        )
    return names


# Auto-detect samples from barcode directories, or use an explicit list.
if config.get("samples"):
    SAMPLES = _check_names(config["samples"])
else:
    barcode_dirs = sorted(glob.glob(os.path.join(INPUT_DIR, "barcode*")))
    SAMPLES = _check_names(
        [os.path.basename(d) for d in barcode_dirs if os.path.isdir(d)]
    )

if not SAMPLES:
    raise ValueError(f"No barcode directories found in {INPUT_DIR}")

EMU_RANKS = config.get("emu_ranks", ["species", "genus", "phylum"])


def flag(name, default=False):
    """A config value that means yes or no, however it arrived.

    From config.yaml it is a real boolean; from `--config per_read=false` on
    the command line it is the *string* "false", and every non-empty string is
    true in Python. Read literally, an option switched off switched itself on.
    """
    v = config.get(name, default)
    if isinstance(v, str):
        return v.strip().lower() in {"1", "true", "yes", "on"}
    return bool(v)


PER_READ = flag("per_read")
MINIMAP_RANKS = config.get("minimap2_ranks", ["species", "genus", "phylum"])
