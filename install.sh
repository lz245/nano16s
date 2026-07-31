#!/usr/bin/env bash
# =============================================================================
# nano16s — installer
# =============================================================================
# Creates the conda environment and puts `nano16s` on your PATH.
#
#   bash install.sh
#
# Unlike an alias, the command installed here works inside scripts, cron jobs,
# and any non-interactive shell.
# =============================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_NAME="nano16s"

say()  { printf '%s\n' "$*"; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

say "nano16s installer"
say "  source: $ROOT"
say ""

# --- 1. conda -----------------------------------------------------------------
if command -v mamba >/dev/null 2>&1; then
    SOLVER=mamba
elif command -v conda >/dev/null 2>&1; then
    SOLVER=conda
else
    die "conda not found.

Install Miniforge first (free, no licence restrictions):
  https://github.com/conda-forge/miniforge#install

On macOS you can also use Homebrew:
  brew install --cask miniforge"
fi
say "[1/3] Using $SOLVER to build the environment (a few minutes the first time)"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    say "      environment '$ENV_NAME' exists — updating"
    "$SOLVER" env update -n "$ENV_NAME" -f "$ROOT/environment.yml" --prune
else
    "$SOLVER" env create -f "$ROOT/environment.yml" -y
fi

# --- 2. put nano16s on PATH ----------------------------------------------------
ENV_PREFIX="$(conda env list | awk -v n="$ENV_NAME" '$1==n {print $NF}')"
[[ -n "$ENV_PREFIX" && -d "$ENV_PREFIX" ]] || die "could not locate the '$ENV_NAME' environment"

say "[2/3] Linking the nano16s command into $ENV_PREFIX/bin"
chmod +x "$ROOT/bin/nano16s"
ln -sf "$ROOT/bin/nano16s" "$ENV_PREFIX/bin/nano16s"

# --- 3. done -------------------------------------------------------------------
say "[3/3] Done."
say ""
say "Next steps:"
say ""
say "  conda activate $ENV_NAME"
say "  nano16s db build      # download and build the reference database (~10 min, once)"
say "  nano16s test          # verify the install on the bundled demo data (~2 min)"
say ""
say "Then analyse your own run:"
say ""
say "  nano16s -d /path/to/fastq_pass -o my_results"
say ""
