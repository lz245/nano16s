# =============================================================================
# common.smk — configuration access and sample discovery
# =============================================================================
# Included by every top-level Snakefile. Defines no rules.
# =============================================================================

import glob
import os

INPUT_DIR = config["input_dir"]
OUTPUT_DIR = config["output_dir"]

# Auto-detect samples from barcode directories, or use an explicit list.
if config.get("samples"):
    SAMPLES = config["samples"]
else:
    barcode_dirs = sorted(glob.glob(os.path.join(INPUT_DIR, "barcode*")))
    SAMPLES = [os.path.basename(d) for d in barcode_dirs if os.path.isdir(d)]

if not SAMPLES:
    raise ValueError(f"No barcode directories found in {INPUT_DIR}")

EMU_RANKS = config.get("emu_ranks", ["species", "genus", "phylum"])
MINIMAP_RANKS = config.get("minimap2_ranks", ["species", "genus", "phylum"])
