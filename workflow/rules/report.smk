# =============================================================================
# report.smk — self-contained HTML summary
# =============================================================================
# Requires from common.smk: OUTPUT_DIR, SAMPLES, EMU_RANKS
# Requires from preprocess.smk: {OUTPUT_DIR}/preprocessing_summary.csv
# Requires from emu.smk: {OUTPUT_DIR}/07_emu_combined/emu-combined-{rank}.tsv
# =============================================================================


rule report:
    input:
        summary = f"{OUTPUT_DIR}/preprocessing_summary.csv",
        species = f"{OUTPUT_DIR}/07_emu_combined/emu-combined-species.tsv",
        genus   = f"{OUTPUT_DIR}/07_emu_combined/emu-combined-genus.tsv",
        phylum  = f"{OUTPUT_DIR}/07_emu_combined/emu-combined-phylum.tsv",
    output:
        f"{OUTPUT_DIR}/nano16s_report.html"
    params:
        db          = config["emu_db"],
        min_length  = config["min_length"],
        max_length  = config["max_length"],
        min_quality = config["min_quality"],
    script:
        "../scripts/make_report.py"
