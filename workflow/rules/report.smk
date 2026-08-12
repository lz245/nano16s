# =============================================================================
# report.smk — self-contained HTML summaries
# =============================================================================
# Two reports, deliberately separate. `report` covers the biology and is what
# most readers open. `performance_report` covers the run — timings, machine
# use, QC thresholds — and is what you open when something took too long or a
# barcode looks wrong.
#
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
        # Set by the nano16s CLI. Empty when Snakemake is run directly, in which
        # case the report cites the tool by name without a version.
        version     = config.get("nano16s_version", ""),
    script:
        "../scripts/make_report.py"


def _total_cores(wildcards):
    """Cores this run was given.

    Reading workflow.cores raises when --cores is unset, which is exactly the
    case during `snakemake -n`. Snakemake evaluates params callables while
    building the DAG, so an unguarded read breaks dry runs. Any real run has
    --cores, so the number is there whenever it is meaningful; 0 means "not
    known" and the report says so rather than inventing a denominator.
    """
    try:
        return int(workflow.cores)
    except Exception:
        return 0


rule performance_report:
    # Depends on the finished per-sample artefacts rather than on the
    # benchmark files themselves. Benchmark paths are a side effect of a rule,
    # not a declared output, so asking Snakemake to build one as an input is
    # not something the DAG can resolve. Waiting on the real outputs gives the
    # same guarantee — every job has finished — and the script then reads the
    # benchmark directory directly.
    input:
        summary  = f"{OUTPUT_DIR}/preprocessing_summary.csv",
        raw      = expand(
            f"{OUTPUT_DIR}/02_nanostat_raw/{{sample}}/{{sample}}_quality_summary.txt",
            sample=SAMPLES,
        ),
        filtered = expand(
            f"{OUTPUT_DIR}/05_nanostat_filtered/{{sample}}/{{sample}}_filtered_quality_summary.txt",
            sample=SAMPLES,
        ),
        emu      = expand(
            f"{OUTPUT_DIR}/06_emu_output/{{sample}}/{{sample}}_rel-abundance.tsv",
            sample=SAMPLES,
        ),
    output:
        html = f"{OUTPUT_DIR}/performance_report.html",
        csv  = f"{OUTPUT_DIR}/performance_summary.csv",
        json = f"{OUTPUT_DIR}/performance.json",
    params:
        bench_dir    = f"{OUTPUT_DIR}/benchmarks",
        raw_dir      = f"{OUTPUT_DIR}/02_nanostat_raw",
        filtered_dir = f"{OUTPUT_DIR}/05_nanostat_filtered",
        # The real core budget the run was given, not a config value: it is
        # what makes "how much of the machine did we use" answerable.
        cores        = _total_cores,
        min_length   = config["min_length"],
        max_length   = config["max_length"],
        min_quality  = config["min_quality"],
        db           = config["emu_db"],
        version      = config.get("nano16s_version", ""),
    script:
        "../scripts/make_perf_report.py"
