# =============================================================================
# preprocess.smk — merge, QC, trim, filter
# =============================================================================
# Shared by Snakefile.emu and Snakefile.minimap2. Both pipelines write to the
# same 01_-05_ paths, so running the second one reuses this work rather than
# repeating it.
#
# Requires from common.smk: INPUT_DIR, OUTPUT_DIR, SAMPLES
# Produces for classifiers: {OUTPUT_DIR}/04_filtered/{sample}_filtered.fastq.gz
# =============================================================================

import os

# ---------------------------------------------------------------------------
# Rule 01: Merge per-barcode FASTQ files
# ---------------------------------------------------------------------------
rule merge:
    input:
        barcode_dir = os.path.join(INPUT_DIR, "{sample}")
    output:
        f"{OUTPUT_DIR}/01_merged/{{sample}}.fastq.gz"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/merge/{{sample}}.tsv"
    resources:
        cpus_per_task = config["resources"]["merge"]["cpus"],
        mem_mb        = config["resources"]["merge"]["mem_mb"],
        runtime       = config["resources"]["merge"]["time_min"],
    shell:
        """
        # nullglob so an empty directory yields an empty array rather than the
        # literal pattern. Without it the array holds one element -- the
        # unexpanded glob -- and the single-file branch below fails with
        # "cp: cannot stat '.../*.fastq.gz'", which reads like a corrupted path
        # rather than an empty directory.
        shopt -s nullglob
        FASTQ_FILES=( "{input.barcode_dir}"/*.fastq.gz )

        if [ "${{#FASTQ_FILES[@]}}" -eq 0 ]; then
            echo "ERROR: no .fastq.gz files in {input.barcode_dir}" >&2
            echo "  Every barcode directory must contain at least one .fastq.gz file." >&2
            echo "  If this barcode produced no reads, remove the directory and re-run." >&2
            exit 1
        fi

        if [ "${{#FASTQ_FILES[@]}}" -eq 1 ]; then
            cp "${{FASTQ_FILES[0]}}" {output}
        else
            cat "${{FASTQ_FILES[@]}}" > {output}
        fi
        """


# ---------------------------------------------------------------------------
# Rule 02: NanoStat on raw merged reads
# ---------------------------------------------------------------------------
rule nanostat_raw:
    input:
        f"{OUTPUT_DIR}/01_merged/{{sample}}.fastq.gz"
    output:
        f"{OUTPUT_DIR}/02_nanostat_raw/{{sample}}/{{sample}}_quality_summary.txt"
    params:
        outdir = f"{OUTPUT_DIR}/02_nanostat_raw/{{sample}}"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/nanostat_raw/{{sample}}.tsv"
    resources:
        cpus_per_task = config["resources"]["nanostat"]["cpus"],
        mem_mb        = config["resources"]["nanostat"]["mem_mb"],
        runtime       = config["resources"]["nanostat"]["time_min"],
    shell:
        """
        mkdir -p {params.outdir}
        if python3 - "{input}" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt", errors="replace") as handle:
    for i, line in enumerate(handle):
        if i == 1 and line.strip():
            sys.exit(0)
sys.exit(1)
PY
        then
            NanoStat --fastq {input} \
                --name {wildcards.sample}_quality_summary.txt \
                --outdir {params.outdir}
        else
            cat > {output} <<'EOF'
Number of reads: 0
Total bases: 0
Median read length: 0
Median read quality: 0
EOF
        fi
        test -s {output}
        """


# ---------------------------------------------------------------------------
# Rule 03: Porechop_ABI adapter trimming
# ---------------------------------------------------------------------------
rule porechop:
    input:
        f"{OUTPUT_DIR}/01_merged/{{sample}}.fastq.gz"
    output:
        f"{OUTPUT_DIR}/03_trimmed/{{sample}}_trimmed.fastq.gz"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/porechop/{{sample}}.tsv"
    threads:
        config["resources"]["porechop"]["cpus"]
    resources:
        cpus_per_task = config["resources"]["porechop"]["cpus"],
        mem_mb        = config["resources"]["porechop"]["mem_mb"],
        runtime       = config["resources"]["porechop"]["time_min"],
    shell:
        """
        TMPDIR="${{TMPDIR:-/tmp}}"
        mkdir -p "$TMPDIR/porechop_tmp_{wildcards.sample}"
        cd "$TMPDIR/porechop_tmp_{wildcards.sample}"
        porechop_abi -abi \
            -i {input} \
            -o {output} \
            --threads {threads}
        rm -rf "$TMPDIR/porechop_tmp_{wildcards.sample}"
        """


# ---------------------------------------------------------------------------
# Rule 04: Chopper quality + length filtering
# ---------------------------------------------------------------------------
rule chopper:
    input:
        f"{OUTPUT_DIR}/03_trimmed/{{sample}}_trimmed.fastq.gz"
    output:
        f"{OUTPUT_DIR}/04_filtered/{{sample}}_filtered.fastq.gz"
    params:
        min_len = config["min_length"],
        max_len = config["max_length"],
        min_q   = config["min_quality"],
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/chopper/{{sample}}.tsv"
    threads:
        config["resources"]["chopper"]["cpus"]
    resources:
        cpus_per_task = config["resources"]["chopper"]["cpus"],
        mem_mb        = config["resources"]["chopper"]["mem_mb"],
        runtime       = config["resources"]["chopper"]["time_min"],
    shell:
        """
        gunzip -c {input} \
            | chopper \
                -l {params.min_len} \
                --maxlength {params.max_len} \
                -q {params.min_q} \
                --threads {threads} \
            | gzip > {output}
        """


# ---------------------------------------------------------------------------
# Rule 05: NanoStat on filtered reads
# ---------------------------------------------------------------------------
rule nanostat_filtered:
    input:
        f"{OUTPUT_DIR}/04_filtered/{{sample}}_filtered.fastq.gz"
    output:
        f"{OUTPUT_DIR}/05_nanostat_filtered/{{sample}}/{{sample}}_filtered_quality_summary.txt"
    params:
        outdir = f"{OUTPUT_DIR}/05_nanostat_filtered/{{sample}}"
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/nanostat_filtered/{{sample}}.tsv"
    resources:
        cpus_per_task = config["resources"]["nanostat"]["cpus"],
        mem_mb        = config["resources"]["nanostat"]["mem_mb"],
        runtime       = config["resources"]["nanostat"]["time_min"],
    shell:
        """
        mkdir -p {params.outdir}
        if python3 - "{input}" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt", errors="replace") as handle:
    for i, line in enumerate(handle):
        if i == 1 and line.strip():
            sys.exit(0)
sys.exit(1)
PY
        then
            NanoStat --fastq {input} \
                --name {wildcards.sample}_filtered_quality_summary.txt \
                --outdir {params.outdir}
        else
            cat > {output} <<'EOF'
Number of reads: 0
Total bases: 0
Median read length: 0
Median read quality: 0
EOF
        fi
        test -s {output}
        """


# ---------------------------------------------------------------------------
# Rule 05b: Generate preprocessing summary CSV
# ---------------------------------------------------------------------------
rule preprocessing_summary:
    input:
        raw = expand(
            f"{OUTPUT_DIR}/02_nanostat_raw/{{sample}}/{{sample}}_quality_summary.txt",
            sample=SAMPLES
        ),
        filtered = expand(
            f"{OUTPUT_DIR}/05_nanostat_filtered/{{sample}}/{{sample}}_filtered_quality_summary.txt",
            sample=SAMPLES
        ),
    output:
        f"{OUTPUT_DIR}/preprocessing_summary.csv"
    params:
        raw_dir      = f"{OUTPUT_DIR}/02_nanostat_raw",
        filtered_dir = f"{OUTPUT_DIR}/05_nanostat_filtered",
        samples      = SAMPLES,
    run:
        import csv

        def parse_nanostat(filepath):
            """Extract metrics from a NanoStat summary file."""
            metrics = {}
            with open(filepath) as f:
                for line in f:
                    parts = line.strip().split(":")
                    if len(parts) == 2:
                        key = parts[0].strip()
                        val = parts[1].strip().replace(",", "")
                        metrics[key] = val
            return metrics

        header = [
            "barcode",
            "raw_reads", "raw_bases", "raw_median_length", "raw_median_quality",
            "filtered_reads", "filtered_bases", "filtered_median_length", "filtered_median_quality",
        ]

        rows = []
        for sample in params.samples:
            raw_file = os.path.join(
                params.raw_dir, sample, f"{sample}_quality_summary.txt"
            )
            filt_file = os.path.join(
                params.filtered_dir, sample, f"{sample}_filtered_quality_summary.txt"
            )

            raw = parse_nanostat(raw_file)
            filt = parse_nanostat(filt_file)

            rows.append([
                sample,
                raw.get("Number of reads", ""),
                raw.get("Total bases", ""),
                raw.get("Median read length", ""),
                raw.get("Median read quality", ""),
                filt.get("Number of reads", ""),
                filt.get("Total bases", ""),
                filt.get("Median read length", ""),
                filt.get("Median read quality", ""),
            ])

        with open(output[0], "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)
