# =============================================================================
# emu.smk — Emu probabilistic classification
# =============================================================================
# Requires from common.smk: OUTPUT_DIR, SAMPLES, EMU_RANKS
# Requires from preprocess.smk: {OUTPUT_DIR}/04_filtered/{sample}_filtered.fastq.gz
# Requires from config: emu_db
# =============================================================================

rule emu:
    input:
        f"{OUTPUT_DIR}/04_filtered/{{sample}}_filtered.fastq.gz"
    output:
        f"{OUTPUT_DIR}/06_emu_output/{{sample}}/{{sample}}_rel-abundance.tsv"
    params:
        db     = config["emu_db"],
        outdir = f"{OUTPUT_DIR}/06_emu_output/{{sample}}",
    threads:
        config["resources"]["emu"]["cpus"]
    resources:
        cpus_per_task = config["resources"]["emu"]["cpus"],
        mem_mb        = config["resources"]["emu"]["mem_mb"],
        runtime       = config["resources"]["emu"]["time_min"],
    shell:
        """
        mkdir -p {params.outdir}

        # Pre-flight: verify database
        if [ ! -f "{params.db}/taxonomy.tsv" ] || [ ! -f "{params.db}/species_taxid.fasta" ]; then
            echo "ERROR: Emu database not found at {params.db}" >&2
            exit 1
        fi

        if ! python3 - "{input}" <<'PY'
import gzip
import sys

with gzip.open(sys.argv[1], "rt", errors="replace") as handle:
    for i, line in enumerate(handle):
        if i == 1 and line.strip():
            sys.exit(0)
sys.exit(1)
PY
        then
            echo "No filtered reads for {wildcards.sample}; writing empty Emu placeholder."
            : > {output}
            echo "no_filtered_reads" > {params.outdir}/{wildcards.sample}_emu_status.txt
            exit 0
        fi

        emu abundance \
            {input} \
            --db {params.db} \
            --keep-counts \
            --output-dir {params.outdir} \
            --threads {threads}

        # Emu names output with a hash — rename to predictable name
        RESULT=$(ls {params.outdir}/*_rel-abundance*.tsv 2>/dev/null | head -1)
        if [ -z "$RESULT" ]; then
            echo "ERROR: Emu produced no output for {wildcards.sample}" >&2
            exit 1
        fi
        if [ "$RESULT" != "{output}" ]; then
            mv "$RESULT" {output}
        fi

        COUNT_RESULT=$(find {params.outdir} -maxdepth 1 -type f -name '*_counts*.tsv' -size +0c | head -n 1)
        if [ -n "$COUNT_RESULT" ] && [ "$COUNT_RESULT" != "{params.outdir}/{wildcards.sample}_counts.tsv" ]; then
            mv "$COUNT_RESULT" {params.outdir}/{wildcards.sample}_counts.tsv
        fi
        """


# -----------------------------------------------------------------------
# Rule 06b: Combine Emu outputs per rank
# -----------------------------------------------------------------------
rule emu_combine:
    input:
        expand(
            f"{OUTPUT_DIR}/06_emu_output/{{sample}}/{{sample}}_rel-abundance.tsv",
            sample=SAMPLES,
        )
    output:
        rel  = f"{OUTPUT_DIR}/07_emu_combined/emu-combined-{{rank}}.tsv",
        cnts = f"{OUTPUT_DIR}/07_emu_combined/emu-combined-{{rank}}-counts.tsv",
    params:
        emu_dir     = f"{OUTPUT_DIR}/06_emu_output",
        combined_dir = f"{OUTPUT_DIR}/07_emu_combined",
        db          = config["emu_db"],
    shell:
        """
        mkdir -p {params.combined_dir}
        COMBINE_INPUT="{params.combined_dir}/.emu_combine_input_{wildcards.rank}"
        rm -rf "$COMBINE_INPUT"
        mkdir -p "$COMBINE_INPUT"

        # Emu combine-outputs expects per-sample TSVs in the emu_dir root.
        # Create symlinks from per-barcode subdirs for both abundance and counts.
        for TSV in {params.emu_dir}/*/*_rel-abundance*.tsv {params.emu_dir}/*/*_counts*.tsv; do
            [ -s "$TSV" ] || continue
            BASENAME=$(basename "$TSV")
            ln -sf "$TSV" "$COMBINE_INPUT/$BASENAME"
        done

        REL_COUNT=$(find "$COMBINE_INPUT" -maxdepth 1 -type l -name '*_rel-abundance*.tsv' | wc -l | tr -d ' ')

        if [ "$REL_COUNT" -eq 0 ]; then
            echo "# no non-empty Emu rel-abundance files available for {wildcards.rank}" > {output.rel}
            echo "# no non-empty Emu count files available for {wildcards.rank}" > {output.cnts}
            exit 0
        fi

        # Read counts are NOT a separate file. With --keep-counts, Emu writes an
        # "estimated counts" column inside each *_rel-abundance.tsv, and
        # `combine-outputs --counts` reads that column. Gating this on the
        # existence of a *_counts*.tsv file — as an earlier version did — meant
        # the counts tables were always empty placeholders.
        #
        # Run from combined_dir so Emu's output location is deterministic even
        # across Emu versions that write to the current working directory.
        (
            cd {params.combined_dir}
            rm -f emu-combined-{wildcards.rank}.tsv emu-combined-{wildcards.rank}-counts.tsv
            emu combine-outputs "$COMBINE_INPUT" {wildcards.rank}
            emu combine-outputs "$COMBINE_INPUT" {wildcards.rank} --counts
        )

        if [ ! -f {output.rel} ] && [ -f "$COMBINE_INPUT/emu-combined-{wildcards.rank}.tsv" ]; then
            mv "$COMBINE_INPUT/emu-combined-{wildcards.rank}.tsv" {output.rel}
        fi

        if [ ! -f {output.cnts} ] && [ -f "$COMBINE_INPUT/emu-combined-{wildcards.rank}-counts.tsv" ]; then
            mv "$COMBINE_INPUT/emu-combined-{wildcards.rank}-counts.tsv" {output.cnts}
        fi

        test -s {output.rel}
        test -s {output.cnts}

        # A counts table that is only the placeholder comment means the
        # estimated-counts column was missing upstream. Fail loudly rather than
        # shipping an empty table the README promises is populated.
        #
        # Two quite different causes land here, so name both. The second is
        # easy to misread as the first: on WSL2 the guest clock drifts from the
        # host and can resynchronise mid-run, leaving an output file with a
        # timestamp behind its own input. Snakemake treats that as a corrupted
        # build and deletes the output -- after this rule has already written
        # it -- so the run fails reporting missing counts when the real cause
        # was the clock.
        if head -1 {output.cnts} | grep -q '^#'; then
            echo "ERROR: the counts table for {wildcards.rank} is empty." >&2
            echo "" >&2
            echo "  Two things cause this:" >&2
            echo "" >&2
            echo "  1. Emu ran without --keep-counts, so no estimated-counts" >&2
            echo "     column was written." >&2
            echo "" >&2
            echo "  2. On WSL2, the clock drifted during the run and Snakemake" >&2
            echo "     removed the output as suspected clock skew. Look further" >&2
            echo "     up the log for 'has older modification time'. If it is" >&2
            echo "     there, run 'wsl --shutdown' from PowerShell and re-run;" >&2
            echo "     completed work is kept, so it finishes quickly." >&2
            exit 1
        fi
        """
