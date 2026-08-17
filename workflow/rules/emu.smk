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
    benchmark:
        f"{OUTPUT_DIR}/benchmarks/emu/{{sample}}.tsv"
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

        # Clear results from any earlier attempt before Emu runs, so whatever
        # is here afterwards is unambiguously this job's own output.
        #
        # Emu names its file after the input, giving
        # barcodeNN_filtered.fastq_rel-abundance.tsv, which the rename below
        # turns into the declared name. An attempt interrupted between those
        # two steps left Emu's name behind -- and the rule does not re-run once
        # its declared output exists, so the stray survived every later run.
        # emu_combine then read it as an extra sample: 32 columns in a
        # 24-barcode table, each duplicate a complete second abundance profile
        # that disagreed with the real one, with nothing to say so.
        rm -f {params.outdir}/*_rel-abundance*.tsv {params.outdir}/*_counts*.tsv

        emu abundance \
            {input} \
            --db {params.db} \
            --keep-counts \
            --output-dir {params.outdir} \
            --threads {threads}

        # A glob into an array, not `ls ... | head -1`. Snakemake runs shell
        # bodies under `set -euo pipefail`, so `ls` finding nothing exits 2, the
        # substitution fails, and the rule dies at the assignment -- before the
        # message below can explain why. That message has never been reachable.
        shopt -s nullglob
        REL=( {params.outdir}/*_rel-abundance*.tsv )
        if [ "${{#REL[@]}}" -eq 0 ]; then
            echo "ERROR: Emu produced no output for {wildcards.sample}" >&2
            echo "  Expected a *_rel-abundance.tsv under {params.outdir}" >&2
            exit 1
        fi
        if [ "${{REL[0]}}" != "{output}" ]; then
            mv -f "${{REL[0]}}" {output}
        fi

        CNT=( {params.outdir}/*_counts*.tsv )
        if [ "${{#CNT[@]}}" -gt 0 ] \
           && [ "${{CNT[0]}}" != "{params.outdir}/{wildcards.sample}_counts.tsv" ]; then
            mv -f "${{CNT[0]}}" {params.outdir}/{wildcards.sample}_counts.tsv
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

        # Emu combine-outputs expects per-sample TSVs in one directory, so link
        # them in from the per-barcode subdirectories.
        #
        # The rule's declared inputs, not a wildcard glob. `*/*_rel-abundance*`
        # matched anything that happened to be sitting in a barcode directory,
        # and every match became a column in the combined table. These paths are
        # exactly one per barcode, by construction, so a stray file can no
        # longer become a sample.
        REL_COUNT=0
        for TSV in {input}; do
            [ -s "$TSV" ] || continue
            ln -sf "$TSV" "$COMBINE_INPUT/$(basename "$TSV")"
            REL_COUNT=$(( REL_COUNT + 1 ))
            CNT="$(dirname "$TSV")/$(basename "$TSV" _rel-abundance.tsv)_counts.tsv"
            if [ -s "$CNT" ]; then
                ln -sf "$CNT" "$COMBINE_INPUT/$(basename "$CNT")"
            fi
        done

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
