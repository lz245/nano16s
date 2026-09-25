# =============================================================================
# emu.smk — Emu probabilistic classification
# =============================================================================
# Requires from common.smk: OUTPUT_DIR, SAMPLES, EMU_RANKS
# Requires from preprocess.smk: {OUTPUT_DIR}/04_filtered/{sample}_filtered.fastq.gz
# Requires from config: emu_db
# =============================================================================

import os

def classifier_input(wildcards):
    """The reads Emu classifies: subsampled when max_reads is set, else filtered.

    A function rather than a fixed path, so the extra stage joins the DAG only
    when it is switched on. With max_reads at 0 the graph is identical to the
    one that ran before.
    """
    if int(config.get("max_reads", 0) or 0) > 0:
        return f"{OUTPUT_DIR}/04b_subsampled/{wildcards.sample}_subsampled.fastq.gz"
    return f"{OUTPUT_DIR}/04_filtered/{wildcards.sample}_filtered.fastq.gz"


rule emu:
    input:
        classifier_input
    output:
        f"{OUTPUT_DIR}/06_emu_output/{{sample}}/{{sample}}_rel-abundance.tsv"
    params:
        db     = config["emu_db"],
        outdir = f"{OUTPUT_DIR}/06_emu_output/{{sample}}",
        # --keep-read-assignments makes Emu write its read-by-taxon
        # distribution, which rule per_read turns into one line per read. It
        # costs no measurable time (Emu does the same work either way) but the
        # file is reads x taxa and can reach gigabytes for one deep barcode, so
        # it is written only when asked for, and deleted once converted.
        extra  = "--keep-read-assignments" if PER_READ else "",
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
        mkdir -p "{params.outdir}"

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
            : > "{output}"
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
        rm -f {params.outdir}/*_rel-abundance*.tsv {params.outdir}/*_counts*.tsv \
               {params.outdir}/*_read-assignment-distributions.tsv

        emu abundance \
            "{input}" \
            --db "{params.db}" \
            --keep-counts \
            {params.extra} \
            --output-dir "{params.outdir}" \
            --threads {threads}

        # `*_rel-abundance.tsv`, matching the full table exactly, NOT
        # `*_rel-abundance*.tsv`.
        #
        # Emu always writes <input>_rel-abundance.tsv, and additionally writes
        # <input>_rel-abundance-threshold-<min-abundance>.tsv whenever any taxon
        # falls below --min-abundance (default 0.0001). The wildcard matched
        # both, and "-threshold-" sorts before ".tsv", so the *thresholded*
        # table was the one promoted to the declared output while the full
        # result was left behind. A barcode with any taxon under 0.01% shipped a
        # thresholded, re-normalised table -- 160 taxa instead of 165, with
        # every abundance shifted -- while its neighbours shipped full ones. The
        # barcodes within a single run were therefore not comparable, and
        # nothing said so.
        #
        # A glob into an array, not `ls ... | head -1`: Snakemake runs shell
        # bodies under `set -euo pipefail`, so `ls` finding nothing exits 2, the
        # substitution fails, and the rule dies at the assignment -- before the
        # message below can explain why. That message has never been reachable.
        shopt -s nullglob
        REL=( {params.outdir}/*_rel-abundance.tsv )
        if [ "${{#REL[@]}}" -eq 0 ]; then
            echo "ERROR: Emu produced no output for {wildcards.sample}" >&2
            echo "  Expected a *_rel-abundance.tsv under {params.outdir}" >&2
            exit 1
        fi
        if [ "${{REL[0]}}" != "{output}" ]; then
            mv -f "${{REL[0]}}" {output}
        fi

        # Emu's thresholded table is a legitimate secondary output, so keep it
        # -- but under a predictable name, so it reads as what it is and cannot
        # be mistaken for the result again.
        #
        # The count is checked before the array is expanded, as it is for REL
        # and CNT. Not style: bash before 4.4 treats "${{THRESH[@]}}" as an
        # unset variable when the array is empty, and this runs under `set -u`,
        # so the expansion aborts the rule with `THRESH[@]: unbound variable`.
        # macOS ships bash 3.2 as /bin/bash, which is what Snakemake runs, and
        # Emu writes a thresholded table only when some taxon falls below
        # --min-abundance -- so on macOS the rule failed for every barcode that
        # did *not* trigger the threshold, which is the ordinary case.
        THRESH=( {params.outdir}/*_rel-abundance-threshold-*.tsv )
        if [ "${{#THRESH[@]}}" -gt 0 ]; then
            for t in "${{THRESH[@]}}"; do
                keep="{params.outdir}/{wildcards.sample}_rel-abundance-threshold-${{t##*-threshold-}}"
                [ "$t" = "$keep" ] || mv -f "$t" "$keep"
            done
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
        label_script = os.path.join(workflow.basedir, "scripts", "label_unclassified.py"),
    shell:
        """
        mkdir -p "{params.combined_dir}"
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
            cd "{params.combined_dir}"
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

        test -s "{output.rel}"
        test -s "{output.cnts}"

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
        if head -1 "{output.cnts}" | grep -q '^#'; then
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

        # Emu leaves one row with every taxonomy field empty: the reads it
        # could not place. Unlabelled, it reads as a blank line, so summing a
        # column silently includes it and filtering out unnamed rows silently
        # drops it. Name it, and the table adds up in plain sight.
        python3 "{params.label_script}" {wildcards.rank} "{output.rel}" "{output.cnts}"
        """


# ---------------------------------------------------------------------------
# Where every read went, per barcode, and how many taxa it found.
#
# The combined tables say what is in each sample; this says whether it adds up.
# One row per barcode: raw, removed by the filter, given to the classifier,
# classified, unclassified, and the number of species and genera found. The
# `check` column is the arithmetic a reader would otherwise have to do.
#
# per_barcode_taxa.tsv answers the other half: which species, not just how
# many. One row per barcode and species, with reads and share. The combined
# tables hold the same numbers as a grid, which suits a heatmap; this suits a
# person filtering one barcode, or a spreadsheet.
# ---------------------------------------------------------------------------
rule read_accounting:
    input:
        rel = expand(
            f"{OUTPUT_DIR}/06_emu_output/{{sample}}/{{sample}}_rel-abundance.tsv",
            sample=SAMPLES,
        ),
        summary = f"{OUTPUT_DIR}/preprocessing_summary.csv",
    output:
        accounting = f"{OUTPUT_DIR}/07_emu_combined/read_accounting.tsv",
        taxa       = f"{OUTPUT_DIR}/07_emu_combined/per_barcode_taxa.tsv",
    params:
        emu_dir = f"{OUTPUT_DIR}/06_emu_output",
        script  = os.path.join(workflow.basedir, "scripts", "read_accounting.py"),
    shell:
        """
        python3 "{params.script}" "{params.emu_dir}" "{input.summary}" \
            "{output.accounting}" "{output.taxa}"
        """


# ---------------------------------------------------------------------------
# One line per sequencing read: what it was called, and how sure Emu is.
#
# Only built with --per-read. Emu is probabilistic -- it spreads each read over
# the references it matched rather than labelling it -- so this reports the
# taxon holding most of that read's probability, with the probability itself,
# so an ambiguous read is visible as ambiguous rather than silently rounded to
# a species name.
#
# The distribution Emu writes is reads x taxa and can reach gigabytes for one
# barcode. It is converted here and then removed: what remains is a compact
# gzipped table, about a hundred bytes per read.
# ---------------------------------------------------------------------------
rule per_read:
    input:
        abundance = f"{OUTPUT_DIR}/06_emu_output/{{sample}}/{{sample}}_rel-abundance.tsv",
        reads     = classifier_input,
    output:
        f"{OUTPUT_DIR}/08_per_read/{{sample}}_per_read.tsv.gz"
    params:
        db     = config["emu_db"],
        outdir = f"{OUTPUT_DIR}/06_emu_output/{{sample}}",
        script = os.path.join(workflow.basedir, "scripts", "per_read_calls.py"),
    shell:
        """
        # Emu names the file after the input, so find it rather than guess.
        shopt -s nullglob
        DIST=( {params.outdir}/*_read-assignment-distributions.tsv )
        if [ "${{#DIST[@]}}" -eq 0 ]; then
            echo "ERROR: no read assignments for {wildcards.sample}" >&2
            echo "  Emu writes them only when nano16s is run with --per-read." >&2
            echo "  If this run already classified without it, the reads must be" >&2
            echo "  classified again: delete {params.outdir} and re-run with --per-read." >&2
            exit 1
        fi

        python3 "{params.script}" "${{DIST[0]}}" "{params.db}" "{input.reads}" \
            "{wildcards.sample}" "{output}"

        # Converted, so the matrix has served its purpose. Removing it here
        # rather than declaring it temp() keeps the emu rule's output list --
        # and so every existing output directory -- unchanged.
        rm -f "${{DIST[0]}}"
        """
