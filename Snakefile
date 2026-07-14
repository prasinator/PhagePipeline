# =============================================================================
# Therapeutic Phage Candidate Selection Pipeline — Astraphage Innovations
#
# Raw tool outputs  : results/{sample}/{tool}/
# Human-facing      : reports/{run_id}/{sample}/
# Run scores        : reports/{run_id}/phage_scores.tsv
#
# Usage:
#   Single/few : snakemake --cores 8 --config samples="T4" run_id="2026-06-26_test"
#   Batch file : snakemake --cores 8 --config batch="batches/ecoli.txt" run_id="2026-06-26_ecoli"
# =============================================================================

import os
from datetime import date

configfile: "config.yaml"

# ── run_id ────────────────────────────────────────────────────────────────────
RUN_ID = config.get("run_id", date.today().strftime("%Y-%m-%d"))

# ── sample list — batch file takes priority over inline samples ───────────────
if "batch" in config:
    _batch_file = config["batch"]
    if not os.path.exists(_batch_file):
        raise ValueError(f"Batch file not found: {_batch_file}")
    with open(_batch_file) as _fh:
        SAMPLES = [l.strip() for l in _fh if l.strip() and not l.startswith("#")]
    if not SAMPLES:
        raise ValueError(f"Batch file is empty: {_batch_file}")
    print(f"Batch mode: {len(SAMPLES)} samples from {_batch_file}")
elif "samples" in config:
    SAMPLES = [s.strip() for s in config["samples"].split(",") if s.strip()]
else:
    raise ValueError(
        "No samples specified.\n"
        "  Single/few:  --config samples=\"PhageA,PhageB\"\n"
        "  Batch file:  --config batch=\"batches/my_batch.txt\""
    )

# ── paths ─────────────────────────────────────────────────────────────────────
ROOT      = os.path.dirname(os.path.abspath(workflow.snakefile))
RAW_INPUT = os.path.join(ROOT, config["raw_input"])
RESULTS   = os.path.join(ROOT, config["results"])
REPORTS   = os.path.join(ROOT, config["reports"], RUN_ID)
SCRIPTS   = os.path.join(ROOT, config["scripts"])
DATABASES = os.path.join(ROOT, config["databases"])
THREADS   = config["threads"]


# =============================================================================
# rule all
# =============================================================================

rule all:
    input:
        expand(os.path.join(REPORTS, "{sample}", "{sample}_annotated.gbk"),    sample=SAMPLES),
        expand(os.path.join(REPORTS, "{sample}", "{sample}_protein_table.tsv"),sample=SAMPLES),
        expand(os.path.join(REPORTS, "{sample}", "{sample}_genome_map.png"),   sample=SAMPLES),
        expand(os.path.join(REPORTS, "{sample}", "{sample}_tool_summary.tsv"), sample=SAMPLES),
        expand(os.path.join(REPORTS, "{sample}", "{sample}_report.txt"),       sample=SAMPLES),
        os.path.join(REPORTS, "phage_scores.tsv"),
        expand(os.path.join(RESULTS, "{sample}", "phagerbpdetect", "{sample}_rbp_predictions.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS, "{sample}", "phagerbpdetect", "{sample}_rbp_summary.tsv"),      sample=SAMPLES),
        expand(os.path.join(RESULTS, "{sample}", "phagetailfinder", "{sample}_tailfinder.tsv"),      sample=SAMPLES),


# =============================================================================
# Stage 0 — Quality
# =============================================================================

rule checkv:
    input:  genome=os.path.join(RAW_INPUT, "{sample}.fasta")
    output: os.path.join(RESULTS, "{sample}", "checkv", "quality_summary.tsv")
    params:
        script=os.path.join(SCRIPTS, "checkv", "run_checkv.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "checkv"),
        threads=THREADS
    log: os.path.join(RESULTS, "{sample}", "logs", "checkv.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.genome} {wildcards.sample} {params.threads} {params.outdir} &> {log}
        """


# =============================================================================
# Stage 1 — Annotation chain
# =============================================================================

rule pharokka:
    input:
        genome=os.path.join(RAW_INPUT, "{sample}.fasta")
    output:
        gbk=os.path.join(
            RESULTS,
            "{sample}",
            "pharokka",
            "{sample}.gbk"
        ),
        faa=os.path.join(
            RESULTS,
            "{sample}",
            "pharokka",
            "{sample}.faa"
        ),
        ffn=os.path.join(
            RESULTS,
            "{sample}",
            "pharokka",
            "{sample}_cds.ffn"
        ),
        tsv=os.path.join(
            RESULTS,
            "{sample}",
            "pharokka",
            "{sample}_length_gc_cds_density.tsv"
        )
    params:
        script=os.path.join(SCRIPTS, "pharokka", "run_pharokka.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "pharokka"),
        threads=THREADS
    log: os.path.join(RESULTS, "{sample}", "logs", "pharokka.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.genome} {wildcards.sample} {params.threads} {params.outdir} &> {log}
        """


rule phold:
    input:  gbk=os.path.join(RESULTS, "{sample}", "pharokka", "{sample}.gbk")
    output:
        gbk    =os.path.join(RESULTS, "{sample}", "phold", "{sample}.gbk"),
        cds_tsv=os.path.join(RESULTS, "{sample}", "phold", "phold_per_cds_predictions.tsv")
    params:
        script=os.path.join(SCRIPTS, "phage", "run_phold.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "phold")
    log: os.path.join(RESULTS, "{sample}", "logs", "phold.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.gbk} {wildcards.sample} {params.outdir} &> {log}
        """


rule phynteny:
    input:  gbk=os.path.join(RESULTS, "{sample}", "phold", "{sample}.gbk")
    output: directory(os.path.join(RESULTS, "{sample}", "phynteny"))
    params:
        script=os.path.join(SCRIPTS, "phage", "run_phynteny.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "phynteny")
    log: os.path.join(RESULTS, "{sample}", "logs", "phynteny.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.gbk} {wildcards.sample} {params.outdir} &> {log}
        """
 
# =============================================================================
# NEW: PhageDPO
# Input : pharokka .ffn
# Output: HTML + TSV
#
# Depends on pharokka .ffn — runs independently of phold/phynteny
# =============================================================================
 
rule phagedpo:
    input:
        ffn=os.path.join(RESULTS, "{sample}", "pharokka", "{sample}_cds.ffn")
    output:
        html=os.path.join(RESULTS, "{sample}", "phagedpo", "{sample}_cds_output.html"),
        tsv=os.path.join(RESULTS, "{sample}", "phagedpo",  "{sample}_cds_output.tsv")
    params:
        script  =os.path.join(SCRIPTS, "phagedpo", "run_phagedpo.sh"),
        pharokka=lambda wc: os.path.join(RESULTS, wc.sample, "pharokka"),
        outdir  =lambda wc: os.path.join(RESULTS, wc.sample, "phagedpo")
    log: os.path.join(RESULTS, "{sample}", "logs", "phagedpo.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} \
            {params.pharokka} \
            {wildcards.sample} \
            {params.outdir} \
            &> {log}
        """
 
 
# =============================================================================
# NEW: PhageRBPdetect v4
# Input : pharokka .faa
# Output: per-protein predictions TSV + summary TSV
#
# NOTE: inference serialised via flock — safe with --cores N
# =============================================================================
 
rule phagerbpdetect:
    input:
        faa=os.path.join(RESULTS, "{sample}", "pharokka", "{sample}.faa")
    output:
        predictions=os.path.join(
            RESULTS, "{sample}", "phagerbpdetect", "{sample}_rbp_predictions.tsv"
        ),
        summary=os.path.join(
            RESULTS, "{sample}", "phagerbpdetect", "{sample}_rbp_summary.tsv"
        )
    params:
        script=os.path.join(SCRIPTS, "phagerbpdetect", "run_phagerbpdetect_v4.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "phagerbpdetect")
    log: os.path.join(RESULTS, "{sample}", "logs", "phagerbpdetect.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} \
            {input.faa} \
            {wildcards.sample} \
            {params.outdir} \
            &> {log}
        """
 
 
# =============================================================================
# NEW: PhageTailFinder
# Input : pharokka .faa
# Output: per-protein tail classification TSV
#
# Depends on pharokka .faa (not gbk) — runs independently of phold/phynteny
# =============================================================================
 
rule phagetailfinder:
    input:
        faa=os.path.join(RESULTS, "{sample}", "pharokka", "{sample}.faa")
    output:
        tsv=os.path.join(
            RESULTS, "{sample}", "phagetailfinder", "{sample}_tailfinder.tsv"
        )
    params:
        script=os.path.join(SCRIPTS, "phagetailfinder", "run_phagetailfinder.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "phagetailfinder")
    log: os.path.join(RESULTS, "{sample}", "logs", "phagetailfinder.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} \
            {input.faa} \
            {wildcards.sample} \
            {params.outdir} \
            &> {log}
        """

# =============================================================================
# Stage 2 — Lifestyle
# =============================================================================

rule bacphlip:
    input:  genome=os.path.join(RAW_INPUT, "{sample}.fasta")
    output: os.path.join(RESULTS, "{sample}", "bacphlip", "{sample}.fasta.bacphlip")
    params:
        script=os.path.join(SCRIPTS, "bacphlip", "run_bacphlip.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "bacphlip")
    log: os.path.join(RESULTS, "{sample}", "logs", "bacphlip.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.genome} {wildcards.sample} {params.outdir} &> {log}
        """


# =============================================================================
# Stage 2 — Safety
# =============================================================================

rule abricate_card:
    input:  genome=os.path.join(RAW_INPUT, "{sample}.fasta")
    output: os.path.join(RESULTS, "{sample}", "abricate", "{sample}_card.tsv")
    params:
        script=os.path.join(SCRIPTS, "abricate", "run_abricate.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "abricate"),
        db="card"
    log: os.path.join(RESULTS, "{sample}", "logs", "abricate_card.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.genome} {wildcards.sample} {params.db} {params.outdir} &> {log}
        """


rule abricate_vfdb:
    input:  genome=os.path.join(RAW_INPUT, "{sample}.fasta")
    output: os.path.join(RESULTS, "{sample}", "abricate", "{sample}_vfdb.tsv")
    params:
        script=os.path.join(SCRIPTS, "abricate", "run_abricate.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "abricate"),
        db="vfdb"
    log: os.path.join(RESULTS, "{sample}", "logs", "abricate_vfdb.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.genome} {wildcards.sample} {params.db} {params.outdir} &> {log}
        """


rule rgi:
    input:  genome=os.path.join(RAW_INPUT, "{sample}.fasta")
    output: os.path.join(RESULTS, "{sample}", "rgi", "{sample}_rgi.txt")
    params:
        script=os.path.join(SCRIPTS, "rgi", "run_rgi.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "rgi"),
        threads=THREADS
    log: os.path.join(RESULTS, "{sample}", "logs", "rgi.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.genome} {wildcards.sample} contig {params.threads} {params.outdir} &> {log}
        """


rule acrdb:
    input:  faa=os.path.join(RESULTS, "{sample}", "pharokka", "{sample}.faa")
    output: os.path.join(RESULTS, "{sample}", "acrdb", "{sample}_acrdb_blastp.tsv")
    params:
        script=os.path.join(SCRIPTS, "acrdb", "run_acrdb_blast.sh"),
        outdir=lambda wc: os.path.join(RESULTS, wc.sample, "acrdb"),
        db=os.path.join(ROOT, config["acrdb_db"])
    log: os.path.join(RESULTS, "{sample}", "logs", "acrdb.log")
    shell:
        """
        mkdir -p $(dirname {log})
        {params.script} {input.faa} {wildcards.sample} {params.db} {params.outdir} &> {log}
        """


# =============================================================================
# Stage 3 — Post-processing (reports/{run_id}/{sample}/)
# =============================================================================

rule build_annotated_gb:
    input:
        gbk     =os.path.join(RESULTS, "{sample}", "phold", "{sample}.gbk"),
        cds_tsv =os.path.join(RESULTS, "{sample}", "phold", "phold_per_cds_predictions.tsv"),
        phynteny=os.path.join(RESULTS, "{sample}", "phynteny")
    output: os.path.join(REPORTS, "{sample}", "{sample}_annotated.gbk")
    params:
        script      =os.path.join(SCRIPTS, "build_annotated_gb.py"),
        phynteny_tsv=lambda wc: os.path.join(RESULTS, wc.sample, "phynteny", "phynteny.tsv"),
        outdir      =lambda wc: os.path.join(REPORTS, wc.sample)
    log: os.path.join(RESULTS, "{sample}", "logs", "build_annotated_gb.log")
    shell:
        """
        mkdir -p {params.outdir}
        conda run -n phage python {params.script} \
            {input.gbk} {input.cds_tsv} {params.phynteny_tsv} \
            {wildcards.sample} {params.outdir} &> {log}
        """


rule build_protein_table:
    input:  os.path.join(REPORTS, "{sample}", "{sample}_annotated.gbk")
    output: os.path.join(REPORTS, "{sample}", "{sample}_protein_table.tsv")
    params:
        script=os.path.join(SCRIPTS, "build_protein_table.py"),
        outdir=lambda wc: os.path.join(REPORTS, wc.sample)
    log: os.path.join(RESULTS, "{sample}", "logs", "build_protein_table.log")
    shell:
        """
        conda run -n phage python {params.script} \
            {input} {wildcards.sample} {params.outdir} &> {log}
        """


rule visualize_genome:
    input:
        gbk=os.path.join(REPORTS, "{sample}", "{sample}_annotated.gbk"),
        tsv=os.path.join(REPORTS, "{sample}", "{sample}_protein_table.tsv")
    output:
        png=os.path.join(REPORTS, "{sample}", "{sample}_genome_map.png"),
        svg=os.path.join(REPORTS, "{sample}", "{sample}_genome_map.svg")
    params:
        script    =os.path.join(SCRIPTS, "visualize_genome.py"),
        out_prefix=lambda wc: os.path.join(REPORTS, wc.sample, f"{wc.sample}_genome_map")
    log: os.path.join(RESULTS, "{sample}", "logs", "visualize_genome.log")
    shell:
        """
        conda run -n phage python {params.script} \
            {input.gbk} {input.tsv} {params.out_prefix} &> {log}
        """


rule summarize_sample:
    input:
        checkv  =os.path.join(RESULTS, "{sample}", "checkv",   "quality_summary.tsv"),
        bacphlip=os.path.join(RESULTS, "{sample}", "bacphlip", "{sample}.fasta.bacphlip"),
        rgi     =os.path.join(RESULTS, "{sample}", "rgi",      "{sample}_rgi.txt"),
        card    =os.path.join(RESULTS, "{sample}", "abricate", "{sample}_card.tsv"),
        vfdb    =os.path.join(RESULTS, "{sample}", "abricate", "{sample}_vfdb.tsv"),
        acrdb   =os.path.join(RESULTS, "{sample}", "acrdb",    "{sample}_acrdb_blastp.tsv"),
        protein =os.path.join(REPORTS, "{sample}", "{sample}_protein_table.tsv"),
        pharokka=os.path.join(RESULTS, "{sample}", "pharokka", "{sample}_length_gc_cds_density.tsv")
    output:
        tsv=os.path.join(REPORTS, "{sample}", "{sample}_tool_summary.tsv"),
        txt=os.path.join(REPORTS, "{sample}", "{sample}_report.txt")
    params:
        script=os.path.join(SCRIPTS, "summarize_sample.py"),
        sample="{sample}"
    log: os.path.join(RESULTS, "{sample}", "logs", "summarize_sample.log")
    shell:
        """
        conda run -n phage python {params.script} \
            {params.sample} \
            {input.checkv} \
            {input.bacphlip} \
            {input.rgi} \
            {input.card} \
            {input.vfdb} \
            {input.acrdb} \
            {input.protein} \
            {input.pharokka} \
            {output.tsv} \
            {output.txt} &> {log}
        """


# =============================================================================
# Stage 4 — Scoring
# =============================================================================

rule score_and_rank:
    input:
        expand(os.path.join(REPORTS, "{sample}", "{sample}_tool_summary.tsv"), sample=SAMPLES)
    output:
        os.path.join(REPORTS, "phage_scores.tsv")
    params:
        script     =os.path.join(SCRIPTS, "score_and_rank.py"),
        reports_dir=REPORTS,
        samples    =",".join(SAMPLES)
    log: os.path.join(RESULTS, "logs", f"{RUN_ID}_score_and_rank.log")
    shell:
        """
        mkdir -p $(dirname {log})
        conda run -n phage python {params.script} \
            {params.reports_dir} \
            "{params.samples}" \
            {output} &> {log}
        """
