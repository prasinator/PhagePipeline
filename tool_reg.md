# Astraphage Tool Registry

## Pharokka

Status:
Verified

Purpose:
Primary phage genome annotation.

Input:
FASTA

Important Outputs:

{sample}.gbk

{sample}.gff

phanotate.faa

phanotate.ffn

{sample}_cds_final_merged_output.tsv

{sample}_cds_functions.tsv

Downstream:
PHOLD
PhageTailFinder
PhageRBPdetect
PhageDPO
RGI
Reporting

---

## PHOLD

Status:
Verified

Purpose:
Structural homology-based annotation refinement.

Input:
Pharokka protein FASTA (FAA)

Important Outputs:

phold.gbk

phold_per_cds_predictions.tsv

phold_all_cds_functions.tsv

Downstream:
Phynteny
Reporting

Observed Validation:

T7 comparison:

61 CDS

55 exact matches

2 refined annotations

3 hypothetical proteins annotated

1 broader annotation

0 major contradictions

---

## Phynteny

Status:
Testing

Purpose:
Synteny-based functional prediction.

Planned Input:
phold.gbk

Expected Outputs:
Functional predictions for hypothetical proteins.

Downstream:
Reporting

---

## PhageTailFinder

Status:
Testing

Purpose:
Tail protein prediction.

Input:
Pharokka FAA

Outputs:
Tail protein predictions.

Notes:
Observed false negative on T4.

Use as supplementary evidence.

---

## PhageRBPdetect

Status:
Testing

Purpose:
Receptor-binding protein prediction.

Input:
Pharokka FAA

Outputs:
RBP predictions.

---

## PhageDPO

Status:
Testing

Purpose:
Depolymerase prediction.

Input:
Pharokka FAA

Outputs:
DPO predictions.

Observed:
Accurate true-negative results on T4/T7 validation.

---

## CheckV

Status:
Verified

Purpose:
Genome quality assessment.

Input:
FASTA

Outputs:
Completeness and contamination metrics.

---

## Bacphlip

Status:
Verified

Purpose:
Lifestyle prediction.

Input:
FASTA

Outputs:
Virulent / Temperate prediction scores.

---

## RGI

Status:
Verified

Purpose:
AMR gene detection.

Input:
Pharokka FAA

Outputs:
AMR predictions.

---

## Abricate

Status:
Verified

Purpose:
Resistance and virulence screening.

Input:
FASTA

Outputs:
Database screening results.

---

## AcrDB

Status:
Verified

Purpose:
Anti-CRISPR detection.

Input:
FASTA or predicted proteins depending on workflow implementation.

Outputs:
Anti-CRISPR hits.

Notes:
Weak hits should be filtered using significance thresholds.

---

## Future Integrations

Planned:

* DefenseFinder
* InterProScan

Do not assume these tools are integrated until explicitly implemented.
