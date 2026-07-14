# PhagePipeline

**Therapeutic Phage Candidate Selection Pipeline**  
*Astraphage Innovations*

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-blue.svg)](https://snakemake.readthedocs.io)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org)
[![License](https://img.shields.io/badge/license-proprietary-red.svg)](#)

---

## Overview

PhagePipeline is an end-to-end Snakemake workflow for the systematic evaluation of bacteriophage genome candidates for therapeutic applications. Given one or more raw phage genome FASTA files, the pipeline executes a staged sequence of bioinformatics tools covering genome quality assessment, functional annotation, lifestyle prediction, safety screening, and host-interaction characterisation. All evidence is synthesised into a quantitative scoring system that ranks candidates as **PASS**, **REVIEW**, **FAIL**, or **REJECTED**.

The pipeline is designed for internal use at Astraphage Innovations and is structured to be shareable with collaborators under appropriate terms.

---

## Table of Contents

- [Overview](#overview)
- [Workflow Summary](#workflow-summary)
- [Repository Structure](#repository-structure)
- [Quick Start](#quick-start)
- [Scoring System](#scoring-system)
- [Output Files](#output-files)
- [Tool Summary](#tool-summary)
- [Documentation](#documentation)
- [Known Limitations](#known-limitations)

---

## Workflow Summary

The pipeline executes across five stages:

```mermaid
flowchart TD
    A[raw_input/{sample}.fasta] --> B0[Stage 0: Quality]
    A --> B1[Stage 1: Annotation]
    A --> B2a[Stage 2: Lifestyle]
    A --> B2b[Stage 2: Safety]

    B0 --> C0[CheckV — quality_summary.tsv]

    B1 --> C1[Pharokka — .gbk, .faa, _cds.ffn, metrics.tsv]
    C1 --> C2[PHOLD — refined .gbk, per_cds_predictions.tsv]
    C2 --> C3[Phynteny — phynteny.tsv]

    C1 --> C4[PhageDPO — depolymerase HTML+TSV]
    C1 --> C5[PhageRBPdetect v4 — RBP predictions TSV]
    C1 --> C6[PhageTailFinder — tail classification TSV]

    B2a --> D1[BACPHLIP — lifestyle .bacphlip]

    B2b --> E1[RGI — AMR hits .txt]
    B2b --> E2[Abricate CARD — AMR screen TSV]
    B2b --> E3[Abricate VFDB — virulence screen TSV]
    A --> E4[AcrDB — anti-CRISPR blastp TSV]

    C0 & C3 & C2 --> F1[Stage 3: Post-processing]
    D1 & E1 & E2 & E3 & E4 --> F1

    F1 --> G1[build_annotated_gb → _annotated.gbk]
    G1 --> G2[build_protein_table → _protein_table.tsv]
    G2 --> G3[visualize_genome → _genome_map.png/.svg]
    G2 --> G4[summarize_sample → _tool_summary.tsv, _report.txt]

    G4 --> H[Stage 4: Scoring]
    H --> I[phage_scores.tsv]
```

---

## Repository Structure

```
PhagePipeline/
├── Snakefile                    # Master workflow definition
├── config.yaml                  # Global configuration (paths, thresholds, weights)
├── tool_reg.md                  # Historical tool registry (see Known Limitations)
│
├── raw_input/                   # Input FASTA files — one per phage sample
├── held_fastas/                 # Staging area for genomes not yet in the pipeline
├── test_data/                   # Sample data (lambda_proteins.faa)
│
├── results/                     # Raw tool outputs — results/{sample}/{tool}/
│   └── {sample}/
│       ├── checkv/
│       ├── pharokka/
│       ├── phold/
│       ├── phynteny/
│       ├── phagedpo/
│       ├── phagerbpdetect/
│       ├── phagetailfinder/
│       ├── bacphlip/
│       ├── rgi/
│       ├── abricate/
│       ├── acrdb/
│       └── logs/
│
├── reports/                     # Human-facing outputs — reports/{run_id}/{sample}/
│   └── {run_id}/
│       ├── {sample}/
│       │   ├── {sample}_annotated.gbk
│       │   ├── {sample}_protein_table.tsv
│       │   ├── {sample}_genome_map.png
│       │   ├── {sample}_genome_map.svg
│       │   ├── {sample}_tool_summary.tsv
│       │   └── {sample}_report.txt
│       └── phage_scores.tsv     # Cross-sample ranking table
│
├── scripts/                     # Shell wrappers and Python post-processing
│   ├── checkv/run_checkv.sh
│   ├── pharokka/run_pharokka.sh
│   ├── phage/run_phold.sh
│   ├── phage/run_phynteny.sh
│   ├── bacphlip/run_bacphlip.sh
│   ├── rgi/run_rgi.sh
│   ├── abricate/run_abricate.sh
│   ├── acrdb/run_acrdb_blast.sh
│   ├── phagedpo/run_phagedpo.sh
│   ├── phagerbpdetect/run_phagerbpdetect_v4.sh
│   ├── phagetailfinder/run_phagetailfinder.sh
│   ├── build_annotated_gb.py
│   ├── build_protein_table.py
│   ├── visualize_genome.py
│   ├── summarize_sample.py
│   └── score_and_rank.py
│
├── databases/                   # External reference databases (not tracked in git)
│   ├── checkv-db/
│   ├── pharokka_db/
│   ├── phold_db/
│   ├── phynteny_models/
│   └── acrdb_db/
│
├── tools/                       # Vendored third-party tools (git submodules / clones)
│   ├── PhageRBPdetection/
│   ├── PhageTailFinder/
│   └── phagedpo/
│
├── configs/                     # Conda environment YAML files
│   ├── checkv.yaml
│   ├── pharokka.yml
│   └── env_phage_ml.yml
│
└── env/                         # Additional conda environment YAML files
    ├── phold_env.yml
    └── phynteny_env.yml
```

---

## Quick Start

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) / [Mamba](https://mamba.readthedocs.io/)
- [Snakemake](https://snakemake.readthedocs.io/) ≥ 7.0
- All conda environments created (see [docs/Installation.md](docs/Installation.md))
- All databases downloaded and placed under `databases/` (see [docs/Installation.md](docs/Installation.md))

### Run — Single Sample

```bash
snakemake --cores 8 \
    --config samples="PhageA" run_id="2026-07-14_PhageA"
```

### Run — Multiple Samples (Inline)

```bash
snakemake --cores 16 \
    --config samples="PhageA,PhageB,PhageC" run_id="2026-07-14_batch1"
```

### Run — Batch File

Create a plain-text file with one sample name per line (no `.fasta` extension, lines starting with `#` are ignored):

```
# batches/my_batch.txt
PhageA
PhageB
PhageC
```

```bash
snakemake --cores 16 \
    --config batch="batches/my_batch.txt" run_id="2026-07-14_batch1"
```

### Input Requirement

Place one FASTA file per sample in `raw_input/`:

```
raw_input/
├── PhageA.fasta
├── PhageB.fasta
└── PhageC.fasta
```

Each file must be a valid FASTA with a `>` header on the first line.

---

## Scoring System

Each phage candidate is scored on a 100-point scale across five dimensions:

| Dimension | Max Points | Primary Evidence Source |
|-----------|-----------|------------------------|
| Quality   | 20        | CheckV genome quality + completeness |
| Lifestyle | 30        | BACPHLIP × annotation markers (integrase/recombinase/transposase) |
| Safety    | 30        | RGI + Abricate-CARD (AMR) + Abricate-VFDB (virulence) + AcrDB (anti-CRISPR) |
| Lysis     | 10        | Named lysis proteins (holin, endolysin, lysin, spanin) × lysis category count |
| Host      | 10        | Tail, head/packaging, and connector structural proteins |
| **Total** | **100**   | |

### Status Thresholds

| Score Range | Status    | Meaning |
|-------------|-----------|---------|
| ≥ 80        | **PASS**  | Recommended for further development |
| 50 – 79     | **REVIEW**| Requires manual expert review |
| < 50        | **FAIL**  | Not recommended |
| —           | **REJECTED** | Hard disqualifier triggered (see below) |

### Hard Disqualifiers

The following automatically set score to 0 and status to **REJECTED**, regardless of other scores:

- `LifestyleConsensus = TEMPERATE` — phage shows evidence of lysogeny
- `AMRConsensus = DISQUALIFY` — antimicrobial resistance genes confirmed by RGI (Strict) and/or Abricate-CARD
- `VirulenceConsensus = DISQUALIFY` — virulence factors confirmed by Abricate-VFDB

---

## Output Files

| File | Location | Description |
|------|----------|-------------|
| `phage_scores.tsv` | `reports/{run_id}/` | Cross-sample ranking table with scores, status, and notes per dimension |
| `{sample}_report.txt` | `reports/{run_id}/{sample}/` | Human-readable candidate report with verdict statements |
| `{sample}_tool_summary.tsv` | `reports/{run_id}/{sample}/` | Machine-readable flat table of all tool metrics, values, and flags |
| `{sample}_annotated.gbk` | `reports/{run_id}/{sample}/` | Enriched GenBank file with PHOLD confidence + Phynteny synteny annotations |
| `{sample}_protein_table.tsv` | `reports/{run_id}/{sample}/` | One-row-per-CDS protein table with evidence-priority annotation |
| `{sample}_genome_map.png` | `reports/{run_id}/{sample}/` | Functional genome map (200 dpi PNG) coloured by category |
| `{sample}_genome_map.svg` | `reports/{run_id}/{sample}/` | Scalable vector version of the genome map |
| `{sample}_rbp_predictions.tsv` | `results/{sample}/phagerbpdetect/` | Per-protein receptor-binding protein predictions |
| `{sample}_rbp_summary.tsv` | `results/{sample}/phagerbpdetect/` | RBP count summary for the sample |
| `{sample}_tailfinder.tsv` | `results/{sample}/phagetailfinder/` | Per-protein tail protein classifications |

---

## Tool Summary

| Tool | Version | Purpose | Status |
|------|---------|---------|--------|
| CheckV | 1.0.1 | Genome completeness and contamination assessment | Production |
| Pharokka | 1.9.0 | Primary phage genome annotation (PHANOTATE + PHROG) | Production |
| PHOLD | latest | Structural homology-based annotation refinement (ESM-2) | Production |
| Phynteny | latest | Synteny-based prediction for hypothetical proteins | Experimental |
| BACPHLIP | latest | ML-based virulent / temperate lifestyle prediction | Production |
| RGI | latest | Antimicrobial resistance gene detection (CARD) | Production |
| Abricate | latest | Resistance and virulence gene screening (CARD, VFDB) | Production |
| AcrDB | latest | Anti-CRISPR protein detection via BLASTP | Production |
| PhageDPO | — | Depolymerase prediction (SVM-based) | Experimental |
| PhageRBPdetect v4 | — | Receptor-binding protein prediction (ESM-2 fine-tuned) | Experimental |
| PhageTailFinder | — | Tail protein classification (HMM + clustering) | Experimental |

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/Installation.md](docs/Installation.md) | Environment setup, database downloads, prerequisites |
| [docs/User_Guide.md](docs/User_Guide.md) | Running the pipeline, configuration, output interpretation |
| [docs/Developer_Guide.md](docs/Developer_Guide.md) | Codebase architecture, adding tools, scoring logic |
| [docs/Troubleshooting.md](docs/Troubleshooting.md) | Common errors and their solutions |

---

## Known Limitations

- `tool_reg.md` in the repository root was written at an earlier stage of development. The documentation in `docs/` is the current authoritative reference.
- Three tools under `tools/` (PhageRBPdetect v4, PhageTailFinder, PhageDPO) are cloned/vendored third-party projects and are marked **Experimental** in this documentation. See [docs/Troubleshooting.md](docs/Troubleshooting.md) for notes on required model/database files.
- Phynteny integration is also **Experimental**; the `phynteny_transformer` command and output format should be verified against the installed version.
- The `env_rgi`, `env_bacphlip`, `env_abricate`, and `env_acrdb` conda environments do not have YAML files in this repository. See [docs/Installation.md](docs/Installation.md) for manual creation instructions.
- BACPHLIP requires the genome FASTA to be present inside the output directory; the wrapper script handles this copy automatically.

---

*Astraphage Innovations — Internal Use*
