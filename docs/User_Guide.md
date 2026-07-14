# User Guide — PhagePipeline

*Astraphage Innovations · Therapeutic Phage Candidate Selection Pipeline*

---

## Table of Contents

1. [Overview](#1-overview)
2. [Prerequisites](#2-prerequisites)
3. [Input Requirements](#3-input-requirements)
4. [Running the Pipeline](#4-running-the-pipeline)
   - 4.1 [Single Sample](#41-single-sample)
   - 4.2 [Multiple Samples — Inline](#42-multiple-samples--inline)
   - 4.3 [Batch File Mode](#43-batch-file-mode)
   - 4.4 [Dry-Run](#44-dry-run)
   - 4.5 [Resuming a Failed Run](#45-resuming-a-failed-run)
5. [Configuration Reference](#5-configuration-reference)
6. [Understanding the Stages](#6-understanding-the-stages)
   - 6.1 [Stage 0 — Quality](#61-stage-0--quality)
   - 6.2 [Stage 1 — Annotation Chain](#62-stage-1--annotation-chain)
   - 6.3 [Stage 2 — Lifestyle and Safety](#63-stage-2--lifestyle-and-safety)
   - 6.4 [Stage 3 — Post-processing and Report Generation](#64-stage-3--post-processing-and-report-generation)
   - 6.5 [Stage 4 — Scoring and Ranking](#65-stage-4--scoring-and-ranking)
7. [Understanding the Outputs](#7-understanding-the-outputs)
   - 7.1 [phage_scores.tsv](#71-phage_scorestsv)
   - 7.2 [Sample Report (.txt)](#72-sample-report-txt)
   - 7.3 [Tool Summary (.tsv)](#73-tool-summary-tsv)
   - 7.4 [Annotated GenBank (.gbk)](#74-annotated-genbank-gbk)
   - 7.5 [Protein Table (.tsv)](#75-protein-table-tsv)
   - 7.6 [Genome Map (.png / .svg)](#76-genome-map-png--svg)
   - 7.7 [PhageRBPdetect Outputs](#77-phagerbpdetect-outputs)
   - 7.8 [PhageTailFinder Output](#78-phagetailfinder-output)
   - 7.9 [PhageDPO Output](#79-phagedpo-output)
8. [Scoring System — Detailed](#8-scoring-system--detailed)
   - 8.1 [Quality Score (max 20)](#81-quality-score-max-20)
   - 8.2 [Lifestyle Score (max 30)](#82-lifestyle-score-max-30)
   - 8.3 [Safety Score (max 30)](#83-safety-score-max-30)
   - 8.4 [Lysis Score (max 10)](#84-lysis-score-max-10)
   - 8.5 [Host Score (max 10)](#85-host-score-max-10)
   - 8.6 [Hard Disqualifiers](#86-hard-disqualifiers)
   - 8.7 [Status Classification](#87-status-classification)
9. [Flag Interpretation](#9-flag-interpretation)
10. [Parallel Execution Considerations](#10-parallel-execution-considerations)

---

## 1. Overview

PhagePipeline takes raw bacteriophage genome FASTA files as input and produces a ranked, evidence-backed report for each candidate. The pipeline integrates ten external bioinformatics tools organised into four sequential stages, followed by a deterministic scoring engine.

The two main output artefacts are:

- **`phage_scores.tsv`** — a single cross-sample table ranking all candidates in a run
- **`{sample}_report.txt`** — a human-readable narrative report for each individual phage

---

## 2. Prerequisites

- Snakemake ≥ 7.0 installed in the base or a dedicated `snakemake` conda environment
- All conda environments created (see [Installation.md](Installation.md))
- All external databases placed under `databases/` (see [Installation.md](Installation.md))
- Raw genome FASTA files present in `raw_input/`

---

## 3. Input Requirements

| Requirement | Detail |
|-------------|--------|
| Location | `raw_input/{sample}.fasta` |
| Format | Standard FASTA — first byte must be `>` |
| Naming | File stem (without `.fasta`) is the sample name used throughout the pipeline |
| Contig model | BACPHLIP is designed for single-contig, complete phage genomes. Multi-contig inputs may process but lifestyle prediction reliability is reduced. |

**Example:**

```
raw_input/
├── EcoPhage_K1.fasta
├── PAO1_lytic_phage.fasta
└── T7_wild_type.fasta
```

You would then run with `--config samples="EcoPhage_K1,PAO1_lytic_phage,T7_wild_type"`.

---

## 4. Running the Pipeline

All commands must be executed from the PhagePipeline project root directory.

### 4.1 Single Sample

```bash
snakemake --cores 8 \
    --config samples="T7_wild_type" run_id="2026-07-14_T7"
```

### 4.2 Multiple Samples — Inline

```bash
snakemake --cores 16 \
    --config samples="EcoPhage_K1,PAO1_lytic_phage,T7_wild_type" \
             run_id="2026-07-14_ecoli_batch"
```

Sample names are comma-separated with no spaces.

### 4.3 Batch File Mode

Create a plain-text batch file (lines beginning with `#` are treated as comments and skipped):

```bash
# batches/ecoli_candidates.txt
EcoPhage_K1
PAO1_lytic_phage
T7_wild_type
```

```bash
snakemake --cores 16 \
    --config batch="batches/ecoli_candidates.txt" \
             run_id="2026-07-14_ecoli_batch"
```

> **Note:** When both `batch` and `samples` are provided in `--config`, the batch file takes priority.

### 4.4 Dry-Run

Preview which rules would execute without actually running them:

```bash
snakemake --cores 8 \
    --config samples="T7_wild_type" run_id="2026-07-14_T7" \
    --dry-run
```

### 4.5 Resuming a Failed Run

Snakemake automatically tracks completed outputs. Re-running the same command will skip completed steps and resume from the point of failure:

```bash
snakemake --cores 8 \
    --config samples="T7_wild_type" run_id="2026-07-14_T7"
```

To force re-execution of specific rules:

```bash
snakemake --cores 8 \
    --config samples="T7_wild_type" run_id="2026-07-14_T7" \
    --forcerun score_and_rank
```

---

## 5. Configuration Reference

All defaults are defined in `config.yaml`. Any parameter can be overridden at the command line with `--config key=value`.

### Directory Paths

| Key | Default | Description |
|-----|---------|-------------|
| `raw_input` | `raw_input` | Directory containing input FASTA files |
| `results` | `results` | Root directory for raw tool outputs |
| `reports` | `reports` | Root directory for human-facing outputs |
| `scripts` | `scripts` | Shell wrappers and Python post-processing scripts |
| `databases` | `databases` | External reference database root |
| `tools` | `tools` | Vendored third-party tools |

### Compute

| Key | Default | Description |
|-----|---------|-------------|
| `threads` | `8` | Thread count passed to tools that support parallelism |

### AcrDB Thresholds

| Key | Default | Description |
|-----|---------|-------------|
| `acrdb_high_evalue` | `1.0e-5` | BLASTP e-value below which hits are HIGH confidence anti-CRISPR (−15 pts) |
| `acrdb_medium_evalue` | `1.0e-3` | BLASTP e-value below which hits are MEDIUM confidence (−7 pts); above this = WEAK, no penalty |

### BACPHLIP Threshold

| Key | Default | Description |
|-----|---------|-------------|
| `bacphlip_temperate_threshold` | `0.6` | Temperate probability ≥ this value triggers TEMPERATE call (hard disqualifier) |

> **Note:** The threshold value in `config.yaml` is loaded for documentation purposes. The actual decision logic is implemented directly in `scripts/summarize_sample.py` with the same value (0.6).

### Scoring Weights

| Key | Default | Dimension |
|-----|---------|-----------|
| `score_quality` | `20` | Genome quality (CheckV) |
| `score_lifestyle` | `30` | Lifestyle prediction (BACPHLIP + markers) |
| `score_safety` | `30` | Safety screening (AMR + virulence + anti-CRISPR) |
| `score_lysis` | `10` | Lysis machinery evidence |
| `score_host` | `10` | Structural host-interaction proteins |

> **Note:** The scoring weights in `config.yaml` are for reference. The scoring logic in `scripts/score_and_rank.py` uses hardcoded maximum values consistent with these settings.

### Pass/Review Thresholds

| Key | Default | Description |
|-----|---------|-------------|
| `score_pass` | `80` | Total score ≥ this → PASS |
| `score_review` | `50` | Total score ≥ this → REVIEW; below → FAIL |

### Database Paths

| Key | Default | Description |
|-----|---------|-------------|
| `acrdb_db` | `databases/acrdb_db/122_KnownAcr/Known_Acr.faa` | Path to the known anti-CRISPR protein FASTA used as a BLAST database |

---

## 6. Understanding the Stages

### 6.1 Stage 0 — Quality

**Rule:** `checkv`

CheckV performs an end-to-end quality assessment of each phage genome. It estimates completeness and contamination using a database of complete viral genomes.

| Output | Path |
|--------|------|
| `quality_summary.tsv` | `results/{sample}/checkv/` |

Key columns consumed downstream:

- `checkv_quality` — categorical quality tier
- `completeness` — percentage completeness estimate
- `contamination` — percentage host contamination estimate
- `viral_genes`, `host_genes` — gene count evidence

The CheckV database at `databases/checkv-db/` must exist before running.

---

### 6.2 Stage 1 — Annotation Chain

This stage runs in a strict dependency chain: **Pharokka → PHOLD → Phynteny**.  
PhageDPO, PhageRBPdetect, and PhageTailFinder branch off Pharokka outputs independently.

#### Pharokka

Performs primary gene calling (PHANOTATE) and functional annotation against the PHROG database.

| Output | Description |
|--------|-------------|
| `{sample}.gbk` | Annotated GenBank file — input to PHOLD |
| `{sample}.faa` | Predicted protein FASTA — input to PHOLD, AcrDB, PhageRBPdetect, PhageTailFinder |
| `{sample}_cds.ffn` | CDS nucleotide FASTA — input to PhageDPO |
| `{sample}_length_gc_cds_density.tsv` | Genome metrics (length, GC%, CDS density) |

The pharokka wrapper (`scripts/pharokka/run_pharokka.sh`) copies `phanotate.faa` → `{sample}.faa` and `phanotate.ffn` → `{sample}_cds.ffn` for consistent downstream naming.

#### PHOLD

Refines Pharokka annotations using structural homology via ESM-2 protein language model embeddings.

| Output | Description |
|--------|-------------|
| `{sample}.gbk` | PHOLD-refined GenBank file |
| `phold_per_cds_predictions.tsv` | Per-CDS annotation confidence scores |

PHOLD overwrites hypothetical protein annotations where structural homology provides functional evidence. The wrapper copies `phold.gbk` → `{sample}.gbk`.

#### Phynteny

Uses synteny context (gene neighbourhood) to predict functions for proteins that remain unannotated after PHOLD. Status: **Experimental**.

| Output | Description |
|--------|-------------|
| `phynteny.tsv` (within output directory) | Synteny-based category and confidence per CDS |

---

### 6.3 Stage 2 — Lifestyle and Safety

These rules run in parallel with each other and with Stage 1.

#### BACPHLIP (Lifestyle)

Predicts phage lifestyle using a random forest model trained on HMMER protein family profiles.

| Output | Description |
|--------|-------------|
| `{sample}.fasta.bacphlip` | Tab-separated file with Virulent and Temperate probability scores |

The wrapper (`scripts/bacphlip/run_bacphlip.sh`) copies the input FASTA into the output directory because BACPHLIP requires the input file to reside in its working directory.

**Decision logic in summarize_sample.py:**

| Temperate Probability | Call |
|-----------------------|------|
| ≥ 0.6 | Temperate |
| 0.4 – 0.59 | Ambiguous |
| < 0.4 | Virulent |

#### RGI (AMR)

Detects antimicrobial resistance genes using DIAMOND alignment against the CARD database.

| Output | Description |
|--------|-------------|
| `{sample}_rgi.txt` | Tab-separated table of AMR hits with drug class and cut-off level |

RGI is invoked in `contig` mode on the raw genome FASTA using `--local` (pre-loaded CARD database required).

#### Abricate — CARD and VFDB

Runs nucleotide-level screening against two databases:

- **CARD** — antimicrobial resistance gene families
- **VFDB** — virulence factor database

| Output | Description |
|--------|-------------|
| `{sample}_card.tsv` | CARD screening results |
| `{sample}_vfdb.tsv` | VFDB screening results |

#### AcrDB — Anti-CRISPR Detection

Performs BLASTP of Pharokka-predicted proteins against a curated set of 122 known anti-CRISPR proteins (`databases/acrdb_db/122_KnownAcr/Known_Acr.faa`).

| Output | Description |
|--------|-------------|
| `{sample}_acrdb_blastp.tsv` | BLASTP output in tabular format 6 |

Hits are classified by e-value:

| E-value Range | Confidence | Safety Penalty |
|---------------|------------|----------------|
| < 1×10⁻⁵ | HIGH | −15 pts |
| 1×10⁻⁵ to 1×10⁻³ | MEDIUM | −7 pts |
| ≥ 1×10⁻³ | WEAK | 0 pts |

---

### 6.4 Stage 3 — Post-processing and Report Generation

These rules produce the human-facing output for each sample.

#### build_annotated_gb.py

Merges the PHOLD-refined GenBank with per-CDS confidence scores and Phynteny synteny predictions. Each CDS feature gains three new qualifiers:

- `annotation_confidence` — from PHOLD (`high`, `medium`, `low`, `none`)
- `phynteny_category` — functional category from Phynteny
- `phynteny_support` — confidence tier (`high`, `medium`, `low`)

Phynteny data is joined by genomic coordinates (start, end, strand) with ±1 bp tolerance.

**Output:** `reports/{run_id}/{sample}/{sample}_annotated.gbk`

#### build_protein_table.py

Converts the annotated GenBank into a flat TSV with one row per CDS. A priority system selects the best available annotation:

| Priority | Condition | Evidence Source |
|----------|-----------|-----------------|
| 1 | PHOLD + Phynteny agree on category | `PHOLD+PHYNTENY` |
| 2 | PHOLD provides a non-unknown function | `PHOLD` |
| 3 | PHOLD unknown, Phynteny has a prediction | `PHYNTENY` |
| 4 | Neither tool provides annotation | `none` (hypothetical protein) |

**Output columns:** `protein_id`, `start`, `end`, `strand`, `product`, `category`, `confidence`, `evidence_source`

**Output:** `reports/{run_id}/{sample}/{sample}_protein_table.tsv`

#### visualize_genome.py

Generates a linear functional genome map using `pygenomeviz`. CDS features are coloured by functional category:

| Category | Colour |
|----------|--------|
| head and packaging | Blue (#4E79A7) |
| tail | Orange (#F28E2B) |
| connector | Red (#E15759) |
| lysis | Teal (#76B7B2) |
| DNA, RNA and nucleotide metabolism | Green (#59A14F) |
| transcription regulation | Yellow (#EDC948) |
| integration and excision | Purple (#B07AA1) |
| moron, auxiliary metabolic gene and host takeover | Pink (#FF9DA7) |
| other | Brown (#9C755F) |
| unknown function | Light grey (#BAB0AC) |

**Outputs:** `{sample}_genome_map.png` (200 dpi) and `{sample}_genome_map.svg`

#### summarize_sample.py

Aggregates outputs from all previous tools into a flat evidence table and a narrative report. It performs cross-tool consensus analysis:

- **LifestyleConsensus**: BACPHLIP call crossed against integrase/recombinase/transposase keywords from the protein table
- **AMRConsensus**: RGI Strict hits crossed against Abricate-CARD hits
- **VirulenceConsensus**: Abricate-VFDB hits crossed against toxin keyword matches from the protein table
- **AcrConsensus**: AcrDB e-value-tiered hits crossed against anti-CRISPR keyword matches from the protein table
- **LysisConsensus**: Named lysis proteins (holin, endolysin, lysin, spanin, anti-holin) crossed against lysis category count

**Outputs:**
- `{sample}_tool_summary.tsv` — machine-readable flat table with columns: `sample`, `tool`, `metric`, `value`, `detail`, `flag`
- `{sample}_report.txt` — human-readable narrative with section headers for Genome, Annotation, Lifestyle, Safety, Anti-CRISPR, and Lysis

---

### 6.5 Stage 4 — Scoring and Ranking

**Rule:** `score_and_rank`

Reads all `{sample}_tool_summary.tsv` files from the run and computes a numeric score for each sample. Results are sorted by status priority (PASS → REVIEW → FAIL → REJECTED) and then by descending total score within each status group.

**Output:** `reports/{run_id}/phage_scores.tsv`

---

## 7. Understanding the Outputs

### 7.1 phage_scores.tsv

The primary ranking output. Columns:

| Column | Description |
|--------|-------------|
| `sample` | Sample name |
| `status` | PASS / REVIEW / FAIL / REJECTED |
| `total_score` | Sum of all dimension scores (0 = REJECTED) |
| `quality_score` | Points from the Quality dimension (max 20) |
| `lifestyle_score` | Points from the Lifestyle dimension (max 30) |
| `safety_score` | Points from the Safety dimension (max 30) |
| `lysis_score` | Points from the Lysis dimension (max 10) |
| `host_score` | Points from the Host dimension (max 10) |
| `rejection_reasons` | Semicolon-separated list of hard disqualifiers (if any) |
| `quality_note` | Evidence notes for the quality dimension |
| `lifestyle_note` | Evidence notes for the lifestyle dimension |
| `safety_note` | Evidence notes for the safety dimension |
| `lysis_note` | Evidence notes for the lysis dimension |
| `host_note` | Evidence notes for the host dimension |

### 7.2 Sample Report (.txt)

A plain-text report with clearly delineated sections. Each section ends with a `>>> VERDICT` line indicating the consensus interpretation. Flags are displayed inline where relevant (e.g., `⚠ LOW — check gene calling`).

### 7.3 Tool Summary (.tsv)

A flat table where each row represents a single metric from a specific tool. Columns:

| Column | Description |
|--------|-------------|
| `sample` | Sample name |
| `tool` | Tool that generated the metric |
| `metric` | Metric identifier (e.g., `checkv_quality`, `bacphlip_call`) |
| `value` | Metric value |
| `detail` | Supporting detail string |
| `flag` | One of: `DISQUALIFY`, `WARN`, `ERROR`, or empty |

### 7.4 Annotated GenBank (.gbk)

Standard GenBank format. CDS features carry additional qualifiers beyond the Pharokka/PHOLD defaults:

- `annotation_confidence` — PHOLD confidence tier
- `phynteny_category` — Phynteny functional category (if matched)
- `phynteny_confidence` — Raw Phynteny confidence score
- `phynteny_score` — Phynteny raw score
- `phynteny_support` — Phynteny confidence tier (high/medium/low)

This file can be opened directly in Geneious, Artemis, SnapGene, or any GenBank-compatible viewer.

### 7.5 Protein Table (.tsv)

Columns:

| Column | Description |
|--------|-------------|
| `protein_id` | Locus tag from the GenBank feature |
| `start` | CDS start position (1-based) |
| `end` | CDS end position (1-based) |
| `strand` | `+` or `-` |
| `product` | Human-readable product name |
| `category` | Functional category (PHROG / PHOLD / Phynteny vocabulary) |
| `confidence` | Evidence confidence tier (high/medium/low/none) |
| `evidence_source` | Which tool(s) provided the annotation |

### 7.6 Genome Map (.png / .svg)

A linear genome visualisation generated by `pygenomeviz`. Features are coloured by functional category. A legend is placed below the genome track. The map title includes the phage name and genome length in bp.

The SVG version can be edited in Inkscape or Illustrator for publication-quality figures.

### 7.7 PhageRBPdetect Outputs

Two files per sample in `results/{sample}/phagerbpdetect/`:

**`{sample}_rbp_predictions.tsv`** — Per-protein predictions:

| Column | Description |
|--------|-------------|
| `protein_id` | Protein identifier |
| `rbp_prediction` | 1 = RBP, 0 = non-RBP |
| `rbp_score` | Model confidence score |
| `rbp_label` | Human-readable label (RBP / non-RBP) |

**`{sample}_rbp_summary.tsv`** — Sample-level summary:

| Column | Description |
|--------|-------------|
| `sample` | Sample name |
| `total_proteins` | Total proteins screened |
| `rbp_count` | Number predicted as RBP |
| `non_rbp_count` | Number predicted as non-RBP |
| `rbp_fraction` | Fraction of proteins predicted as RBP |

> **Note:** PhageRBPdetect v4 requires the ESM-2 fine-tuned model directory `tools/PhageRBPdetection/data/RBPdetect_v4_ESMfine`. If absent, the rule will fail at runtime. Inference is serialised via `flock` to prevent data corruption when running multiple samples in parallel.

### 7.8 PhageTailFinder Output

`results/{sample}/phagetailfinder/{sample}_tailfinder.tsv`

| Column | Description |
|--------|-------------|
| `phage_content` | PhageTailFinder phage content classification |
| `phage_id` | Phage identifier derived from filename |
| `protein_id` | Protein identifier |
| `protein_size` | Protein length |
| `protein_count` | Protein count in cluster |
| `start_index` | Cluster start index |
| `end_index` | Cluster end index |
| `protein_def` | Protein family definition |
| `tail_or_not` | Classification flag (1 = tail, 0 = non-tail) |
| `sample` | Sample name |
| `is_tail` | Boolean derived from `tail_or_not` |

> **Note:** PhageTailFinder requires two Pfam database files (`tools/PhageTailFinder/dbs/tail_pfam` and `tools/PhageTailFinder/dbs/nontail_pfam`). If absent, the rule will fail at runtime. A known false-negative rate on T4 phage has been observed; use as supplementary evidence.

### 7.9 PhageDPO Output

Two files per sample in `results/{sample}/phagedpo/`:

- **`{sample}_cds_output.html`** — Interactive HTML report from PhageDPO (SVM-based depolymerase prediction)
- **`{sample}_cds_output.tsv`** — Tabular version of the HTML report, extracted using `pandas.read_html()`

> **Note:** PhageDPO is an Experimental tool. Accurate true-negative results have been observed on T4 and T7 validation data. The tool operates on CDS nucleotide FASTA (`{sample}_cds.ffn`). The wrapper copies the file into the output directory because PhageDPO derives output file names from the directory it scans.

---

## 8. Scoring System — Detailed

### 8.1 Quality Score (max 20)

| CheckV Quality | Base Points |
|----------------|-------------|
| Complete | 20 |
| High-quality | 20 |
| Medium-quality | 14 |
| Low-quality | 7 |
| Not-determined | 0 |

**Penalties:**
- Contamination > 5%: −5 pts (minimum score: 0)

### 8.2 Lifestyle Score (max 30)

| Consensus | Confidence | Points | Disqualifier |
|-----------|------------|--------|--------------|
| TEMPERATE | Any | 0 | **YES** |
| LYTIC | HIGH | 30 | No |
| LYTIC_WITH_MARKERS | LOW | 15 | No |
| AMBIGUOUS | LOW | 15 | No |
| UNKNOWN | NONE | 0 | No |

**Consensus logic:**

| BACPHLIP | Annotation Markers | Consensus | Confidence |
|----------|--------------------|-----------|------------|
| Temperate | Any | TEMPERATE | HIGH (if markers) / MEDIUM (if no markers) |
| Virulent | None | LYTIC | HIGH |
| Virulent | Present | LYTIC_WITH_MARKERS | LOW |
| Ambiguous | Any | AMBIGUOUS | LOW |

Annotation markers = proteins matching integrase, recombinase, or transposase keyword sets.

### 8.3 Safety Score (max 30)

Starting at 30 pts, the following deductions and disqualifiers apply:

**AMR:**
- RGI Strict hit(s) **or** Abricate-CARD hit(s) → **DISQUALIFY** (score = 0, status = REJECTED)
- Both RGI and CARD confirm AMR → note: "CONFIRMED by RGI+CARD"

**Virulence:**
- Abricate-VFDB hit(s) → **DISQUALIFY**
- Annotation toxin keywords only (no VFDB) → MANUAL_REVIEW, −10 pts

**Anti-CRISPR (AcrDB):**
- HIGH confidence hit(s) (e-value < 1×10⁻⁵): −15 pts
- MEDIUM confidence hit(s) (1×10⁻⁵ ≤ e-value < 1×10⁻³): −7 pts
- WEAK hits: no penalty
- These are penalties only, not disqualifiers

Minimum safety score after penalties: 0.

### 8.4 Lysis Score (max 10)

| Evidence | Points |
|----------|--------|
| Named lysis protein(s) detected (holin, endolysin, lysin, spanin, anti-holin) | 10 |
| Lysis category CDS present but no named proteins | 5 |
| No lysis evidence detected | 0 |

Only `high`-confidence named lysis proteins (from the protein table) are counted toward named protein detection.

### 8.5 Host Score (max 10)

| Category Present | Points |
|-----------------|--------|
| Tail proteins (category = "tail") | +5 |
| Head and packaging proteins | +3 |
| Connector proteins | +2 |

Maximum total capped at 10.

### 8.6 Hard Disqualifiers

When any of the following is true, total score is set to 0 and status to **REJECTED**:

1. `LifestyleConsensus.consensus` = `TEMPERATE`
2. `AMRConsensus.amr_status` = `DISQUALIFY`
3. `VirulenceConsensus.virulence_status` = `DISQUALIFY`

Multiple disqualifiers may apply simultaneously and are all reported in the `rejection_reasons` column.

### 8.7 Status Classification

After scoring:

| Condition | Status |
|-----------|--------|
| Any hard disqualifier | REJECTED |
| Total ≥ 80 | PASS |
| Total ≥ 50 | REVIEW |
| Total < 50 | FAIL |

---

## 9. Flag Interpretation

Flags appear in the `flag` column of `{sample}_tool_summary.tsv` and drive the scoring logic:

| Flag | Meaning | Action |
|------|---------|--------|
| `DISQUALIFY` | Hard disqualifier triggered | Score → 0, Status → REJECTED |
| `WARN` | Soft signal requiring attention | Does not automatically reject; influences scoring |
| `ERROR` | A tool output was missing or empty | Review logs; the pipeline may have failed for this tool |
| *(empty)* | Normal, no concern | — |

---

## 10. Parallel Execution Considerations

- Most pipeline stages run fully in parallel across samples and tools.
- **PhageRBPdetect v4** uses `flock` (file-based locking) to serialise inference across samples. This is because the inference script writes to hardcoded shared files (`tools/PhageRBPdetection/data/sequences.fasta` and `data/predictions.csv`). Samples queue on the lock; only the inference step itself is serialised — pre- and post-processing remain parallel.
- If you observe a hung pipeline with multiple samples, check that no stale `.phagerbp.lock` file exists in `tools/PhageRBPdetection/data/`.

---

*Astraphage Innovations — Internal Use*
