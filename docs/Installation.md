# Installation Guide — PhagePipeline

*Astraphage Innovations · Therapeutic Phage Candidate Selection Pipeline*

---

## Table of Contents

1. [System Requirements](#1-system-requirements)
2. [Clone the Repository](#2-clone-the-repository)
3. [Conda Environment Setup](#3-conda-environment-setup)
   - 3.1 [Environments with YAML Files](#31-environments-with-yaml-files)
   - 3.2 [Environments without YAML Files](#32-environments-without-yaml-files)
   - 3.3 [Environment Summary Table](#33-environment-summary-table)
4. [Database Setup](#4-database-setup)
   - 4.1 [CheckV Database](#41-checkv-database)
   - 4.2 [Pharokka Database](#42-pharokka-database)
   - 4.3 [PHOLD Database](#43-phold-database)
   - 4.4 [Phynteny Models](#44-phynteny-models)
   - 4.5 [AcrDB Known Anti-CRISPR Database](#45-acrdb-known-anti-crispr-database)
   - 4.6 [RGI / CARD Database](#46-rgi--card-database)
5. [Vendored Tool Setup](#5-vendored-tool-setup)
   - 5.1 [PhageRBPdetect v4 — Model Download](#51-phagerbpdetect-v4--model-download)
   - 5.2 [PhageTailFinder — Pfam Database](#52-phagetailfinder--pfam-database)
6. [Verifying the Installation](#6-verifying-the-installation)

---

## 1. System Requirements

| Requirement | Minimum | Recommended |
|-------------|---------|-------------|
| OS | Linux (x86_64) | Ubuntu 22.04 LTS |
| CPU cores | 4 | 8–16 |
| RAM | 16 GB | 32–64 GB |
| Storage | 50 GB free | 100+ GB (databases are large) |
| GPU | Not required for core pipeline | CUDA-compatible GPU for faster ML tool inference |
| Conda / Mamba | conda ≥ 23.x or mamba ≥ 1.5 | mamba (faster solver) |
| Snakemake | ≥ 7.0 | ≥ 7.32 |
| Git | ≥ 2.30 | — |
| `flock` | Required (Linux standard) | — |

> **Note:** The pipeline has not been tested on macOS or Windows. The `flock` command required for PhageRBPdetect concurrency safety is a Linux-only utility.

---

## 2. Clone the Repository

```bash
git clone https://github.com/astraphage/PhagePipeline.git
cd PhagePipeline
```

The `tools/` directory contains vendored third-party repositories. These are included as standalone git repositories (not as submodules in the current configuration). Verify they are present:

```bash
ls tools/
# Expected: PhageRBPdetection  PhageTailFinder  phagedpo
```

If any are missing, they must be cloned separately. Refer to the relevant upstream repositories.

---

## 3. Conda Environment Setup

The pipeline uses multiple separate conda environments. Some have reproducible YAML files in this repository; others must be created manually.

### 3.1 Environments with YAML Files

The following environments can be created from YAML files tracked in the repository:

#### env_checkv

```bash
conda env create -f configs/checkv.yaml
```

| YAML path | `configs/checkv.yaml` |
|-----------|----------------------|
| Env name | `env_checkv` |
| Key packages | CheckV 1.0.1, DIAMOND 2.2.0, HMMER 3.3 |

---

#### env_pharokka

```bash
conda env create -f configs/pharokka.yml
```

| YAML path | `configs/pharokka.yml` |
|-----------|----------------------|
| Env name | `env_pharokka` |
| Key packages | Pharokka 1.9.0, Python 3.10 |

---

#### pholdENV

```bash
conda env create -f env/phold_env.yml
```

| YAML path | `env/phold_env.yml` |
|-----------|----------------------|
| Env name | `pholdENV` |
| Key packages | phold (latest at time of install), Python 3.11 |

---

#### phyntenyENV

```bash
conda env create -f env/phynteny_env.yml
```

| YAML path | `env/phynteny_env.yml` |
|-----------|----------------------|
| Env name | `phyntenyENV` |
| Key packages | phynteny_transformer (latest at time of install), Python 3.10 |

---

#### env_phage_ml

This is a fully pinned environment for the three ML-based vendored tools (PhageDPO, PhageRBPdetect v4, PhageTailFinder) and all Python post-processing scripts that require PyTorch and Hugging Face.

```bash
conda env create -f configs/env_phage_ml.yml
```

| YAML path | `configs/env_phage_ml.yml` |
|-----------|---------------------------|
| Env name | `env_phage_ml` |
| Key packages | PyTorch 2.12.0, Transformers 5.10.2, XGBoost 3.2.0, scikit-learn 1.7.2, BioPython 1.87, pandas 2.3.3, matplotlib 3.10.9 |
| GPU support | NVIDIA CUDA 13 via pip packages (`cuda-toolkit`, `nvidia-cudnn-cu13`, etc.) |

> **Note:** The `env_phage_ml.yml` is a fully pinned lock-file-style export. If you are installing on a machine without a CUDA-compatible GPU, the CUDA packages will still install but GPU acceleration will not be available. The pipeline will use CPU inference in that case (significantly slower for PhageRBPdetect v4).

---

### 3.2 Environments without YAML Files

The following environments are referenced in shell wrappers but **do not have YAML files in this repository**. They must be created manually. The exact package versions to install are not specified here because the repository does not contain this information — consult the tool's official documentation for current installation instructions.

#### env_bacphlip

Used by: `scripts/bacphlip/run_bacphlip.sh`

```bash
conda create -n env_bacphlip python=3.x
conda activate env_bacphlip
# Install bacphlip per its official documentation
# https://github.com/adamhockenberry/bacphlip
```

BACPHLIP also requires `hmmsearch` to be available in PATH. The wrapper locates it via `which hmmsearch` and passes the full path using `--local_hmmsearch`.

---

#### env_rgi

Used by: `scripts/rgi/run_rgi.sh`

The pipeline invokes RGI with DIAMOND alignment in `--local` mode, which requires the CARD database to be pre-loaded. See also [Section 4.6 — RGI / CARD Database](#46-rgi--card-database).

```bash
conda create -n env_rgi -c conda-forge -c bioconda rgi
# Then load the CARD database (see Section 4.6)
```

---

#### env_abricate

Used by: `scripts/abricate/run_abricate.sh`

```bash
conda create -n env_abricate -c conda-forge -c bioconda abricate
```

Abricate uses its own bundled databases (CARD, VFDB, etc.) accessed via `--db`. Ensure these are available or run `abricate --setupdb` after installation.

---

#### env_acrdb

Used by: `scripts/acrdb/run_acrdb_blast.sh`

The AcrDB step uses BLASTP against a FASTA database. The environment requires BLAST+ to be installed.

```bash
conda create -n env_acrdb -c conda-forge -c bioconda blast
```

The FASTA database at `databases/acrdb_db/122_KnownAcr/Known_Acr.faa` must be formatted with `makeblastdb` before the first run (see [Section 4.5](#45-acrdb-known-anti-crispr-database)).

---

#### phage (Python post-processing environment)

Used by: All Python post-processing rules in the Snakefile (`build_annotated_gb`, `build_protein_table`, `visualize_genome`, `summarize_sample`, `score_and_rank`)

This environment requires at minimum:

- Python 3.x
- BioPython
- pygenomeviz
- matplotlib

```bash
conda create -n phage python=3.10
conda activate phage
pip install biopython pygenomeviz matplotlib
```

> **Important:** No YAML file for this environment exists in the repository. The exact package versions required are not specified. If you encounter import errors, refer to each Python script's import statements to determine the required packages.

---

### 3.3 Environment Summary Table

| Env Name | YAML File Exists | YAML Path | Used By |
|----------|-----------------|-----------|---------|
| `env_checkv` | ✅ | `configs/checkv.yaml` | `run_checkv.sh` |
| `env_pharokka` | ✅ | `configs/pharokka.yml` | `run_pharokka.sh` |
| `pholdENV` | ✅ | `env/phold_env.yml` | `run_phold.sh` |
| `phyntenyENV` | ✅ | `env/phynteny_env.yml` | `run_phynteny.sh` |
| `env_phage_ml` | ✅ | `configs/env_phage_ml.yml` | `run_phagedpo.sh`, `run_phagerbpdetect_v4.sh`, `run_phagetailfinder.sh` |
| `env_bacphlip` | ❌ | — | `run_bacphlip.sh` |
| `env_rgi` | ❌ | — | `run_rgi.sh` |
| `env_abricate` | ❌ | — | `run_abricate.sh` |
| `env_acrdb` | ❌ | — | `run_acrdb_blast.sh` |
| `phage` | ❌ | — | `build_annotated_gb.py`, `build_protein_table.py`, `visualize_genome.py`, `summarize_sample.py`, `score_and_rank.py` |

---

## 4. Database Setup

All databases must be placed under the `databases/` directory. This directory is excluded from git tracking (see `.gitignore`).

### 4.1 CheckV Database

The `run_checkv.sh` wrapper expects the database at `databases/checkv-db/`.

```bash
conda run -n env_checkv \
    checkv download_database databases/checkv-db
```

Alternatively:

```bash
conda activate env_checkv
checkv download_database databases/checkv-db
```

Expected directory after download: `databases/checkv-db/` containing HMM and DIAMOND database files.

---

### 4.2 Pharokka Database

The `run_pharokka.sh` wrapper expects the database at `databases/pharokka_db/`.

```bash
conda run -n env_pharokka \
    pharokka install_databases -o databases/pharokka_db
```

---

### 4.3 PHOLD Database

The `run_phold.sh` wrapper expects the database at `databases/phold_db/`.

```bash
conda run -n pholdENV \
    phold install_db -d databases/phold_db
```

---

### 4.4 Phynteny Models

The `run_phynteny.sh` wrapper expects the model weights at `databases/phynteny_models/models/`.

```bash
# Refer to the phynteny_transformer documentation for the current
# model download procedure. The expected path is:
# databases/phynteny_models/models/
```

> **Note:** The exact download command for Phynteny model weights is not specified in this repository. Refer to the upstream `phynteny_transformer` documentation.

---

### 4.5 AcrDB Known Anti-CRISPR Database

The pipeline uses the 122 known anti-CRISPR proteins FASTA as a BLAST database. The raw FASTA file is present at `databases/acrdb_db/122_KnownAcr/Known_Acr.faa`. Before the first run, it must be formatted with `makeblastdb`:

```bash
conda run -n env_acrdb \
    makeblastdb \
    -in databases/acrdb_db/122_KnownAcr/Known_Acr.faa \
    -dbtype prot \
    -title "Known_Acr"
```

This creates the BLAST index files (`.phr`, `.pin`, `.psq`, etc.) in the same directory. The BLAST step in the pipeline will fail if these files are absent.

> **Note:** A larger `acr_mge_db` BLAST database is also present in `databases/acrdb_db/` but is **not** referenced by the current pipeline. Only the `Known_Acr.faa` database (122 proteins) is used.

---

### 4.6 RGI / CARD Database

RGI is run in `--local` mode which requires the CARD database to be pre-loaded into the `env_rgi` environment:

```bash
conda activate env_rgi
# Download the CARD database
wget https://card.mcmaster.ca/latest/data -O card.tar.bz2
tar -xjf card.tar.bz2
rgi load --card_json card.json --local
conda deactivate
```

The exact download URL and CARD database version to use are not specified in this repository. Refer to the [CARD RGI documentation](https://github.com/arpcard/rgi) for current instructions.

---

## 5. Vendored Tool Setup

### 5.1 PhageRBPdetect v4 — Model Download

The ESM-2 fine-tuned model is required at:

```
tools/PhageRBPdetection/data/RBPdetect_v4_ESMfine/
```

This model directory is **not included in the repository** (too large for git). Obtain it from the PhageRBPdetection upstream repository or from internal Astraphage storage:

```bash
# Check the tools/PhageRBPdetection/README.md for current download instructions
cat tools/PhageRBPdetection/README.md
```

After downloading, verify the directory exists:

```bash
ls tools/PhageRBPdetection/data/RBPdetect_v4_ESMfine/
```

If the directory is absent, the `phagerbpdetect` rule will fail at runtime.

---

### 5.2 PhageTailFinder — Pfam Database

PhageTailFinder requires two Pfam database files:

```
tools/PhageTailFinder/dbs/tail_pfam
tools/PhageTailFinder/dbs/nontail_pfam
```

These files are **not included in the repository**. Obtain them from the PhageTailFinder upstream repository or from internal Astraphage storage:

```bash
cat tools/PhageTailFinder/README.md
```

After obtaining the files, verify they are present:

```bash
ls -lh tools/PhageTailFinder/dbs/
```

If the files are absent, the `phagetailfinder` rule will fail at the validation step with an explicit error message.

---

## 6. Verifying the Installation

After completing all setup steps, verify the installation with a dry-run. You will need at least one FASTA file in `raw_input/`:

```bash
# Place a test genome
cp test_data/lambda_proteins.faa raw_input/  # Note: this is a protein FASTA, not a genome
# For a real verification, place a complete genome FASTA in raw_input/

# Dry-run to check the DAG resolves
snakemake --cores 1 \
    --config samples="MyPhage" run_id="install_test" \
    --dry-run
```

Check that:
- The dry-run completes without errors
- The expected rules appear in the dry-run output
- No `KeyError` or `ValueError` is raised by the Snakefile

To verify individual conda environments are correctly set up:

```bash
conda run -n env_checkv checkv --version
conda run -n env_pharokka pharokka.py --version
conda run -n pholdENV phold --version
conda run -n env_phage_ml python -c "import torch; print(torch.__version__)"
```

---

*Astraphage Innovations — Internal Use*
