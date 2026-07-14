# Troubleshooting Guide — PhagePipeline

*Astraphage Innovations · Therapeutic Phage Candidate Selection Pipeline*

---

## Table of Contents

1. [General Debugging Strategy](#1-general-debugging-strategy)
2. [Snakemake / Workflow Errors](#2-snakemake--workflow-errors)
3. [Stage 0 — CheckV Errors](#3-stage-0--checkv-errors)
4. [Stage 1 — Annotation Errors](#4-stage-1--annotation-errors)
   - 4.1 [Pharokka](#41-pharokka)
   - 4.2 [PHOLD](#42-phold)
   - 4.3 [Phynteny](#43-phynteny)
   - 4.4 [PhageDPO](#44-phagedpo)
   - 4.5 [PhageRBPdetect v4](#45-phagerbpdetect-v4)
   - 4.6 [PhageTailFinder](#46-phagetailfinder)
5. [Stage 2 — Safety / Lifestyle Errors](#5-stage-2--safety--lifestyle-errors)
   - 5.1 [BACPHLIP](#51-bacphlip)
   - 5.2 [RGI](#52-rgi)
   - 5.3 [Abricate](#53-abricate)
   - 5.4 [AcrDB](#54-acrdb)
6. [Stage 3 — Post-processing Errors](#6-stage-3--post-processing-errors)
7. [Stage 4 — Scoring Errors](#7-stage-4--scoring-errors)
8. [Conda Environment Issues](#8-conda-environment-issues)
9. [Database Issues](#9-database-issues)
10. [Parallel Execution Issues](#10-parallel-execution-issues)
11. [Output Interpretation Warnings](#11-output-interpretation-warnings)

---

## 1. General Debugging Strategy

### Step 1: Read the Log File

Every rule writes its stdout and stderr to a log file. This is the first place to look:

```
results/{sample}/logs/{tool}.log
results/logs/{run_id}_score_and_rank.log
```

```bash
cat results/MySample/logs/pharokka.log
```

### Step 2: Run with Verbose Snakemake Output

```bash
snakemake --cores 8 \
    --config samples="MySample" run_id="debug" \
    --verbose 2>&1 | tee snakemake_debug.log
```

### Step 3: Test the Wrapper in Isolation

Every shell wrapper can be run independently of Snakemake. This isolates whether the problem is in Snakemake's invocation or in the tool itself:

```bash
bash scripts/pharokka/run_pharokka.sh raw_input/MySample.fasta MySample 4
```

### Step 4: Check Conda Environment

Verify the tool is available in the expected environment:

```bash
conda run -n env_pharokka which pharokka.py
conda run -n env_phage_ml python -c "import torch; print(torch.__version__)"
```

---

## 2. Snakemake / Workflow Errors

### Error: `ValueError: No samples specified`

**Cause:** Neither `samples` nor `batch` was provided in `--config`.

**Fix:**
```bash
snakemake --cores 8 --config samples="MySample" run_id="test"
# OR
snakemake --cores 8 --config batch="batches/my_batch.txt" run_id="test"
```

---

### Error: `ValueError: Batch file not found: batches/my_batch.txt`

**Cause:** The batch file path provided does not exist relative to the project root.

**Fix:** Verify the path:
```bash
ls batches/my_batch.txt
```
Paths in `--config batch=...` are relative to the directory from which you run Snakemake.

---

### Error: `ValueError: Batch file is empty`

**Cause:** The batch file exists but contains no non-comment, non-empty lines.

**Fix:** Ensure at least one uncommented sample name exists in the file.

---

### Error: `MissingInputException` for `raw_input/{sample}.fasta`

**Cause:** The FASTA file for one or more samples is not present in `raw_input/`.

**Fix:**
```bash
ls raw_input/
# Confirm {sample}.fasta exists exactly (case-sensitive)
```

---

### Snakemake runs but produces no output, exits with no error

**Cause:** All target files already exist and Snakemake considers the run complete.

**Fix:** Use `--forcerun` to re-run specific rules, or delete the output files you want to regenerate:
```bash
snakemake --cores 8 --config samples="MySample" run_id="2026-07-14" \
    --forcerun summarize_sample score_and_rank
```

---

## 3. Stage 0 — CheckV Errors

### Error: `Input does not look like FASTA`

**Cause:** The input FASTA file does not start with `>`. This check is performed by the shell wrapper.

**Fix:** Verify the FASTA format:
```bash
head -c 1 raw_input/MySample.fasta
# Should print: >
```

---

### Error: CheckV fails with database-related messages

**Cause:** The CheckV database at `databases/checkv-db/` is missing or incomplete.

**Fix:**
```bash
conda run -n env_checkv \
    checkv download_database databases/checkv-db
```

---

### Warning flag: `checkv_quality = Low-quality` or `Not-determined` in reports

**Cause:** CheckV could not assign a quality tier. Common for fragmentary or unusual phages.

**Impact:** Quality score = 7 (Low-quality) or 0 (Not-determined). This reduces but does not automatically fail the candidate.

**Action:** Review the genome manually; verify the FASTA is complete and belongs to a single phage.

---

## 4. Stage 1 — Annotation Errors

### 4.1 Pharokka

#### Error: Pharokka fails with database path error

**Cause:** `databases/pharokka_db/` is missing or empty.

**Fix:**
```bash
conda run -n env_pharokka \
    pharokka install_databases -o databases/pharokka_db
```

---

#### Error: Expected output `{sample}.faa` or `{sample}_cds.ffn` not found

**Cause:** The `cp` normalisation steps at the end of `run_pharokka.sh` rely on Pharokka writing `phanotate.faa` and `phanotate.ffn`. If Pharokka uses a different gene caller or a newer output file naming convention, these files may not exist.

**Fix:** Inspect the Pharokka output directory:
```bash
ls results/MySample/pharokka/
```
If the files have different names, update the `cp` commands in `scripts/pharokka/run_pharokka.sh`.

---

### 4.2 PHOLD

#### Error: `Input GenBank file not found`

**Cause:** The Pharokka rule did not complete successfully, or the `.gbk` output is missing.

**Fix:** Confirm the Pharokka step completed:
```bash
ls results/MySample/pharokka/MySample.gbk
cat results/MySample/logs/pharokka.log
```

---

#### Error: PHOLD model / database errors

**Cause:** `databases/phold_db/` is missing or incomplete.

**Fix:**
```bash
conda run -n pholdENV \
    phold install_db -d databases/phold_db
```

---

#### Error: `phold.gbk` not found after PHOLD run

**Cause:** PHOLD always writes `phold.gbk` (no prefix). If the run failed silently, this file will be absent. The `cp` at the end of `run_phold.sh` will then fail.

**Fix:** Check the PHOLD log for errors:
```bash
cat results/MySample/logs/phold.log
```

---

### 4.3 Phynteny

#### Error: `phynteny_transformer: command not found`

**Cause:** The `phyntenyENV` conda environment is not correctly set up, or the `phynteny_transformer` package is not installed.

**Fix:**
```bash
conda env create -f env/phynteny_env.yml
conda run -n phyntenyENV which phynteny_transformer
```

---

#### Missing Phynteny model files

**Cause:** The Phynteny model weights are not present at `databases/phynteny_models/models/`.

**Fix:** Refer to the upstream `phynteny_transformer` documentation for the model download procedure.

---

#### Phynteny output not matching expected format

**Cause:** Phynteny is **Experimental**. Output column names or file names may differ between versions.

**Impact:** The `build_annotated_gb.py` script joins Phynteny data on `(start, end, strand)` coordinates and expects columns `phynteny_category`, `phynteny_confidence`, and `phynteny_score`. If the column names differ, Phynteny enrichment will silently be absent from the annotated GBK.

**Action:** Inspect the Phynteny output:
```bash
head -5 results/MySample/phynteny/phynteny.tsv
```

---

### 4.4 PhageDPO

#### Error: `ERROR: input FFN not found`

**Cause:** `{sample}_cds.ffn` was not produced by the Pharokka wrapper.

**Fix:** Check that Pharokka completed and that the `cp` in `run_pharokka.sh` succeeded:
```bash
ls results/MySample/pharokka/MySample_cds.ffn
```

---

#### Error: `ERROR: PhageDPO did not write expected HTML`

**Cause:** `phagedpo_cli.py` may have failed (sequence format error, parsing issue), or the filename stem does not match what PhageDPO expects.

**Fix:** Run PhageDPO manually with additional verbosity:
```bash
cat results/MySample/logs/phagedpo.log
```
PhageDPO derives output file names from the input FASTA stem. The wrapper copies `{sample}_cds.ffn` as `{sample}_cds.fasta` into the output directory, so the expected HTML is `{sample}_cds_output.html`.

---

#### Error: TSV conversion fails (pandas.read_html)

**Cause:** The HTML output from PhageDPO is empty or malformed; `pandas.read_html` finds no tables.

**Fix:** Open the HTML file manually and inspect its contents:
```bash
cat results/MySample/phagedpo/MySample_cds_output.html | head -50
```

---

### 4.5 PhageRBPdetect v4

#### Error: `ERROR: v4 model dir missing: tools/PhageRBPdetection/data/RBPdetect_v4_ESMfine`

**Cause:** The ESM-2 fine-tuned model has not been downloaded.

**Fix:** Obtain the model from the PhageRBPdetection upstream repository or internal storage:
```bash
cat tools/PhageRBPdetection/README.md
```

---

#### Error: `ERROR: predictions.csv not written`

**Cause:** The PhageRBPdetect v4 inference script failed. This may be due to CUDA errors, missing dependencies, or a corrupted model.

**Fix:**
```bash
cat results/MySample/logs/phagerbpdetect.log
```
Look for CUDA / PyTorch error messages. Verify `env_phage_ml` is correctly set up:
```bash
conda run -n env_phage_ml python -c "import torch; print(torch.cuda.is_available())"
```

---

#### Pipeline appears hung with multiple samples

**Cause:** A stale lock file exists from a previous failed run, causing all samples to block indefinitely on `flock`.

**Fix:**
```bash
rm -f tools/PhageRBPdetection/data/.phagerbp.lock
```
Then re-run the pipeline.

---

#### Silent wrong results with multiple parallel samples

**Cause:** If the lock mechanism fails (e.g., `flock` not available), multiple samples may concurrently overwrite `data/sequences.fasta` and `data/predictions.csv`.

**Verification:** Ensure `flock` is available:
```bash
which flock
```
`flock` is part of the `util-linux` package on Linux. If absent, install it:
```bash
sudo apt-get install util-linux
```

---

### 4.6 PhageTailFinder

#### Error: `ERROR: PhageTailFinder DB files missing`

**Cause:** `tools/PhageTailFinder/dbs/tail_pfam` or `tools/PhageTailFinder/dbs/nontail_pfam` does not exist.

**Fix:** Obtain the database files from the PhageTailFinder upstream repository or internal storage:
```bash
cat tools/PhageTailFinder/README.md
```

---

#### Error: per-protein output not found after PhageTailFinder run

**Cause:** `predict.py` derives output file names from the **input filename stem**. If the formatted input file is not named `{PREFIX}.faa`, the output will not be named `{PREFIX}_prot_result_table.txt`.

**Fix:** Inspect the contents of the output directory:
```bash
find results/MySample/phagetailfinder/ -type f | sort
cat results/MySample/logs/phagetailfinder.log
```

---

#### Known Limitation: False Negatives on T4 Phage

PhageTailFinder has a known false-negative rate on T4 phage (noted in `tool_reg.md`). This is a property of the tool's classifier, not a pipeline error. Use PhageTailFinder results as supplementary evidence alongside the annotation categories from PHOLD.

---

## 5. Stage 2 — Safety / Lifestyle Errors

### 5.1 BACPHLIP

#### Error: `BACPHLIP` output is `MISSING` or `EMPTY` in the report

**Cause:** BACPHLIP failed to run or produced an empty output file.

**Fix:**
```bash
cat results/MySample/logs/bacphlip.log
ls -lh results/MySample/bacphlip/
```

BACPHLIP requires:
1. The input FASTA to be copied inside the output directory (handled by `run_bacphlip.sh`)
2. `hmmsearch` available in PATH within the `env_bacphlip` environment

```bash
conda run -n env_bacphlip which hmmsearch
```

---

#### BACPHLIP output format parsing failure

The script handles two known BACPHLIP output formats (with and without a leading index column). If a different format is encountered, the BACPHLIP call will be reported as `UNKNOWN` lifestyle.

**Diagnosis:**
```bash
cat results/MySample/bacphlip/MySample.fasta.bacphlip
```
The expected format is a TSV with headers `Virulent` and `Temperate`. BACPHLIP debug output is printed to stderr (visible in the Snakemake run output and the log file).

---

### 5.2 RGI

#### Error: `CARD database not loaded` or similar

**Cause:** RGI is run with `--local`, which requires the CARD database to be loaded into the `env_rgi` environment.

**Fix:**
```bash
conda activate env_rgi
# Download and load CARD database (see Installation.md Section 4.6)
rgi load --card_json card.json --local
conda deactivate
```

---

#### Error: RGI fails with DIAMOND error

**Cause:** DIAMOND may not be installed in `env_rgi`, or the version is incompatible.

**Fix:**
```bash
conda run -n env_rgi which diamond
conda run -n env_rgi diamond --version
```

---

### 5.3 Abricate

#### Error: `Database not found`

**Cause:** The Abricate built-in databases (CARD, VFDB) are not installed or were not set up.

**Fix:**
```bash
conda run -n env_abricate abricate --setupdb
```

---

#### Abricate produces no hits but a non-zero exit code

**Cause:** Abricate writes results to stdout and may produce an empty file with a header only when no hits are found. This is normal behaviour.

**Impact:** The `summarize_sample.py` script correctly handles empty (header-only) Abricate TSVs.

---

### 5.4 AcrDB

#### Error: `BLAST: No alias or index file found`

**Cause:** The AcrDB FASTA has not been formatted with `makeblastdb`.

**Fix:**
```bash
conda run -n env_acrdb \
    makeblastdb \
    -in databases/acrdb_db/122_KnownAcr/Known_Acr.faa \
    -dbtype prot \
    -title "Known_Acr"
```

---

## 6. Stage 3 — Post-processing Errors

All post-processing scripts run in the `phage` conda environment.

#### Error: `ModuleNotFoundError: No module named 'Bio'`

**Cause:** The `phage` conda environment does not have BioPython installed.

**Fix:**
```bash
conda activate phage
pip install biopython
```

---

#### Error: `ModuleNotFoundError: No module named 'pygenomeviz'`

**Cause:** `pygenomeviz` is not installed in the `phage` environment.

**Fix:**
```bash
conda activate phage
pip install pygenomeviz
```

---

#### `build_annotated_gb.py` — Phynteny join produces 0 matches

**Cause:** The `(start, end, strand)` coordinate keys in `phynteny.tsv` do not match those in the PHOLD GBK, and the ±1 bp fuzzy match also fails.

**Diagnosis:** Check the coordinate formats:
```bash
head -5 results/MySample/phynteny/phynteny.tsv
```
Verify the expected columns are `start`, `end`, and `strand`.

**Impact:** The annotated GBK will have `phynteny_category = "no_match"` for all CDS. The pipeline will not fail; the protein table will use PHOLD annotations only.

---

#### `build_protein_table.py` — All proteins show `evidence_source = none`

**Cause:** PHOLD may have annotated all proteins as `"unknown function"` and Phynteny enrichment may have produced no matches.

**Impact:** The lysis and host scores may be 0. Review the PHOLD and Phynteny logs.

---

#### `visualize_genome.py` — `ERROR: no records found in {gbk}`

**Cause:** The annotated GBK file is empty or malformed.

**Fix:** Inspect the file:
```bash
head -20 reports/2026-07-14/MySample/MySample_annotated.gbk
```

---

#### `summarize_sample.py` — `[BACPHLIP DEBUG]` lines in log

**Cause:** This is expected. The script prints debug lines for BACPHLIP output parsing to stderr. These will appear in the `summarize_sample.log` file and are informational only.

---

## 7. Stage 4 — Scoring Errors

#### Error: `WARNING: summary not found for {sample}`

**Cause:** `score_and_rank.py` cannot find the `_tool_summary.tsv` for a sample. This means `summarize_sample` did not complete for that sample.

**Fix:** Check the summarize_sample log:
```bash
cat results/MySample/logs/summarize_sample.log
```

---

#### Sample shows `status = FAIL` with score 0 but no rejection_reasons

**Cause:** All dimension scores are 0 (e.g., CheckV quality = Not-determined, lifestyle = UNKNOWN, safety = 0 due to no hits, lysis = 0, host = 0). This is a legitimate FAIL, not a scoring bug.

**Action:** Review the `_report.txt` to understand which evidence is missing.

---

#### Changing `score_pass` / `score_review` in `config.yaml` has no effect

**Cause:** These thresholds are hardcoded in `score_and_rank.py` (lines 29–30) and are not read from `config.yaml` at runtime.

**Fix:** Update the constants directly in `scripts/score_and_rank.py`:
```python
SCORE_PASS   = 80   # change this value
SCORE_REVIEW = 50   # change this value
```

This is a known inconsistency. See the Developer Guide for details.

---

## 8. Conda Environment Issues

#### Error: `conda run -n ENV_NAME python` — environment not found

**Cause:** The conda environment has not been created.

**Fix:** Create the environment:
```bash
conda env create -f configs/checkv.yaml   # example
```
Or check which environments exist:
```bash
conda env list
```

---

#### Error: `conda: command not found` when running from Snakemake

**Cause:** Snakemake is running in a shell where conda is not initialised.

**Fix:** Ensure conda is initialised in the shell before running Snakemake:
```bash
source ~/.bashrc   # or ~/.bash_profile
conda activate base
snakemake ...
```

---

#### `env_phage_ml` installation is very slow or fails on CUDA packages

**Cause:** The `configs/env_phage_ml.yml` includes a large number of pinned CUDA and PyTorch packages. Installation may be slow and can fail if package versions are no longer available.

**Fix options:**
- Use `mamba` instead of `conda` for faster solving:
  ```bash
  mamba env create -f configs/env_phage_ml.yml
  ```
- If specific CUDA packages are unavailable, manually remove version pins from the YAML and re-attempt.

---

## 9. Database Issues

#### Databases directory is large / out of disk space

**Cause:** Reference databases for CheckV, Pharokka, and PHOLD are several gigabytes each.

**Fix:** Ensure sufficient disk space (100+ GB recommended) before starting database downloads.

---

#### CheckV completeness estimates are unreliable

**Cause:** CheckV relies on comparisons to a database of complete viral genomes. Novel phages with no close relatives may receive low or indeterminate completeness scores.

**Impact:** Quality score contribution is reduced; this is expected for novel phages and does not indicate a pipeline error.

---

## 10. Parallel Execution Issues

#### Multiple samples fail PhageRBPdetect simultaneously

**Cause:** Likely a stale lock file from a previous run.

**Fix:**
```bash
rm -f tools/PhageRBPdetection/data/.phagerbp.lock
```

---

#### Race condition on `abricate` output files

**Cause:** Both `abricate_card` and `abricate_vfdb` write to the same output directory (`results/{sample}/abricate/`) but with different filenames. Snakemake tracks these by the declared output filenames; no actual race condition should occur. If outputs are missing, check the log files.

---

## 11. Output Interpretation Warnings

### Ambiguous lifestyle (`LYTIC_WITH_MARKERS` or `AMBIGUOUS`)

A phage classified as `LYTIC_WITH_MARKERS` or `AMBIGUOUS` should be reviewed by a domain expert before advancing in the therapeutic development process. Annotation markers (integrase, recombinase, transposase keywords) can be false positives in purely lytic phages; the combination of evidence should be assessed in biological context.

### Annotation toxin keywords without VFDB confirmation

If the report shows `virulence_status = MANUAL_REVIEW` (annotation toxin keywords found, but no VFDB hit), this is a soft flag. Annotators should verify whether the flagged protein is a genuine toxin or a false-positive keyword match. A 10-point safety penalty is applied automatically.

### Low lysis score on known lytic phages

If a known lytic phage receives a lysis score of 5 (INFERRED) or 0, check the protein table to confirm that lysis proteins are present under their expected category. PHOLD/Phynteny may have annotated them under different category names, or the annotation confidence may be `low` rather than `high`, which is the threshold for a protein to be counted as a named lysis protein.

### High AcrDB weak hit count

The `weak_hits` count (e-value ≥ 1×10⁻³) carries no scoring penalty. A large number of weak hits does not indicate anti-CRISPR activity and can be safely ignored.

---

*Astraphage Innovations — Internal Use*
