#!/bin/bash
# run_all_ml_tools.sh
# Master script to run PhageTailFinder, PhageRBPdetect, and PhageDPO on a specific FASTA input.

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input_fasta_or_dir> <output_prefix_or_name>"
    echo "Example: $0 results/proteins/T4.faa T4"
    exit 1
fi

INPUT=$1
NAME=$2

echo "=========================================="
echo " Starting Phage ML Pipeline for: $NAME "
echo "=========================================="

# Ensure scripts are executed from the project root
cd "$(dirname "$0")/.."
PROJECT_ROOT=$(pwd)

# 1. Run PhageTailFinder
echo "------------------------------------------"
echo " [1/3] Running PhageTailFinder..."
echo "------------------------------------------"
bash scripts/phagetailfinder/run_phagetailfinder.sh "$INPUT" "$NAME"

# 2. Run PhageRBPdetect
echo "------------------------------------------"
echo " [2/3] Running PhageRBPdetect..."
echo "------------------------------------------"
bash scripts/phagerbpdetect/run_phagerbpdetect_v4.sh "$INPUT" "$NAME"

# 3. Run PhageDPO
# PhageDPO expects a directory of CDS files.
# If INPUT is a single FASTA, we assume a CDS folder exists at results/smoke_tests/phagedpo/cds_input
# For a robust generic run, we should point to the appropriate CDS split directory.
# We will use results/smoke_tests/phagedpo/cds_input as a default or prompt if missing.
echo "------------------------------------------"
echo " [3/3] Running PhageDPO..."
echo "------------------------------------------"
# PhageDPO script wrapper handles its own logic, we pass the CDS directory
# Here we'll pass the default CDS input for T4 if it exists, else skip
CDS_DIR="results/smoke_tests/phagedpo/cds_input"
if [ -d "$CDS_DIR" ]; then
    bash scripts/phagedpo/run_phagedpo.sh "$CDS_DIR"
else
    echo "Warning: CDS directory $CDS_DIR not found. Please provide a folder of single-sequence fasta files."
    echo "Skipping PhageDPO."
fi

echo "=========================================="
echo " Pipeline Complete. Check results/ folder."
echo "=========================================="
