#!/bin/bash
set -euo pipefail

GENOME=${1:?Usage: scripts/bacphlip/run_bacphlip.sh <single-contig-genome.fasta> <prefix> [outdir]}
PREFIX=${2:?Usage: scripts/bacphlip/run_bacphlip.sh <single-contig-genome.fasta> <prefix> [outdir]}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
HMMSEARCH=$(which hmmsearch)
# 3rd arg overrides default -- Snakemake passes results/{sample}/bacphlip/
OUTDIR=${3:-"$ROOT/results/bacphlip/${PREFIX}"}
INPUT_COPY="$OUTDIR/${PREFIX}.fasta"

if [[ "$(head -c 1 "$GENOME")" != ">" ]]; then
  echo "Input does not look like FASTA: $GENOME" >&2
  exit 2
fi

mkdir -p "$OUTDIR"
# BACPHLIP requires input fasta to be inside the output directory
cp "$GENOME" "$INPUT_COPY"

conda run -n env_bacphlip \
  bacphlip \
  -i "$INPUT_COPY" \
  -f \
  --local_hmmsearch "$HMMSEARCH"