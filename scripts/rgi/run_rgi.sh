#!/bin/bash
set -euo pipefail

INPUT=${1:?Usage: scripts/rgi/run_rgi.sh <input.fasta> <prefix> [contig|protein] [threads] [outdir]}
PREFIX=${2:?Usage: scripts/rgi/run_rgi.sh <input.fasta> <prefix> [contig|protein] [threads] [outdir]}
TYPE=${3:-contig}
THREADS=${4:-2}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
# 5th arg overrides default -- Snakemake passes results/{sample}/rgi/
OUTDIR=${5:-"$ROOT/results/rgi/${PREFIX}"}

if [[ "$(head -c 1 "$INPUT")" != ">" ]]; then
  echo "Input does not look like FASTA: $INPUT" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n env_rgi \
  rgi main \
  -i "$INPUT" \
  -o "$OUTDIR/${PREFIX}_rgi" \
  -t "$TYPE" \
  -a DIAMOND \
  -n "$THREADS" \
  --local \
  --clean
