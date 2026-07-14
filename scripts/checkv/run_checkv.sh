#!/bin/bash
set -euo pipefail

GENOME=${1:?Usage: scripts/checkv/run_checkv.sh <genome.fasta> <prefix> [threads] [outdir]}
PREFIX=${2:?Usage: scripts/checkv/run_checkv.sh <genome.fasta> <prefix> [threads] [outdir]}
THREADS=${3:-2}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
DB="$ROOT/databases/checkv-db"
# 4th arg overrides default -- Snakemake passes results/{sample}/checkv/
OUTDIR=${4:-"$ROOT/results/checkv/${PREFIX}"}

if [[ "$(head -c 1 "$GENOME")" != ">" ]]; then
  echo "Input does not look like FASTA: $GENOME" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n env_checkv \
  checkv end_to_end "$GENOME" "$OUTDIR" \
  -d "$DB" \
  -t "$THREADS"
