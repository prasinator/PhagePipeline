#!/bin/bash
set -euo pipefail

GENOME=${1:?Usage: scripts/abricate/run_abricate.sh <genome.fasta> <prefix> [db] [outdir]}
PREFIX=${2:?Usage: scripts/abricate/run_abricate.sh <genome.fasta> <prefix> [db] [outdir]}
DB=${3:-card}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
# 4th arg overrides default -- Snakemake passes results/{sample}/abricate/
OUTDIR=${4:-"$ROOT/results/abricate"}

if [[ "$(head -c 1 "$GENOME")" != ">" ]]; then
  echo "Input does not look like FASTA: $GENOME" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n env_abricate \
  abricate \
  --db "$DB" \
  "$GENOME" \
  > "$OUTDIR/${PREFIX}_${DB}.tsv"