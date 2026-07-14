#!/bin/bash
set -euo pipefail

GENOME=${1:?Usage: scripts/pharokka/run_pharokka.sh <genome.fasta> <prefix> [threads] [outdir]}
PREFIX=${2:?Usage: scripts/pharokka/run_pharokka.sh <genome.fasta> <prefix> [threads] [outdir]}
THREADS=${3:-2}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
DB="$ROOT/databases/pharokka_db"
# 4th arg overrides default -- Snakemake passes results/{sample}/pharokka/
OUTDIR=${4:-"$ROOT/results/pharokka/${PREFIX}"}

if [[ "$(head -c 1 "$GENOME")" != ">" ]]; then
  echo "Input does not look like FASTA: $GENOME" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n env_pharokka \
  pharokka.py \
  -i "$GENOME" \
  -o "$OUTDIR" \
  -d "$DB" \
  -t "$THREADS" \
  -p "$PREFIX" \
  -f

# Normalize pharokka's fixed-name outputs to prefix-named files so downstream
# Snakemake rules can reference them as {sample}.faa / {sample}_cds.ffn.
# phanotate.faa  -> {PREFIX}.faa   (used by: acrdb, phagedpo trigger)
# phanotate.ffn  -> {PREFIX}_cds.ffn  (used by: phagedpo as input;
#   _cds suffix is intentional -- phagedpo derives output name from input
#   stem: {PREFIX}_cds.ffn → {PREFIX}_cds_output.html, matching Snakemake)
cp "$OUTDIR/phanotate.faa" "$OUTDIR/${PREFIX}.faa"
cp "$OUTDIR/phanotate.ffn" "$OUTDIR/${PREFIX}_cds.ffn"