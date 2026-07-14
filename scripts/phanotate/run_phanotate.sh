#!/bin/bash
set -euo pipefail

GENOME=${1:?Usage: scripts/phanotate/run_phanotate.sh <genome.fasta> <prefix>}
PREFIX=${2:?Usage: scripts/phanotate/run_phanotate.sh <genome.fasta> <prefix>}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
OUTDIR="$ROOT/results/phanotate"

if [[ "$(head -c 1 "$GENOME")" != ">" ]]; then
  echo "Input does not look like FASTA: $GENOME" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n env_pharokka \
  phanotate.py "$GENOME" \
  -o "$OUTDIR/${PREFIX}_phanotate.tsv" \
  -f tabular

conda run -n env_pharokka \
  phanotate.py "$GENOME" \
  -o "$OUTDIR/${PREFIX}_phanotate.faa" \
  -f faa
