#!/bin/bash
set -euo pipefail
GBK_IN=${1:?Usage: run_phynteny.sh <input.gbk> <prefix> [outdir]}
PREFIX=${2:?Usage: run_phynteny.sh <input.gbk> <prefix> [outdir]}
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MODELS_PATH="$ROOT/databases/phynteny_models/models"
OUTDIR=${3:-"$ROOT/results/phynteny/${PREFIX}"}
if [ ! -f "$GBK_IN" ]; then
  echo "Input GenBank file not found: $GBK_IN" >&2
  exit 2
fi
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/esm"
conda run -n phyntenyENV \
  phynteny_transformer \
  -m "$MODELS_PATH" \
  -o "$OUTDIR" \
  -f "$GBK_IN"