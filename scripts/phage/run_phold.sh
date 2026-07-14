#!/bin/bash
set -euo pipefail

GBK_IN=${1:?Usage: scripts/phold/run_phold.sh <input.gbk> <prefix> [outdir]}
PREFIX=${2:?Usage: scripts/phold/run_phold.sh <input.gbk> <prefix> [outdir]}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
DB_PATH="$ROOT/databases/phold_db"
# 3rd arg overrides default -- Snakemake passes results/{sample}/phold/
OUTDIR=${3:-"$ROOT/results/phold/${PREFIX}"}

if [ ! -f "$GBK_IN" ]; then
  echo "Input GenBank file not found: $GBK_IN" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n pholdENV \
  phold run \
  -i "$GBK_IN" \
  -o "$OUTDIR" \
  -d "$DB_PATH" \
  -t 8 \
  -f

# phold always writes phold.gbk (no prefix flag). Copy to {PREFIX}.gbk so
# Snakemake and phynteny can find it as {sample}.gbk.
# Verified: phold.gbk present in results/phold/T7/ and results/phold/Lambda/
cp "$OUTDIR/phold.gbk" "$OUTDIR/${PREFIX}.gbk"