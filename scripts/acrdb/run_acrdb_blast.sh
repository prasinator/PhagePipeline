#!/bin/bash
set -euo pipefail

PROTEINS=${1:?Usage: scripts/acrdb/run_acrdb_blast.sh <proteins.faa> <prefix> [database-prefix] [outdir]}
PREFIX=${2:?Usage: scripts/acrdb/run_acrdb_blast.sh <proteins.faa> <prefix> [database-prefix] [outdir]}

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
DB=${3:-"$ROOT/databases/acrdb_db/122_KnownAcr/Known_Acr.faa"}
# 4th arg overrides default -- Snakemake passes results/{sample}/acrdb/
OUTDIR=${4:-"$ROOT/results/acrdb"}

if [[ "$(head -c 1 "$PROTEINS")" != ">" ]]; then
  echo "Input does not look like FASTA: $PROTEINS" >&2
  exit 2
fi

mkdir -p "$OUTDIR"

conda run -n env_acrdb \
  blastp \
  -query "$PROTEINS" \
  -db "$DB" \
  -out "$OUTDIR/${PREFIX}_acrdb_blastp.tsv" \
  -outfmt 6 \
  -num_threads 2