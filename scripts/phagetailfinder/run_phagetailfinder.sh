#!/bin/bash
set -euo pipefail

# =============================================================================
# run_phagetailfinder.sh
#
# Usage:
#   run_phagetailfinder.sh <proteins.faa> <prefix> <outdir>
#
# predict.py derives phageid from the INPUT FILENAME STEM, not from headers.
# So the input faa must be named {PREFIX}.faa for outputs to be named
# {PREFIX}_prot_result_table.txt.
#
# Output locations (testone mode):
#   {outdir}/{prefix}_prot_result_table.txt        ← per-protein (main)
#   {outdir}/each_phage_result/{prefix}_result_table.txt ← per-phage summary
# =============================================================================

PROTEINS=${1:?Usage: run_phagetailfinder.sh <proteins.faa> <prefix> <outdir>}
PROTEINS=$(realpath "$PROTEINS")
PREFIX=${2:?Usage: run_phagetailfinder.sh <proteins.faa> <prefix> <outdir>}
OUTDIR=${3:?Usage: run_phagetailfinder.sh <proteins.faa> <prefix> <outdir>}
OUTDIR=$(realpath -m "$OUTDIR")

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
TOOL="$ROOT/tools/PhageTailFinder"
FORMAT_SCRIPT="$ROOT/scripts/phagetailfinder/format_fasta.py"

# Input to predict.py must be named {PREFIX}.faa so phageid = PREFIX
FORMATTED="$OUTDIR/${PREFIX}.faa"

PROT_OUT="$OUTDIR/${PREFIX}_prot_result_table.txt"
TSV_OUT="$OUTDIR/${PREFIX}_tailfinder.tsv"

# ── validate ──────────────────────────────────────────────────────────────────
if [ ! -f "$PROTEINS" ]; then
    echo "ERROR: input .faa not found: $PROTEINS" >&2
    exit 2
fi
if [[ "$(head -c 1 "$PROTEINS")" != ">" ]]; then
    echo "ERROR: input does not look like FASTA: $PROTEINS" >&2
    exit 2
fi
if [ ! -f "$TOOL/dbs/tail_pfam" ] || [ ! -f "$TOOL/dbs/nontail_pfam" ]; then
    echo "ERROR: PhageTailFinder DB files missing in $TOOL/dbs/" >&2
    exit 2
fi

mkdir -p "$OUTDIR"

# ── format FASTA headers ──────────────────────────────────────────────────────
# Output file is named {PREFIX}.faa so predict.py extracts phageid = PREFIX
# from the filename stem → output: {PREFIX}_prot_result_table.txt
echo "[$PREFIX] Formatting FASTA headers..."
conda run -n env_phage_ml \
    python "$FORMAT_SCRIPT" "$PROTEINS" "$FORMATTED" "$PREFIX"

if [ ! -f "$FORMATTED" ]; then
    echo "ERROR: format_fasta.py did not produce: $FORMATTED" >&2
    exit 2
fi

# ── run PhageTailFinder ───────────────────────────────────────────────────────
# MUST cd to scripts/ — predict.py uses os.path.abspath("..") for model/db
echo "[$PREFIX] Running PhageTailFinder..."
cd "$TOOL/scripts"
conda run -n env_phage_ml \
    python predict.py \
        -i "$FORMATTED" \
        -o "$OUTDIR"

# ── verify output ─────────────────────────────────────────────────────────────
if [ ! -f "$PROT_OUT" ]; then
    echo "ERROR: per-protein output not found: $PROT_OUT" >&2
    echo "       Contents of $OUTDIR:" >&2
    find "$OUTDIR" -type f | sort >&2
    exit 2
fi

# ── convert to pipeline-standard TSV ─────────────────────────────────────────
# ── convert to pipeline-standard TSV ─────────────────────────────────────────
CONVERT_SCRIPT=$(mktemp /tmp/phagetailfinder_convert_XXXXXX.py)
cat > "$CONVERT_SCRIPT" <<'PYEOF'
import sys
import pandas as pd
prot_src = sys.argv[1]
tsv_dst  = sys.argv[2]
sample   = sys.argv[3]
df = pd.read_csv(prot_src, sep="\t")
df = df.rename(columns={
    "PhageContent":      "phage_content",
    "PhageId":           "phage_id",
    "ProteinId":         "protein_id",
    "ProteinSize":       "protein_size",
    "ProteinCount":      "protein_count",
    "ProteinStartIndex": "start_index",
    "ProteinEndIndex":   "end_index",
    "ProteinDef":        "protein_def",
    "TailOrNotTail":     "tail_or_not",
})
df["sample"]  = sample
df["is_tail"] = df["tail_or_not"].astype(int) == 1
n_total = len(df)
n_tail  = int(df["is_tail"].sum())
df.to_csv(tsv_dst, sep="\t", index=False)
print(f"Proteins: {n_total}  |  Tail proteins: {n_tail}")
print(f"TSV written: {tsv_dst}")
PYEOF

conda run -n env_phage_ml \
    python "$CONVERT_SCRIPT" "$PROT_OUT" "$TSV_OUT" "$PREFIX"

rm -f "$CONVERT_SCRIPT"

if [ ! -f "$TSV_OUT" ]; then
    echo "ERROR: TSV conversion failed, file not found: $TSV_OUT" >&2
    exit 2
fi

echo "PhageTailFinder done: $PREFIX"
echo "  Raw per-protein : $PROT_OUT"
echo "  Pipeline TSV    : $TSV_OUT"