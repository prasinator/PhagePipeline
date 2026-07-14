#!/bin/bash
set -euo pipefail

# =============================================================================
# run_phagerbpdetect_v4.sh
#
# Usage:
#   run_phagerbpdetect_v4.sh <proteins.faa> <prefix> <outdir>
#
# CONCURRENCY SAFETY
# ------------------
# PhageRBPdetect_v4_inference.py hardcodes two shared files:
#   data/sequences.fasta   (input,  overwritten each run)
#   data/predictions.csv   (output, overwritten each run)
#
# If Snakemake runs multiple samples in parallel, these files get clobbered
# producing silent wrong results — no crash, no warning, wrong therapeutics.
#
# Fix: flock on a lockfile in the tool's data/ dir. Each sample acquires the
# lock, stages its input, runs inference, copies output, then releases.
# Other samples block at flock until the lock is free.
# This serializes only the inference step; all other pipeline stages remain
# fully parallel.
# =============================================================================

PROTEINS=${1:?Usage: run_phagerbpdetect_v4.sh <proteins.faa> <prefix> <outdir>}
PROTEINS=$(realpath "$PROTEINS")
PREFIX=${2:?Usage: run_phagerbpdetect_v4.sh <proteins.faa> <prefix> <outdir>}
OUTDIR=${3:?Usage: run_phagerbpdetect_v4.sh <proteins.faa> <prefix> <outdir>}
OUTDIR=$(realpath -m "$OUTDIR")

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
TOOL="$ROOT/tools/PhageRBPdetection"
DATA="$TOOL/data"
MODEL="$DATA/RBPdetect_v4_ESMfine"
INFER="$TOOL/PhageRBPdetect_v4/PhageRBPdetect_v4_inference.py"

SEQ_DST="$DATA/sequences.fasta"
PRED_SRC="$DATA/predictions.csv"
LOCKFILE="$DATA/.phagerbp.lock"

# ── validate ──────────────────────────────────────────────────────────────────
if [ ! -f "$PROTEINS" ]; then
    echo "ERROR: input .faa not found: $PROTEINS" >&2
    exit 2
fi
if [[ "$(head -c 1 "$PROTEINS")" != ">" ]]; then
    echo "ERROR: input does not look like FASTA: $PROTEINS" >&2
    exit 2
fi
if [ ! -d "$MODEL" ]; then
    echo "ERROR: v4 model dir missing: $MODEL" >&2
    exit 2
fi
if [ ! -f "$INFER" ]; then
    echo "ERROR: inference script not found: $INFER" >&2
    exit 2
fi

mkdir -p "$OUTDIR"

# ── acquire lock, run inference, release lock ─────────────────────────────────
# flock -x: exclusive lock (blocks until acquired)
# flock releases automatically when the subshell exits (success or error)
echo "[$PREFIX] Waiting for PhageRBPdetect lock..."
(
    flock -x 200

    echo "[$PREFIX] Lock acquired — staging input"
    cp "$PROTEINS" "$SEQ_DST"

    echo "[$PREFIX] Running PhageRBPdetect v4 inference..."
    cd "$TOOL/PhageRBPdetect_v4"
    conda run -n env_phage_ml \
        python PhageRBPdetect_v4_inference.py

    if [ ! -f "$PRED_SRC" ]; then
        echo "ERROR: [$PREFIX] predictions.csv not written at $PRED_SRC" >&2
        exit 2
    fi

    echo "[$PREFIX] Copying output"
    cp "$PRED_SRC" "$OUTDIR/${PREFIX}_predictions.csv"

    # clean up shared files before releasing lock
    rm -f "$SEQ_DST" "$PRED_SRC"

    echo "[$PREFIX] Lock released"

) 200>"$LOCKFILE"

# ── post-process outside the lock (pure per-sample work, safe to parallelize) ─
# ── post-process outside the lock (pure per-sample work, safe to parallelize) ─
CONVERT_SCRIPT=$(mktemp /tmp/phagerbpdetect_convert_XXXXXX.py)
cat > "$CONVERT_SCRIPT" <<'PYEOF'
import sys
import pandas as pd
csv_in   = sys.argv[1]
tsv_pred = sys.argv[2]
tsv_summ = sys.argv[3]
sample   = sys.argv[4]
df = pd.read_csv(csv_in)
df_out = df.rename(columns={
    "protein_name": "protein_id",
    "preds":        "rbp_prediction",
    "score":        "rbp_score"
})
df_out["rbp_label"] = df_out["rbp_prediction"].map({1: "RBP", 0: "non-RBP"})
df_out.to_csv(tsv_pred, sep="\t", index=False)
n_total = len(df_out)
n_rbp   = int(df_out["rbp_prediction"].sum())
summary = pd.DataFrame([{
    "sample":         sample,
    "total_proteins": n_total,
    "rbp_count":      n_rbp,
    "non_rbp_count":  n_total - n_rbp,
    "rbp_fraction":   round(n_rbp / n_total, 4) if n_total else 0,
}])
summary.to_csv(tsv_summ, sep="\t", index=False)
print(f"Proteins: {n_total}  |  RBPs predicted: {n_rbp}")
print(f"Predictions TSV : {tsv_pred}")
print(f"Summary TSV     : {tsv_summ}")
PYEOF

conda run -n env_phage_ml \
    python "$CONVERT_SCRIPT" \
        "$OUTDIR/${PREFIX}_predictions.csv" \
        "$OUTDIR/${PREFIX}_rbp_predictions.tsv" \
        "$OUTDIR/${PREFIX}_rbp_summary.tsv" \
        "$PREFIX"

rm -f "$CONVERT_SCRIPT"

if [ ! -f "$OUTDIR/${PREFIX}_rbp_predictions.tsv" ] || [ ! -f "$OUTDIR/${PREFIX}_rbp_summary.tsv" ]; then
    echo "ERROR: RBP TSV conversion failed" >&2
    exit 2
fi

echo "PhageRBPdetect v4 done: $PREFIX"
echo "  CSV (raw)  : $OUTDIR/${PREFIX}_predictions.csv"
echo "  TSV (pred) : $OUTDIR/${PREFIX}_rbp_predictions.tsv"
echo "  TSV (summ) : $OUTDIR/${PREFIX}_rbp_summary.tsv"