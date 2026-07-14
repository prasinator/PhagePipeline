#!/bin/bash
set -euo pipefail

# =============================================================================
# run_phagedpo.sh
#
# Usage:
#   run_phagedpo.sh <pharokka_outdir> <prefix> <outdir>
#
# phagedpo_cli.py writes HTML as:
#   os.path.splitext(self.seqfile)[0] + '_output.html'
# where self.seqfile = full path to the .ffn file it finds in INPUT_DIR.
#
# So if it scans INPUT_DIR and finds Lambda_cds.ffn, it writes:
#   INPUT_DIR/Lambda_cds_output.html
#
# Strategy: copy .ffn into OUTDIR so the HTML lands in OUTDIR, not pharokka/.
# =============================================================================

INPUT_DIR=${1:?Usage: run_phagedpo.sh <pharokka_outdir> <prefix> <outdir>}
INPUT_DIR=$(realpath "$INPUT_DIR")
PREFIX=${2:?Usage: run_phagedpo.sh <pharokka_outdir> <prefix> <outdir>}
OUTDIR=${3:?Usage: run_phagedpo.sh <pharokka_outdir> <prefix> <outdir>}
OUTDIR=$(realpath -m "$OUTDIR")

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
TOOL="$ROOT/tools/phagedpo"

CDS_FFN="$INPUT_DIR/${PREFIX}_cds.ffn"
COPIED_FASTA="$OUTDIR/${PREFIX}_cds.fasta"
HTML_OUT="$OUTDIR/${PREFIX}_cds_output.html"
TSV_OUT="$OUTDIR/${PREFIX}_cds_output.tsv"

# ── validate ──────────────────────────────────────────────────────────────────
if [ ! -f "$CDS_FFN" ]; then
    echo "ERROR: input FFN not found: $CDS_FFN" >&2
    exit 2
fi

mkdir -p "$OUTDIR"

# ── copy .ffn into OUTDIR so HTML writes there ────────────────────────────────
# phagedpo_cli.py scans -i <dir> for .fasta files, finds our .ffn,
# then writes HTML next to it: OUTDIR/{PREFIX}_cds_output.html
cp "$CDS_FFN" "$COPIED_FASTA"

# ── run phagedpo ──────────────────────────────────────────────────────────────
echo "[$PREFIX] Running PhageDPO..."
cd "$TOOL"
conda run -n env_phage_ml \
    python phagedpo_cli.py \
    -i "$OUTDIR"

# ── verify HTML was written ───────────────────────────────────────────────────
if [ ! -f "$HTML_OUT" ]; then
    echo "ERROR: PhageDPO did not write expected HTML: $HTML_OUT" >&2
    echo "       Files in $OUTDIR:" >&2
    ls -la "$OUTDIR" >&2
    exit 2
fi

# ── convert HTML → TSV ───────────────────────────────────────────────────────
# ── convert HTML → TSV ───────────────────────────────────────────────────────
CONVERT_SCRIPT=$(mktemp /tmp/phagedpo_html2tsv_XXXXXX.py)
cat > "$CONVERT_SCRIPT" <<'PYEOF'
import sys
import pandas as pd

html_path = sys.argv[1]
tsv_path  = sys.argv[2]

tables = pd.read_html(html_path)
if not tables:
    print(f"ERROR: no tables found in {html_path}", file=sys.stderr)
    sys.exit(1)
df = tables[0]
df.to_csv(tsv_path, sep="\t", index=False)
print(f"TSV written: {tsv_path}  ({len(df)} rows)")
PYEOF

conda run -n env_phage_ml \
    python "$CONVERT_SCRIPT" "$HTML_OUT" "$TSV_OUT"

rm -f "$CONVERT_SCRIPT"

if [ ! -f "$TSV_OUT" ]; then
    echo "ERROR: TSV conversion failed, file not found: $TSV_OUT" >&2
    exit 2
fi

# ── clean up copied .ffn ─────────────────────────────────────────────────────
rm -f "$COPIED_FASTA"

echo "PhageDPO done: $PREFIX"
echo "  HTML: $HTML_OUT"
echo "  TSV : $TSV_OUT"