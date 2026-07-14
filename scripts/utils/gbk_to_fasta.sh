#!/bin/bash
set -euo pipefail

GBK=${1:?Usage: scripts/utils/gbk_to_fasta.sh <input.gbk> <output.fasta> [header]}
OUT=${2:?Usage: scripts/utils/gbk_to_fasta.sh <input.gbk> <output.fasta> [header]}
HEADER=${3:-$(basename "$GBK" .gbk)}

awk -v header="$HEADER" '
  BEGIN { print ">" header }
  /^ORIGIN/ { flag=1; next }
  /^\/\// { flag=0 }
  flag {
    gsub(/[0-9 ]/, "")
    print toupper($0)
  }
' "$GBK" > "$OUT"
