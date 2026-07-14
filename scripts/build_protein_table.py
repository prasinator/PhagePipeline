#!/usr/bin/env python3
"""
build_protein_table.py
Produces a flat, biologically focused protein table from the annotated GBK.

One row per CDS. Columns reflect the single best annotation available,
with evidence_source indicating which tool provided it.

Priority logic per CDS:
  1. PHOLD hit (function != "unknown function")
     → evidence_source = PHOLD
  2. PHOLD + PHYNTENY agree on category
     → evidence_source = PHOLD+PHYNTENY
  3. PHOLD unknown, phynteny has prediction
     → evidence_source = PHYNTENY
  4. Nothing useful
     → product = hypothetical protein, category = unknown function,
       confidence = none, evidence_source = none

Usage:
    python build_protein_table.py <sample_annotated.gbk> <sample> <outdir>

Output:
    <outdir>/<sample>_protein_table.tsv
"""

import sys
import os
from collections import Counter
from Bio import SeqIO

PHOLD_UNKNOWN = "unknown function"


def get_qualifier(feat, key, default=""):
    vals = feat.qualifiers.get(key, [default])
    return vals[0].strip() if vals else default


def phold_confidence(feat):
    conf = get_qualifier(feat, "annotation_confidence", "")
    if conf in ("high", "medium", "low"):
        return conf
    if conf == "pharokka":
        return "medium"
    return "none"


def parse_cds(feat):
    protein_id = get_qualifier(feat, "locus_tag") or get_qualifier(feat, "ID")
    start      = int(feat.location.start) + 1
    end        = int(feat.location.end)
    strand     = "+" if feat.location.strand == 1 else "-"
    product    = get_qualifier(feat, "product", "hypothetical protein")

    phold_func   = get_qualifier(feat, "function", PHOLD_UNKNOWN)
    phynteny_cat = get_qualifier(feat, "phynteny_category", "NA")
    phynteny_sup = get_qualifier(feat, "phynteny_support", "none")

    phold_known   = phold_func != PHOLD_UNKNOWN
    phynteny_known = phynteny_cat not in ("NA", "no_match", "")

    if phold_known and phynteny_known:
        if phold_func == phynteny_cat:
            category        = phold_func
            confidence      = phold_confidence(feat)
            evidence_source = "PHOLD+PHYNTENY"
        else:
            category        = phold_func
            confidence      = phold_confidence(feat)
            evidence_source = "PHOLD"

    elif phold_known:
        category        = phold_func
        confidence      = phold_confidence(feat)
        evidence_source = "PHOLD"

    elif phynteny_known:
        category        = phynteny_cat
        confidence      = phynteny_sup
        evidence_source = "PHYNTENY"

    else:
        category        = PHOLD_UNKNOWN
        confidence      = "none"
        evidence_source = "none"

    return {
        "protein_id":      protein_id,
        "start":           start,
        "end":             end,
        "strand":          strand,
        "product":         product,
        "category":        category,
        "confidence":      confidence,
        "evidence_source": evidence_source,
    }


def main():
    if len(sys.argv) < 4:
        print("Usage: build_protein_table.py <annotated.gbk> <sample> <outdir>")
        sys.exit(1)

    gbk_path = sys.argv[1]
    sample   = sys.argv[2]
    outdir   = sys.argv[3]

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{sample}_protein_table.tsv")

    records = list(SeqIO.parse(gbk_path, "genbank"))

    columns = ["protein_id", "start", "end", "strand",
               "product", "category", "confidence", "evidence_source"]

    rows = []
    for record in records:
        for feat in record.features:
            if feat.type != "CDS":
                continue
            rows.append(parse_cds(feat))

    with open(outpath, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in columns) + "\n")

    sources = Counter(r["evidence_source"] for r in rows)
    cats    = Counter(r["category"] for r in rows)

    print(f"Written {len(rows)} proteins to {outpath}")
    print("\nEvidence source breakdown:")
    for src, n in sorted(sources.items(), key=lambda x: -x[1]):
        print(f"  {src:20s} {n}")
    print("\nCategory breakdown:")
    for cat, n in sorted(cats.items(), key=lambda x: -x[1]):
        print(f"  {cat:45s} {n}")


if __name__ == "__main__":
    main()
    