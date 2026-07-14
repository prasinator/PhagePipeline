#!/usr/bin/env python3
"""
build_annotated_gb.py
Enriches phold GBK with phold confidence + phynteny synteny-based predictions.

Two enrichments applied:
  1. annotation_confidence from phold_per_cds_predictions.tsv (joined by cds_id/locus_tag)
  2. phynteny_category, phynteny_confidence, phynteny_support from phynteny.tsv
     (joined by coordinates)

Usage:
    python build_annotated_gb.py <phold.gbk> <phold_per_cds.tsv> <phynteny.tsv> <sample> <outdir>

Output:
    <outdir>/<sample>_annotated.gbk
"""

import sys
import os
import csv
from Bio import SeqIO

CONFIDENCE_HIGH   = 0.6
CONFIDENCE_MEDIUM = 0.4


def support_tier(conf):
    if conf >= CONFIDENCE_HIGH:
        return "high"
    elif conf >= CONFIDENCE_MEDIUM:
        return "medium"
    else:
        return "low"


def load_phold_by_cds_id(phold_tsv):
    # Keyed by cds_id (matches locus_tag in GBK)
    cds_map = {}
    with open(phold_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cds_id = row.get("cds_id", "").strip()
            if cds_id:
                cds_map[cds_id] = row
    return cds_map


def load_phynteny_by_coords(phynteny_tsv):
    # Keyed by (start, end, strand)
    coord_map = {}
    with open(phynteny_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                start  = int(row["start"])
                end    = int(row["end"])
                strand = row["strand"].strip()
            except (KeyError, ValueError):
                continue
            coord_map[(start, end, strand)] = row
    return coord_map


def strand_str(location):
    if location.strand == 1:
        return "+"
    elif location.strand == -1:
        return "-"
    return "?"


def enrich_record(record, phold_map, phynteny_map):
    phold_matched    = 0
    phynteny_matched = 0
    phynteny_unmatched = 0

    for feat in record.features:
        if feat.type != "CDS":
            continue

        locus_tag = feat.qualifiers.get("locus_tag", [None])[0]

        # --- Phold confidence (join by locus_tag) ---
        phold_row = phold_map.get(locus_tag)
        if phold_row:
            phold_matched += 1
            conf = phold_row.get("annotation_confidence", "").strip()
            feat.qualifiers["annotation_confidence"] = [conf if conf else "none"]
        else:
            feat.qualifiers["annotation_confidence"] = ["none"]

        # --- Phynteny (join by coordinates) ---
        start  = int(feat.location.start)
        end    = int(feat.location.end)
        strand = strand_str(feat.location)

        key = (start, end, strand)
        row = phynteny_map.get(key)

        if row is None:
            for ds in (-1, 1):
                for de in (-1, 1):
                    candidate = (start + ds, end + de, strand)
                    if candidate in phynteny_map:
                        row = phynteny_map[candidate]
                        break
                if row:
                    break

        if row is None:
            phynteny_unmatched += 1
            feat.qualifiers["phynteny_category"]   = ["no_match"]
            feat.qualifiers["phynteny_confidence"] = ["NA"]
            feat.qualifiers["phynteny_support"]    = ["none"]
            continue

        phynteny_matched += 1
        phynteny_cat = row.get("phynteny_category", "NA").strip()
        conf_str     = row.get("phynteny_confidence", "NA").strip()
        score_str    = row.get("phynteny_score", "NA").strip()

        if phynteny_cat in ("NA", ""):
            feat.qualifiers["phynteny_category"]   = ["NA"]
            feat.qualifiers["phynteny_confidence"] = ["NA"]
            feat.qualifiers["phynteny_support"]    = ["none"]
            continue

        try:
            tier = support_tier(float(conf_str))
        except ValueError:
            tier = "unknown"

        feat.qualifiers["phynteny_category"]   = [phynteny_cat]
        feat.qualifiers["phynteny_confidence"] = [conf_str]
        feat.qualifiers["phynteny_score"]      = [score_str]
        feat.qualifiers["phynteny_support"]    = [tier]

    print(f"  Phold join:    {phold_matched} matched")
    print(f"  Phynteny join: {phynteny_matched} matched, {phynteny_unmatched} unmatched")
    return record


def main():
    if len(sys.argv) < 6:
        print("Usage: build_annotated_gb.py <phold.gbk> <phold_per_cds.tsv> <phynteny.tsv> <sample> <outdir>")
        sys.exit(1)

    phold_gbk    = sys.argv[1]
    phold_tsv    = sys.argv[2]
    phynteny_tsv = sys.argv[3]
    sample       = sys.argv[4]
    outdir       = sys.argv[5]

    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"{sample}_annotated.gbk")

    phold_map    = load_phold_by_cds_id(phold_tsv)
    phynteny_map = load_phynteny_by_coords(phynteny_tsv)
    print(f"Loaded {len(phold_map)} phold entries (cds_id-keyed)")
    print(f"Loaded {len(phynteny_map)} phynteny entries (coord-keyed)")

    records = list(SeqIO.parse(phold_gbk, "genbank"))
    print(f"Loaded {len(records)} records from phold GBK")

    enriched = []
    for record in records:
        enriched.append(enrich_record(record, phold_map, phynteny_map))

    SeqIO.write(enriched, outpath, "genbank")
    print(f"Written enriched GBK: {outpath}")


if __name__ == "__main__":
    main()
    