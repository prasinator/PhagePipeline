import sys
import os
import csv
import glob

# =============================================================================
# summarize_sample.py
# Cross-checking version — BACPHLIP x annotation markers, RGI x CARD,
# VFDB x annotation toxins, AcrDB x annotation keywords, lysis named products
# =============================================================================

if len(sys.argv) != 12:
    print(
        "usage: summarize_sample.py <sample> <checkv_tsv> <bacphlip> "
        "<rgi_txt> <card_tsv> <vfdb_tsv> <acrdb_tsv> <protein_table_tsv> "
        "<pharokka_tsv> <output_tsv> <output_txt>"
    )
    sys.exit(1)

sample        = sys.argv[1]
checkv_path   = sys.argv[2]
bacphlip_path = sys.argv[3]
rgi_path      = sys.argv[4]
card_path     = sys.argv[5]
vfdb_path     = sys.argv[6]
acrdb_path    = sys.argv[7]
protein_path  = sys.argv[8]
pharokka_path = sys.argv[9]
output_tsv    = sys.argv[10]
output_txt    = sys.argv[11]

ACRDB_HIGH   = 1e-5
ACRDB_MEDIUM = 1e-3

INTEGRASE_KEYWORDS   = {"integrase"}
RECOMBINASE_KEYWORDS = {"recombinase"}
TRANSPOSASE_KEYWORDS = {"transposase"}
LYSIS_NAMED_KEYWORDS = {"holin", "endolysin", "lysin", "spanin", "anti-holin"}
TOXIN_KEYWORDS       = {"toxin", "cytotoxin", "enterotoxin", "exotoxin"}
ANTICRISPR_KEYWORDS  = {"anti-crispr", "anticrispr", "acr"}

rows = []

def add(tool, metric, value, detail="", flag=""):
    rows.append({"sample": sample, "tool": tool, "metric": metric,
                 "value": str(value), "detail": detail, "flag": flag})

# ── parse protein table first (used by multiple checks) ──────────────────────
cat_counts       = {}
named_lysis      = []
integrase_hits   = []
recombinase_hits = []
transposase_hits = []
toxin_hits_annot = []
anticrispr_annot = []
total_cds        = 0
hypothetical_cds = 0

if os.path.exists(protein_path):
    with open(protein_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cat     = row.get("category", "unknown function")
            product = row.get("product",  "").lower()
            conf    = row.get("confidence", "")
            total_cds += 1
            cat_counts[cat] = cat_counts.get(cat, 0) + 1
            if "hypothetical" in product or product == "":
                hypothetical_cds += 1
            if any(k in product for k in INTEGRASE_KEYWORDS):
                integrase_hits.append(row.get("product",""))
            if any(k in product for k in RECOMBINASE_KEYWORDS):
                recombinase_hits.append(row.get("product",""))
            if any(k in product for k in TRANSPOSASE_KEYWORDS):
                transposase_hits.append(row.get("product",""))
            if cat == "lysis" and conf == "high":
                if any(k in product for k in LYSIS_NAMED_KEYWORDS):
                    named_lysis.append(row.get("product",""))
            if any(k in product for k in TOXIN_KEYWORDS):
                toxin_hits_annot.append(row.get("product",""))
            if any(k in product for k in ANTICRISPR_KEYWORDS):
                anticrispr_annot.append(row.get("product",""))

# ── pharokka summary TSV ──────────────────────────────────────────────────────
pharokka_data = {}
if os.path.exists(pharokka_path):
    with open(pharokka_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pharokka_data = row
            break

coding_density = pharokka_data.get("cds_coding_density", "")
_gc_raw        = pharokka_data.get("gc_perc", "")
gc_content     = f"{float(_gc_raw)*100:.2f}" if _gc_raw else ""
taxonomy       = pharokka_data.get("phage_description", "")  # not in this TSV, will be empty

# ── CheckV ────────────────────────────────────────────────────────────────────
checkv_data = {}
if os.path.exists(checkv_path):
    with open(checkv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            checkv_data = row
            break
    quality       = checkv_data.get("checkv_quality", "Not-determined")
    completeness  = checkv_data.get("completeness",  "0")
    contamination = checkv_data.get("contamination", "0")
    viral_genes   = checkv_data.get("viral_genes",   "0")
    host_genes    = checkv_data.get("host_genes",    "0")
    genome_length = checkv_data.get("contig_length", "")
    add("CheckV", "genome_length",  genome_length, detail="bp")
    add("CheckV", "checkv_quality", quality,
        flag="WARN" if quality in ("Low-quality","Not-determined") else "")
    add("CheckV", "completeness",   completeness, detail="%",
        flag="WARN" if float(completeness or 0) < 70 else "")
    add("CheckV", "contamination",  contamination, detail="%",
        flag="WARN" if float(contamination or 0) > 5 else "")
    add("CheckV", "viral_genes",    viral_genes)
    add("CheckV", "host_genes",     host_genes,
        flag="WARN" if int(host_genes or 0) > 0 else "")
else:
    quality = "Not-determined"; completeness = contamination = "0"; genome_length = ""
    add("CheckV", "status", "MISSING", flag="ERROR")

# ── Pharokka metrics ──────────────────────────────────────────────────────────
if pharokka_data:
    add("Pharokka", "coding_density", coding_density, detail="%",
        flag="WARN" if coding_density and float(coding_density) < 85 else "")
    add("Pharokka", "gc_content",     gc_content,  detail="%")
    add("Pharokka", "hypothetical_cds", hypothetical_cds)
    add("Pharokka", "taxonomy",       taxonomy)
else:
    add("Pharokka", "status", "MISSING", flag="ERROR")

# ── BACPHLIP ──────────────────────────────────────────────────────────────────
bacphlip_lifestyle = ""
bacphlip_temperate = bacphlip_virulent = 0.0
if os.path.exists(bacphlip_path):
    with open(bacphlip_path) as fh:
        lines = [l.strip() for l in fh if l.strip()]

    if len(lines) >= 2:
        headers = lines[0].split("\t")
        values  = lines[1].split("\t")

        # Some BACPHLIP versions write an index column ("0")
        if len(values) == len(headers) + 1:
            values = values[1:]

        print(f"[BACPHLIP DEBUG] headers={headers}", file=sys.stderr)
        print(f"[BACPHLIP DEBUG] values={values}", file=sys.stderr)

        brow = dict(zip(headers, values))

        bacphlip_virulent = float(brow.get("Virulent", 0))
        bacphlip_temperate = float(brow.get("Temperate", 0))

        bacphlip_lifestyle = (
            "Temperate" if bacphlip_temperate >= 0.6 else
            "Ambiguous" if bacphlip_temperate >= 0.4 else
            "Virulent"
        )

        add("BACPHLIP", "virulent_score", f"{bacphlip_virulent:.4f}")
        add("BACPHLIP", "temperate_score", f"{bacphlip_temperate:.4f}")
        add(
            "BACPHLIP",
            "bacphlip_call",
            bacphlip_lifestyle,
            flag="DISQUALIFY" if bacphlip_lifestyle == "Temperate"
                 else "WARN" if bacphlip_lifestyle == "Ambiguous"
                 else ""
        )
    else:
        add("BACPHLIP", "status", "EMPTY", flag="ERROR")
else:
    add("BACPHLIP", "status", "MISSING", flag="ERROR")


# Lifestyle markers from annotation
integration_cat_n = cat_counts.get("integration and excision", 0)
add("LifestyleMarkers", "integrase_proteins",   len(integrase_hits),
    detail="; ".join(integrase_hits)   if integrase_hits   else "none",
    flag="WARN" if integrase_hits else "")
add("LifestyleMarkers", "recombinase_proteins", len(recombinase_hits),
    detail="; ".join(recombinase_hits) if recombinase_hits else "none",
    flag="WARN" if recombinase_hits else "")
add("LifestyleMarkers", "transposase_proteins", len(transposase_hits),
    detail="; ".join(transposase_hits) if transposase_hits else "none",
    flag="WARN" if transposase_hits else "")
add("LifestyleMarkers", "integration_excision_category_count", integration_cat_n,
    flag="WARN" if integration_cat_n > 0 else "")

annot_markers = integrase_hits + recombinase_hits + transposase_hits
if bacphlip_lifestyle == "Temperate":
    lc = "TEMPERATE"; lconf = "HIGH" if annot_markers else "MEDIUM"; lf = "DISQUALIFY"
elif bacphlip_lifestyle == "Virulent" and not annot_markers:
    lc = "LYTIC";               lconf = "HIGH";   lf = ""
elif bacphlip_lifestyle == "Virulent" and annot_markers:
    lc = "LYTIC_WITH_MARKERS";  lconf = "LOW";    lf = "WARN"
elif bacphlip_lifestyle == "Ambiguous":
    lc = "AMBIGUOUS";           lconf = "LOW";    lf = "WARN"
else:
    lc = "UNKNOWN";             lconf = "NONE";   lf = "WARN"

add("LifestyleConsensus", "consensus",  lc,
    detail=f"BACPHLIP={bacphlip_lifestyle}; markers={len(annot_markers)}", flag=lf)
add("LifestyleConsensus", "confidence", lconf)

# ── RGI + Abricate-CARD (AMR cross-check) ────────────────────────────────────
rgi_strict = []; rgi_loose = []
if os.path.exists(rgi_path):
    with open(rgi_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cutoff = row.get("Cut_Off","").strip()
            entry  = f"{row.get('Best_Hit_ARO','')} [{row.get('Drug Class','')}]"
            if cutoff == "Strict": rgi_strict.append(entry)
            elif cutoff in ("Perfect","Loose"): rgi_loose.append(entry)
    add("RGI","strict_hit_count",len(rgi_strict),
        detail="; ".join(rgi_strict) if rgi_strict else "none",
        flag="DISQUALIFY" if rgi_strict else "")
    add("RGI","loose_hit_count", len(rgi_loose),
        detail="; ".join(rgi_loose)  if rgi_loose  else "none",
        flag="WARN" if rgi_loose else "")
else:
    add("RGI","status","MISSING",flag="ERROR")

card_hits = []
if os.path.exists(card_path):
    with open(card_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.strip().split("\t")
            card_hits.append(f"{p[5] if len(p)>5 else ''} ({p[10] if len(p)>10 else ''}% id)")
    add("Abricate-CARD","hit_count",len(card_hits),
        detail="; ".join(card_hits) if card_hits else "none",
        flag="DISQUALIFY" if card_hits else "")
else:
    add("Abricate-CARD","status","MISSING",flag="ERROR")

amr_disq = bool(rgi_strict or card_hits)
amr_note = ("CONFIRMED by RGI+CARD" if rgi_strict and card_hits else
            "flagged by RGI only"   if rgi_strict else
            "flagged by CARD only"  if card_hits  else "clean")
add("AMRConsensus","amr_status","DISQUALIFY" if amr_disq else "CLEAN",
    detail=amr_note, flag="DISQUALIFY" if amr_disq else "")

# ── Abricate-VFDB + annotation toxins (virulence cross-check) ────────────────
vfdb_hits = []
if os.path.exists(vfdb_path):
    with open(vfdb_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip(): continue
            p = line.strip().split("\t")
            vfdb_hits.append(f"{p[5] if len(p)>5 else ''} ({p[10] if len(p)>10 else ''}% id)")
    add("Abricate-VFDB","hit_count",len(vfdb_hits),
        detail="; ".join(vfdb_hits) if vfdb_hits else "none",
        flag="DISQUALIFY" if vfdb_hits else "")
else:
    add("Abricate-VFDB","status","MISSING",flag="ERROR")

add("VirulenceMarkers","toxin_keywords_in_annotation",len(toxin_hits_annot),
    detail="; ".join(toxin_hits_annot) if toxin_hits_annot else "none",
    flag="WARN" if toxin_hits_annot else "")

# VFDB hit = hard evidence → DISQUALIFY
# Annotation keyword only = soft evidence → MANUAL_REVIEW
if vfdb_hits and toxin_hits_annot:
    vf_status = "DISQUALIFY"
    vf_note   = "CONFIRMED by VFDB+annotation"
elif vfdb_hits:
    vf_status = "DISQUALIFY"
    vf_note   = "flagged by VFDB only"
elif toxin_hits_annot:
    vf_status = "MANUAL_REVIEW"
    vf_note   = "annotation keyword only — verify manually"
else:
    vf_status = "CLEAN"
    vf_note   = "clean"

vf_flag = "DISQUALIFY" if vf_status == "DISQUALIFY" else ("WARN" if vf_status == "MANUAL_REVIEW" else "")
add("VirulenceConsensus","virulence_status", vf_status, detail=vf_note, flag=vf_flag)

# ── AcrDB + annotation anti-CRISPR cross-check ───────────────────────────────
acrdb_high_hits = acrdb_medium_hits = acrdb_weak_hits = []
if os.path.exists(acrdb_path):
    acrdb_high_hits = []; acrdb_medium_hits = []; acrdb_weak_hits = []
    with open(acrdb_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"): continue
            p = line.split("\t")
            if len(p) < 11: continue
            try: ev = float(p[10])
            except ValueError: continue
            entry = f"{p[1]} (e={ev:.2g}, {p[2]}% id)"
            if ev < ACRDB_HIGH:       acrdb_high_hits.append(entry)
            elif ev < ACRDB_MEDIUM:   acrdb_medium_hits.append(entry)
            else:                     acrdb_weak_hits.append(entry)
    add("AcrDB","high_confidence_hits",  len(acrdb_high_hits),
        detail="; ".join(acrdb_high_hits[:5])   if acrdb_high_hits   else "none",
        flag="WARN" if acrdb_high_hits else "")
    add("AcrDB","medium_confidence_hits",len(acrdb_medium_hits),
        detail="; ".join(acrdb_medium_hits[:5]) if acrdb_medium_hits else "none",
        flag="WARN" if acrdb_medium_hits else "")
    add("AcrDB","weak_hits",             len(acrdb_weak_hits),
        detail=f"count={len(acrdb_weak_hits)}")
else:
    add("AcrDB","status","MISSING",flag="ERROR")
    acrdb_high_hits = acrdb_medium_hits = acrdb_weak_hits = []

add("AcrDBAnnotation","anticrispr_keywords_in_annotation",len(anticrispr_annot),
    detail="; ".join(anticrispr_annot) if anticrispr_annot else "none",
    flag="WARN" if anticrispr_annot else "")

acr_pos  = bool(acrdb_high_hits or acrdb_medium_hits)
acr_both = acr_pos and bool(anticrispr_annot)
acr_note = ("CONFIRMED by AcrDB+annotation" if acr_both else
            "flagged by AcrDB only"          if acr_pos         else
            "flagged by annotation only"     if anticrispr_annot else
            "no anti-CRISPR detected")
acr_conf = "HIGH" if acr_both else ("MEDIUM" if acr_pos or anticrispr_annot else "NONE")
add("AcrConsensus","acr_status",     acr_note, flag="WARN" if acr_pos or anticrispr_annot else "")
add("AcrConsensus","acr_confidence", acr_conf)

# ── Lysis cross-check ─────────────────────────────────────────────────────────
lysis_cat_n = cat_counts.get("lysis", 0)
add("Lysis","lysis_category_count",lysis_cat_n,
    flag="WARN" if lysis_cat_n == 0 else "")
add("Lysis","named_lysis_proteins",len(named_lysis),
    detail="; ".join(named_lysis) if named_lysis else "none",
    flag="WARN" if not named_lysis else "")

lysis_note = ("CONFIRMED — named lysis proteins detected"          if named_lysis else
              "INFERRED — lysis category only, no named proteins"  if lysis_cat_n > 0 else
              "NOT DETECTED")
add("LysisConsensus","lysis_status",lysis_note,
    flag="" if named_lysis else "WARN")

# ── functional annotation summary ────────────────────────────────────────────
add("FunctionalAnnotation","total_cds",total_cds)
add("FunctionalAnnotation","hypothetical_cds",hypothetical_cds,
    detail=f"{int(hypothetical_cds/total_cds*100) if total_cds else 0}% hypothetical")
for cat in sorted(cat_counts):
    safe = cat.replace(" ","_").replace(",","").replace("/","_")
    add("FunctionalAnnotation",f"cat_{safe}",cat_counts[cat],detail=cat)
add("FunctionalAnnotation","has_tail_proteins",
    "yes" if cat_counts.get("tail",0)>0 else "no",
    flag="" if cat_counts.get("tail",0)>0 else "WARN")
add("FunctionalAnnotation","has_lysis_proteins",
    "yes" if lysis_cat_n>0 else "no",
    flag="" if lysis_cat_n>0 else "WARN")

# =============================================================================
# Write TSV
# =============================================================================
os.makedirs(os.path.dirname(output_tsv) if os.path.dirname(output_tsv) else ".", exist_ok=True)
with open(output_tsv,"w",newline="") as fh:
    writer = csv.DictWriter(fh,fieldnames=["sample","tool","metric","value","detail","flag"],delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

# =============================================================================
# Write human-readable .txt report
# =============================================================================
sep = "=" * 70
lines = [
    sep,
    f"  PHAGE CANDIDATE REPORT — {sample}",
    sep, "",
    "── GENOME ──────────────────────────────────────────────────────────",
    f"  Sample name           : {sample}",
    f"  Genome length         : {genome_length} bp",
    f"  GC content            : {gc_content} %",
    f"  Coding density        : {coding_density} %"
        + (" ⚠ LOW — check gene calling" if coding_density and float(coding_density)<85 else ""),
    f"  CheckV quality        : {quality}",
    f"  Completeness          : {completeness} %",
    f"  Contamination         : {contamination} %",
    f"  Taxonomy (Pharokka)   : {taxonomy if taxonomy else 'N/A'}",
    "",
    "── ANNOTATION ──────────────────────────────────────────────────────",
    f"  Total CDS             : {total_cds}",
    f"  Hypothetical CDS      : {hypothetical_cds} ({int(hypothetical_cds/total_cds*100) if total_cds else 0}%)",
    "  Functional categories :",
]
for cat in sorted(cat_counts, key=lambda c: -cat_counts[c]):
    lines.append(f"    {cat_counts[cat]:>4}  {cat}")

lines += [
    "",
    "── LIFESTYLE ───────────────────────────────────────────────────────",
    f"  BACPHLIP              : {bacphlip_lifestyle} (virulent={bacphlip_virulent:.3f}, temperate={bacphlip_temperate:.3f})",
    f"  Integrase proteins    : {len(integrase_hits)}"
        + (f"  — {', '.join(integrase_hits)}" if integrase_hits else ""),
    f"  Recombinase proteins  : {len(recombinase_hits)}"
        + (f"  — {', '.join(recombinase_hits)}" if recombinase_hits else ""),
    f"  Transposase proteins  : {len(transposase_hits)}"
        + (f"  — {', '.join(transposase_hits)}" if transposase_hits else ""),
    f"  >>> VERDICT           : {lc} (confidence: {lconf})",
    "",
    "── SAFETY: AMR ─────────────────────────────────────────────────────",
    f"  RGI Strict hits       : {len(rgi_strict)}"
        + (f"  — {'; '.join(rgi_strict)}" if rgi_strict else ""),
    f"  Abricate-CARD hits    : {len(card_hits)}"
        + (f"  — {'; '.join(card_hits)}" if card_hits else ""),
    f"  >>> VERDICT           : {amr_note.upper()}",
    "",
    "── SAFETY: VIRULENCE ───────────────────────────────────────────────",
    f"  Abricate-VFDB hits    : {len(vfdb_hits)}"
        + (f"  — {'; '.join(vfdb_hits)}" if vfdb_hits else ""),
    f"  Annotation toxin kw   : {len(toxin_hits_annot)}"
        + (f"  — {'; '.join(toxin_hits_annot)}" if toxin_hits_annot else ""),
    f"  >>> VERDICT           : {vf_note.upper()}",
    "",
    "── ANTI-CRISPR ─────────────────────────────────────────────────────",
    f"  AcrDB HIGH hits       : {len(acrdb_high_hits)}"
        + (f"  — {'; '.join(acrdb_high_hits[:3])}" if acrdb_high_hits else ""),
    f"  AcrDB MEDIUM hits     : {len(acrdb_medium_hits)}"
        + (f"  — {'; '.join(acrdb_medium_hits[:3])}" if acrdb_medium_hits else ""),
    f"  AcrDB WEAK hits       : {len(acrdb_weak_hits)}",
    f"  Annotation acr kw     : {len(anticrispr_annot)}"
        + (f"  — {'; '.join(anticrispr_annot)}" if anticrispr_annot else ""),
    f"  >>> VERDICT           : {acr_note} (confidence: {acr_conf})",
    "",
    "── LYSIS MACHINERY ─────────────────────────────────────────────────",
    f"  Lysis category CDS    : {lysis_cat_n}",
    f"  Named lysis proteins  : {len(named_lysis)}"
        + (f"  — {', '.join(named_lysis)}" if named_lysis else ""),
    f"  >>> VERDICT           : {lysis_note}",
    "",
    sep, "",
]

with open(output_txt,"w") as fh:
    fh.write("\n".join(lines))

flags_list = [r for r in rows if r["flag"] in ("DISQUALIFY","ERROR","WARN")]
print(f"{sample}: {len(rows)} metrics → {output_tsv}")
print(f"         report    → {output_txt}")
disq = [r for r in rows if r["flag"] == "DISQUALIFY"]
if disq:
    for r in disq:
        print(f"  [DISQUALIFY] {r['tool']} / {r['metric']} = {r['value']} | {r['detail']}")
elif flags_list:
    print(f"  {len(flags_list)} warnings")
else:
    print("  no flags")