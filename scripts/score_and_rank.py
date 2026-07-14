import sys
import os
import csv

# =============================================================================
# score_and_rank.py — reads consensus fields from tool_summary.tsv
#
# Scoring (100 pts):
#   Quality   20  — CheckV
#   Lifestyle 30  — LifestyleConsensus (BACPHLIP x annotation markers)
#   Safety    30  — AMRConsensus + VirulenceConsensus + AcrDB tiered
#   Lysis     10  — LysisConsensus (named proteins x category)
#   Host      10  — tail / head / connector proteins
#
# Hard disqualifiers → REJECTED:
#   LifestyleConsensus = TEMPERATE
#   AMRConsensus = DISQUALIFY
#   VirulenceConsensus = DISQUALIFY
# =============================================================================

if len(sys.argv) != 4:
    print("usage: score_and_rank.py <reports_run_dir> <samples_csv> <output_tsv>")
    sys.exit(1)

reports_run_dir = sys.argv[1]
samples         = [s.strip() for s in sys.argv[2].split(",") if s.strip()]
output_tsv      = sys.argv[3]

SCORE_PASS   = 80
SCORE_REVIEW = 50


def load_summary(sample, reports_run_dir):
    path = os.path.join(reports_run_dir, sample, f"{sample}_tool_summary.tsv")
    data = {}
    if not os.path.exists(path):
        return None, path
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            data[(row["tool"], row["metric"])] = (row["value"], row["detail"], row["flag"])
    return data, path


def get(data, tool, metric, default=""):
    e = data.get((tool, metric))
    return e[0] if e else default


def score_quality(data):
    quality = get(data, "CheckV", "checkv_quality", "Not-determined")
    pts = {"Complete":20,"High-quality":20,"Medium-quality":14,"Low-quality":7,"Not-determined":0}.get(quality, 0)
    notes = [f"CheckV: {quality}"]
    try:
        cont = float(get(data,"CheckV","contamination","0") or 0)
        comp = float(get(data,"CheckV","completeness","0") or 0)
        if cont > 5:
            pts = max(0, pts-5)
            notes.append(f"contamination {cont:.1f}% -5")
        notes.append(f"completeness {comp:.1f}%")
    except ValueError:
        pass
    return pts, "; ".join(notes)


def score_lifestyle(data):
    consensus  = get(data, "LifestyleConsensus", "consensus",  "UNKNOWN")
    confidence = get(data, "LifestyleConsensus", "confidence", "NONE")
    detail     = (data.get(("LifestyleConsensus","consensus")) or ("","",""))[1]

    if consensus == "TEMPERATE":
        return 0, True, f"TEMPERATE confirmed (confidence: {confidence}) — {detail}"
    if consensus == "LYTIC" and confidence == "HIGH":
        return 30, False, f"LYTIC confirmed by BACPHLIP + no annotation markers"
    if consensus == "LYTIC_WITH_MARKERS":
        return 15, False, f"LYTIC (BACPHLIP) but annotation markers found — {detail}"
    if consensus == "AMBIGUOUS":
        return 15, False, f"AMBIGUOUS BACPHLIP — {detail}"
    return 0, False, f"UNKNOWN lifestyle — {detail}"


def score_safety(data):
    pts = 30
    disq = False
    reject = []
    notes  = []

    # AMR consensus
    amr_status = get(data,"AMRConsensus","amr_status","CLEAN")
    amr_detail = (data.get(("AMRConsensus","amr_status")) or ("","",""))[1]
    if amr_status == "DISQUALIFY":
        disq = True
        reject.append(f"AMR: {amr_detail}")
    else:
        notes.append(f"AMR: {amr_detail}")

    # Virulence consensus
    vf_status = get(data,"VirulenceConsensus","virulence_status","CLEAN")
    vf_detail = (data.get(("VirulenceConsensus","virulence_status")) or ("","",""))[1]
    if vf_status == "DISQUALIFY":
        disq = True
        reject.append(f"Virulence: {vf_detail}")
    elif vf_status == "MANUAL_REVIEW":
        pts = max(0, pts - 10)   # penalty but not rejection
        notes.append(f"Virulence: MANUAL_REVIEW — {vf_detail}")
    else:
        notes.append(f"Virulence: {vf_detail}")

    # AcrDB tiered penalty (not a disqualifier)
    high_n   = int(get(data,"AcrDB","high_confidence_hits","0")   or 0)
    medium_n = int(get(data,"AcrDB","medium_confidence_hits","0") or 0)
    weak_n   = int(get(data,"AcrDB","weak_hits","0")              or 0)
    acr_conf = get(data,"AcrConsensus","acr_confidence","NONE")
    acr_note = get(data,"AcrConsensus","acr_status","")
    if high_n:
        pts = max(0, pts-15)
        notes.append(f"AcrDB HIGH {high_n} hits ({acr_conf}) -15")
    if medium_n:
        pts = max(0, pts-7)
        notes.append(f"AcrDB MEDIUM {medium_n} hits -7")
    if not high_n and not medium_n:
        notes.append(f"AcrDB: {acr_note}")

    if disq:
        pts = 0

    return pts, disq, "; ".join(reject + notes), reject


def score_lysis(data):
    lysis_note = get(data,"LysisConsensus","lysis_status","NOT DETECTED")
    named_n    = int(get(data,"Lysis","named_lysis_proteins","0") or 0)
    cat_n      = int(get(data,"Lysis","lysis_category_count","0") or 0)
    if named_n > 0:
        return 10, f"CONFIRMED — {named_n} named lysis proteins"
    if cat_n > 0:
        return 5,  f"INFERRED — {cat_n} lysis category CDS, no named proteins"
    return 0, "NOT DETECTED"


def score_host(data):
    pts = 0; notes = []
    tail_n = int(get(data,"FunctionalAnnotation","cat_tail","0") or 0)
    head_n = int(get(data,"FunctionalAnnotation","cat_head_and_packaging","0") or 0)
    conn_n = int(get(data,"FunctionalAnnotation","cat_connector","0") or 0)
    if tail_n: pts += 5; notes.append(f"tail: {tail_n}")
    if head_n: pts += 3; notes.append(f"head/packaging: {head_n}")
    if conn_n: pts += 2; notes.append(f"connector: {conn_n}")
    return min(pts,10), "; ".join(notes) if notes else "no structural proteins"


result_rows = []

for sample in samples:
    data, path = load_summary(sample, reports_run_dir)
    if data is None:
        print(f"WARNING: summary not found for {sample} at {path} — skipping")
        continue
    print(f"scoring {sample} ...")

    quality_pts,  quality_note                          = score_quality(data)
    lifestyle_pts, lifestyle_disq, lifestyle_note       = score_lifestyle(data)
    safety_pts, safety_disq, safety_note, reject_reasons = score_safety(data)
    lysis_pts,    lysis_note                            = score_lysis(data)
    host_pts,     host_note                             = score_host(data)

    disq = lifestyle_disq or safety_disq
    if disq:
        total = 0; status = "REJECTED"
    else:
        total  = quality_pts + lifestyle_pts + safety_pts + lysis_pts + host_pts
        status = "PASS" if total >= SCORE_PASS else ("REVIEW" if total >= SCORE_REVIEW else "FAIL")

    rejection_str = "; ".join(reject_reasons)
    if lifestyle_disq:
        rejection_str = ("LIFESTYLE_TEMPERATE; " + rejection_str).rstrip("; ")

    result_rows.append({
        "sample":            sample,
        "status":            status,
        "total_score":       total,
        "quality_score":     quality_pts,
        "lifestyle_score":   lifestyle_pts,
        "safety_score":      safety_pts,
        "lysis_score":       lysis_pts,
        "host_score":        host_pts,
        "rejection_reasons": rejection_str,
        "quality_note":      quality_note,
        "lifestyle_note":    lifestyle_note,
        "safety_note":       safety_note,
        "lysis_note":        lysis_note,
        "host_note":         host_note,
    })

status_order = {"PASS":0,"REVIEW":1,"FAIL":2,"REJECTED":3}
result_rows.sort(key=lambda r: (status_order.get(r["status"],9), -r["total_score"]))

fieldnames = [
    "sample","status","total_score",
    "quality_score","lifestyle_score","safety_score","lysis_score","host_score",
    "rejection_reasons",
    "quality_note","lifestyle_note","safety_note","lysis_note","host_note",
]
os.makedirs(os.path.dirname(output_tsv) if os.path.dirname(output_tsv) else ".", exist_ok=True)
with open(output_tsv,"w",newline="") as fh:
    writer = csv.DictWriter(fh,fieldnames=fieldnames,delimiter="\t")
    writer.writeheader()
    writer.writerows(result_rows)

print(f"\nscored {len(result_rows)} sample(s)  →  {output_tsv}\n")
hdr = f"{'sample':<22} {'status':<10} {'total':>6} {'qual':>5} {'life':>5} {'safe':>5} {'lys':>5} {'host':>5}  rejection"
print(hdr)
print("-"*len(hdr))
for r in result_rows:
    print(f"{r['sample']:<22} {r['status']:<10} {r['total_score']:>6} "
          f"{r['quality_score']:>5} {r['lifestyle_score']:>5} {r['safety_score']:>5} "
          f"{r['lysis_score']:>5} {r['host_score']:>5}  {r['rejection_reasons']}")