import sys
import csv
import matplotlib.patches as mpatches
from pygenomeviz import GenomeViz
from Bio import SeqIO

# ── inputs ────────────────────────────────────────────────────────────────────
if len(sys.argv) != 4:
    print("usage: visualize_genome.py <annotated.gbk> <protein_table.tsv> <out_prefix>")
    sys.exit(1)

gbk_path   = sys.argv[1]
tsv_path   = sys.argv[2]
out_prefix = sys.argv[3]

# ── category → color ──────────────────────────────────────────────────────────
CATEGORY_COLORS = {
    "head and packaging":                                    "#4E79A7",
    "tail":                                                  "#F28E2B",
    "connector":                                             "#E15759",
    "lysis":                                                 "#76B7B2",
    "DNA, RNA and nucleotide metabolism":                    "#59A14F",
    "transcription regulation":                              "#EDC948",
    "integration and excision":                              "#B07AA1",
    "moron, auxiliary metabolic gene and host takeover":     "#FF9DA7",
    "other":                                                 "#9C755F",
    "unknown function":                                      "#BAB0AC",
    "none":                                                  "#D3D3D3",
}
DEFAULT_COLOR = "#D3D3D3"

# ── load protein table ────────────────────────────────────────────────────────
protein_info = {}
with open(tsv_path) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        protein_info[row["protein_id"]] = (row["category"], row["confidence"])

# ── load genome record ────────────────────────────────────────────────────────
records = list(SeqIO.parse(gbk_path, "genbank"))
if not records:
    print(f"ERROR: no records found in {gbk_path}")
    sys.exit(1)

record     = records[0]
genome_len = len(record.seq)
phage_name = record.name if record.name and record.name != "." else out_prefix.split("/")[-1]

print(f"genome: {phage_name}  length: {genome_len} bp  CDS: {sum(1 for f in record.features if f.type == 'CDS')}")

# ── collect CDS features ──────────────────────────────────────────────────────
features = []
for feat in record.features:
    if feat.type != "CDS":
        continue
    start  = int(feat.location.start) + 1
    end    = int(feat.location.end)
    strand = feat.location.strand

    locus_tag = feat.qualifiers.get("locus_tag", [""])[0]
    product   = feat.qualifiers.get("product", ["hypothetical protein"])[0]

    if locus_tag in protein_info:
        category, confidence = protein_info[locus_tag]
    else:
        category, confidence = "none", "NA"

    color = CATEGORY_COLORS.get(category, DEFAULT_COLOR)
    features.append((start, end, strand, locus_tag, product, category, confidence, color))

# ── build figure ──────────────────────────────────────────────────────────────
gv = GenomeViz(
    fig_width=20,
    fig_track_height=0.8,
    track_align_type="left",
    feature_track_ratio=0.4,
)
gv.set_scale_bar()

track = gv.add_feature_track(phage_name, genome_len)

for (start, end, strand, locus_tag, product, category, confidence, color) in features:
    track.add_feature(
        start  = start,
        end    = end,
        strand = strand,
        label  = "",
        fc     = color,
        ec     = "black",
        lw     = 0.3,
    )

# ── legend ────────────────────────────────────────────────────────────────────
present_cats = sorted(set(f[5] for f in features))
legend_handles = []
for cat in present_cats:
    color = CATEGORY_COLORS.get(cat, DEFAULT_COLOR)
    legend_handles.append(mpatches.Patch(facecolor=color, edgecolor="black", linewidth=0.5, label=cat))

fig = gv.plotfig()
fig.legend(
    handles=legend_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, -0.15),
    ncol=5,
    fontsize=8,
    frameon=True,
    title="Functional Category",
    title_fontsize=9,
)
fig.suptitle(
    f"{phage_name} — Functional Genome Map ({genome_len:,} bp)",
    fontsize=13,
    y=1.5
)

# ── save ──────────────────────────────────────────────────────────────────────
png_path = out_prefix + ".png"
svg_path = out_prefix + ".svg"

fig.savefig(
    png_path,
    dpi=200,
    bbox_inches="tight",
    pad_inches=1.0
)

fig.savefig(
    svg_path,
    bbox_inches="tight",
    pad_inches=1.0
)

print(f"saved: {png_path}")
print(f"saved: {svg_path}")