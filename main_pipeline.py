"""
main_pipeline.py
================
chr20 Variant Annotation Pipeline
Author  : Farhana Sayed
Data    : 1000 Genomes Project chr20 subset
Tool    : Ensembl VEP (web) GRCh38 + Python

Run:
    pip3 install pandas matplotlib --break-system-packages
    python3 scripts/main_pipeline.py

Outputs:
    results/chr20_variants.csv
    results/chr20_plots.png
    results/chr20_report.html
"""

import os
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from datetime import datetime

# ── Paths ─────────────────────────────────────────────────
BASE    = os.path.expanduser("~/variant-annotation-pipeline")
VEP_IN  = os.path.join(BASE, "results", "chr20_annotated.txt")
OUT_CSV = os.path.join(BASE, "results", "chr20_variants.csv")
OUT_PNG = os.path.join(BASE, "results", "chr20_plots.png")
OUT_HTM = os.path.join(BASE, "results", "chr20_report.html")

os.makedirs(os.path.join(BASE, "results"), exist_ok=True)


# ══════════════════════════════════════════════════════════
# 1. PARSE VEP FILE
# ══════════════════════════════════════════════════════════
def parse_vep(path):
    print("\n[1/4] Parsing VEP output...")

    if not os.path.exists(path):
        print(f"  ERROR: File not found: {path}")
        print("  Copy chr20_annotated.txt into results/ folder")
        sys.exit(1)

    rows, header = [], []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("## "):
                continue
            elif line.startswith("#"):
                # Header row — strip leading #
                header = line.lstrip("#").split("\t")
            elif line.strip():
                parts = line.split("\t")
                # Pad if row is shorter than header
                while len(parts) < len(header):
                    parts.append("")
                rows.append(parts[:len(header)])

    if not header:
        print("  ERROR: No header found in VEP file")
        sys.exit(1)

    df = pd.DataFrame(rows, columns=header)

    # Replace VEP missing value markers with empty string
    df = df.replace("-", "")

    print(f"  OK: {len(df)} annotation rows loaded")
    print(f"  OK: {len(df.columns)} columns: {list(df.columns[:8])} ...")
    return df


# ══════════════════════════════════════════════════════════
# 2. CLEAN & PROCESS
# ══════════════════════════════════════════════════════════
def process(df):
    print("\n[2/4] Cleaning and processing...")

    # ── Rename VEP columns to friendly names ──────────────
    rename = {
        "Uploaded_variation": "variant_id",
        "Location":           "location",
        "Allele":             "alt_allele",
        "Consequence":        "consequence",
        "IMPACT":             "impact",
        "SYMBOL":             "gene_symbol",
        "Gene":               "ensembl_gene",
        "Feature_type":       "feature_type",
        "Feature":            "transcript_id",
        "BIOTYPE":            "biotype",
        "EXON":               "exon",
        "INTRON":             "intron",
        "HGVSc":              "hgvsc",
        "HGVSp":              "hgvsp",
        "Existing_variation": "rsid",
        "REF_ALLELE":         "ref_allele",
        "STRAND":             "strand",
        "MANE_SELECT":        "mane_select",
        "SIFT":               "sift",
        "PolyPhen":           "polyphen",
        "AF":                 "af_1kg",
        "CLIN_SIG":           "clinvar_sig",
        "APPRIS":             "appris",
        "TSL":                "tsl",
        "DISTANCE":           "distance",
    }
    df = df.rename(columns={k: v for k, v in rename.items()
                             if k in df.columns})

    # ── Keep MANE_Select transcripts preferentially ────────
    # MANE_Select is the gold-standard single transcript per gene
    if "mane_select" in df.columns:
        mane = df[df["mane_select"].str.strip() != ""]
        non_mane = df[df["mane_select"].str.strip() == ""]
        # Use MANE where available, else keep all
        df = pd.concat([mane, non_mane[
            ~non_mane["variant_id"].isin(mane["variant_id"])
        ]], ignore_index=True)
        print(f"  OK: {len(mane)} MANE Select transcript rows identified")

    # ── Extract SIFT prediction label ─────────────────────
    # VEP format: "deleterious(0.02)" → extract "deleterious"
    if "sift" in df.columns:
        df["sift_pred"] = (df["sift"]
                           .str.extract(r"^([a-z_]+)", expand=False)
                           .fillna(""))
        df["sift_score"] = (df["sift"]
                            .str.extract(r"\(([0-9.]+)\)", expand=False)
                            .fillna(""))

    # ── Extract PolyPhen prediction label ─────────────────
    if "polyphen" in df.columns:
        df["polyphen_pred"] = (df["polyphen"]
                               .str.extract(r"^([a-z_]+)", expand=False)
                               .fillna(""))
        df["polyphen_score"] = (df["polyphen"]
                                .str.extract(r"\(([0-9.]+)\)", expand=False)
                                .fillna(""))

    # ── Classify impact from consequence if IMPACT missing ─
    if "impact" not in df.columns or df["impact"].str.strip().eq("").all():
        df["impact"] = df["consequence"].apply(classify_impact)
        print("  OK: Impact classified from consequence column")

    # ── Add human-readable consequence category ────────────
    df["consequence_category"] = df["consequence"].apply(categorise)

    print(f"  OK: Final dataframe: {len(df)} rows x {len(df.columns)} cols")
    return df


def classify_impact(consequence):
    """Classify variant impact from consequence string."""
    c = str(consequence).lower()
    if any(x in c for x in [
        "stop_gained", "frameshift", "splice_donor",
        "splice_acceptor", "start_lost", "transcript_ablation"
    ]):
        return "HIGH"
    elif any(x in c for x in [
        "missense", "inframe_insertion", "inframe_deletion",
        "protein_altering", "splice_region"
    ]):
        return "MODERATE"
    elif any(x in c for x in [
        "synonymous", "stop_retained", "start_retained"
    ]):
        return "LOW"
    else:
        return "MODIFIER"


def categorise(consequence):
    """Group consequence into broad category for plotting."""
    c = str(consequence).lower()
    if "intergenic" in c:
        return "Intergenic"
    elif "upstream" in c:
        return "Upstream"
    elif "downstream" in c:
        return "Downstream"
    elif "intron" in c:
        return "Intronic"
    elif "missense" in c:
        return "Missense"
    elif "synonymous" in c:
        return "Synonymous"
    elif "frameshift" in c:
        return "Frameshift"
    elif "stop" in c:
        return "Stop variant"
    elif "splice" in c:
        return "Splice variant"
    elif "regulatory" in c:
        return "Regulatory"
    elif "utr" in c:
        return "UTR"
    else:
        return "Other"


# ══════════════════════════════════════════════════════════
# 3. SAVE CSV
# ══════════════════════════════════════════════════════════
def save_csv(df):
    print("\n[3/4] Saving CSV...")

    # Select most useful columns for the CSV output
    keep = [c for c in [
        "variant_id", "location", "alt_allele",
        "gene_symbol", "ensembl_gene", "transcript_id",
        "consequence", "consequence_category", "impact",
        "biotype", "exon", "intron",
        "hgvsc", "hgvsp",
        "sift", "sift_pred", "sift_score",
        "polyphen", "polyphen_pred", "polyphen_score",
        "rsid", "ref_allele", "af_1kg",
        "clinvar_sig", "mane_select", "distance",
    ] if c in df.columns]

    df[keep].to_csv(OUT_CSV, index=False)
    size_kb = os.path.getsize(OUT_CSV) // 1024
    print(f"  OK: Saved {OUT_CSV} ({size_kb} KB, {len(df)} rows)")


# ══════════════════════════════════════════════════════════
# 4. GENERATE PLOTS
# ══════════════════════════════════════════════════════════
def plot(df):
    print("\n[4/4] Generating plots...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "chr20 Variant Annotation Report\n"
        "Data: 1000 Genomes Project  |  Tool: Ensembl VEP v115  |  "
        "Assembly: GRCh38  |  Author: Farhana Sayed",
        fontsize=11, fontweight="bold", y=1.01
    )

    # ── Panel 1: Consequence categories ───────────────────
    ax = axes[0, 0]
    cat_colors = {
        "Intergenic":    "#0f3460",
        "Upstream":      "#16213e",
        "Downstream":    "#1a4a7a",
        "Intronic":      "#2b6cb0",
        "Regulatory":    "#4a90d9",
        "Missense":      "#f39c12",
        "Synonymous":    "#27ae60",
        "Frameshift":    "#e74c3c",
        "Stop variant":  "#c0392b",
        "Splice variant":"#e67e22",
        "UTR":           "#8e44ad",
        "Other":         "#888888",
    }
    if "consequence_category" in df.columns:
        cat = df["consequence_category"].value_counts()
        colors = [cat_colors.get(c, "#888") for c in cat.index]
        bars = ax.barh(cat.index[::-1], cat.values[::-1],
                       color=colors[::-1], edgecolor="white", height=0.65)
        ax.set_title("Consequence categories", fontweight="bold", fontsize=11)
        ax.set_xlabel("Number of annotations")
        for bar, val in zip(bars, cat.values[::-1]):
            ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
                    str(val), va="center", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # ── Panel 2: Impact distribution (donut) ──────────────
    ax = axes[0, 1]
    if "impact" in df.columns:
        ic = df["impact"].replace("", "MODIFIER").value_counts()
        pie_col = {
            "HIGH":     "#e74c3c",
            "MODERATE": "#f39c12",
            "LOW":      "#27ae60",
            "MODIFIER": "#3498db",
        }
        colors = [pie_col.get(i, "#888") for i in ic.index]
        wedges, texts, auto = ax.pie(
            ic.values, labels=None, colors=colors,
            autopct="%1.1f%%", startangle=90,
            pctdistance=0.75,
            wedgeprops={"width": 0.5, "edgecolor": "white", "linewidth": 2}
        )
        for at in auto:
            at.set_fontsize(9)
        ax.set_title("Impact distribution", fontweight="bold", fontsize=11)
        # Centre label
        ax.text(0, 0, f"n={len(df)}\nvariants",
                ha="center", va="center", fontsize=9, fontweight="bold")
        # Legend
        legend_patches = [
            mpatches.Patch(color=pie_col.get(k, "#888"), label=f"{k} ({ic.get(k, 0)})")
            for k in ["HIGH", "MODERATE", "LOW", "MODIFIER"] if k in ic.index
        ]
        ax.legend(handles=legend_patches, loc="lower center",
                  bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9)

    # ── Panel 3: Genes affected ────────────────────────────
    ax = axes[1, 0]
    if "gene_symbol" in df.columns:
        genes = df["gene_symbol"].replace("", "intergenic").value_counts()
        gene_colors = ["#e94560" if g != "intergenic" else "#888"
                       for g in genes.index]
        bars = ax.bar(genes.index, genes.values,
                      color=gene_colors, edgecolor="white", width=0.6)
        ax.set_title("Annotations per gene", fontweight="bold", fontsize=11)
        ax.set_ylabel("Annotation count")
        ax.tick_params(axis="x", rotation=30)
        for bar, val in zip(bars, genes.values):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.1,
                    str(val), ha="center", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        # Add legend
        legend_patches = [
            mpatches.Patch(color="#e94560", label="Protein-coding gene"),
            mpatches.Patch(color="#888",    label="Intergenic"),
        ]
        ax.legend(handles=legend_patches, fontsize=8)

    # ── Panel 4: Biotype breakdown ─────────────────────────
    ax = axes[1, 1]
    if "biotype" in df.columns:
        bt = df["biotype"].replace("", "intergenic").value_counts()
        bt_colors = {
            "protein_coding":               "#0f3460",
            "protein_coding_CDS_not_defined":"#2b6cb0",
            "intergenic":                   "#888",
            "enhancer":                     "#8e44ad",
        }
        colors = [bt_colors.get(b, "#4a90d9") for b in bt.index]
        # Shorten long biotype labels
        labels = [b.replace("protein_coding_CDS_not_defined",
                             "protein_coding\n(CDS not defined)")
                  for b in bt.index]
        ax.barh(labels[::-1], bt.values[::-1],
                color=colors[::-1], edgecolor="white", height=0.6)
        ax.set_title("Transcript biotype", fontweight="bold", fontsize=11)
        ax.set_xlabel("Count")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for i, v in enumerate(bt.values[::-1]):
            ax.text(v + 0.1, i, str(v), va="center", fontsize=9)
    else:
        ax.text(0.5, 0.5, "No biotype data", ha="center", va="center")

    plt.tight_layout(pad=2.5)
    plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close()
    size_kb = os.path.getsize(OUT_PNG) // 1024
    print(f"  OK: Saved {OUT_PNG} ({size_kb} KB)")


# ══════════════════════════════════════════════════════════
# 5. HTML REPORT
# ══════════════════════════════════════════════════════════
def html_report(df):
    print("\n[+] Generating HTML report...")

    # Summary stats
    n_total    = len(df)
    n_genes    = df["gene_symbol"].replace("", pd.NA).dropna().nunique() \
                 if "gene_symbol" in df.columns else 0
    n_unique_v = df["variant_id"].nunique() \
                 if "variant_id" in df.columns else n_total
    n_mane     = len(df[df["mane_select"].str.strip() != ""]) \
                 if "mane_select" in df.columns else 0

    # Build variant table — show key columns only
    show = [c for c in [
        "variant_id", "location", "gene_symbol",
        "consequence", "impact",
        "sift_pred", "polyphen_pred",
        "rsid", "biotype", "mane_select",
    ] if c in df.columns]

    table_html = (df[show]
                  .to_html(index=False, border=0,
                           classes="vtable", escape=True))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>chr20 Annotation Report — Farhana Sayed</title>
<style>
*{{box-sizing:border-box;margin:0;padding:0}}
body{{font-family:'Segoe UI',sans-serif;background:#f4f6fb;color:#1a1a2e}}
.hdr{{background:linear-gradient(135deg,#0f3460,#16213e);color:white;
      padding:40px 48px}}
.hdr h1{{font-size:22px;margin-bottom:8px}}
.hdr p{{font-size:13px;opacity:.78}}
.body{{padding:32px 48px}}
.cards{{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;
        margin-bottom:28px}}
.card{{background:white;border-radius:10px;padding:20px;text-align:center;
       box-shadow:0 2px 8px rgba(0,0,0,.06)}}
.card .num{{font-size:30px;font-weight:700;color:#0f3460}}
.card .lbl{{font-size:11px;color:#888;text-transform:uppercase;
            letter-spacing:.5px;margin-top:4px}}
h2{{font-size:16px;color:#0f3460;border-left:4px solid #e94560;
    padding-left:12px;margin:24px 0 14px}}
img{{width:100%;border-radius:10px;
     box-shadow:0 2px 12px rgba(0,0,0,.08);margin-bottom:24px}}
.vtable{{width:100%;border-collapse:collapse;font-size:12px;
         background:white;border-radius:10px;overflow:hidden;
         box-shadow:0 2px 8px rgba(0,0,0,.06)}}
.vtable th{{background:#0f3460;color:white;padding:9px 12px;
            text-align:left;font-size:11px}}
.vtable td{{padding:8px 12px;border-bottom:1px solid #eee;
            vertical-align:top}}
.vtable tr:nth-child(even){{background:#f8f9fd}}
.vtable tr:hover{{background:#eef2ff}}
.footer{{background:#0f3460;color:white;text-align:center;
         padding:18px;font-size:12px;margin-top:24px;opacity:.9}}
@media(max-width:700px){{
  .cards{{grid-template-columns:1fr 1fr}}
  .hdr,.body{{padding:20px 16px}}
}}
</style>
</head>
<body>

<div class="hdr">
  <h1>🧬 chr20 Variant Annotation Report</h1>
  <p>Data: 1000 Genomes Project &nbsp;|&nbsp;
     Tool: Ensembl VEP v115 &nbsp;|&nbsp;
     Assembly: GRCh38 &nbsp;|&nbsp;
     Generated: {datetime.now().strftime("%d %b %Y %H:%M")}</p>
  <p style="margin-top:6px">
     Author: Farhana Sayed &nbsp;|&nbsp;
     B.Tech Bioinformatics &amp; Data Science, D.Y. Patil &nbsp;|&nbsp;
     Research: Tata Memorial Centre – ACTREC</p>
</div>

<div class="body">

  <div class="cards">
    <div class="card">
      <div class="num">{n_unique_v}</div>
      <div class="lbl">Unique variants</div>
    </div>
    <div class="card">
      <div class="num">{n_total}</div>
      <div class="lbl">Total annotations</div>
    </div>
    <div class="card">
      <div class="num">{n_genes}</div>
      <div class="lbl">Genes affected</div>
    </div>
    <div class="card">
      <div class="num">{n_mane}</div>
      <div class="lbl">MANE Select rows</div>
    </div>
  </div>

  <h2>Summary Plots</h2>
  <img src="chr20_plots.png" alt="Variant Summary Plots">

  <h2>Annotation Table ({n_total} rows)</h2>
  {table_html}

</div>

<div class="footer">
  chr20 Variant Annotation Pipeline &nbsp;|&nbsp;
  Farhana Sayed Portfolio Project 5 &nbsp;|&nbsp;
  Tools: Ensembl VEP · Python · pandas · matplotlib
</div>
</body>
</html>"""

    with open(OUT_HTM, "w") as f:
        f.write(html)
    size_kb = os.path.getsize(OUT_HTM) // 1024
    print(f"  OK: Saved {OUT_HTM} ({size_kb} KB)")


# ══════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print("  chr20 Variant Annotation Pipeline — Farhana Sayed")
    print("  Data   : 1000 Genomes Project chr20 subset")
    print("  Tool   : Ensembl VEP v115 | GRCh38")
    print("  Input  :", VEP_IN)
    print("=" * 60)

    df_raw = parse_vep(VEP_IN)
    df     = process(df_raw)
    save_csv(df)
    plot(df)
    html_report(df)

    # ── Print final summary table ──────────────────────────
    print("\n" + "=" * 60)
    print("  PIPELINE COMPLETE — RESULTS SUMMARY")
    print("=" * 60)
    print(f"  Unique variants   : {df['variant_id'].nunique()}")
    print(f"  Total annotations : {len(df)}")

    if "impact" in df.columns:
        print("\n  Impact breakdown:")
        for level, color in [("HIGH","!!"), ("MODERATE","! "),
                              ("LOW","  "), ("MODIFIER","  ")]:
            n = len(df[df["impact"] == level])
            if n:
                print(f"    {color} {level:<10} : {n}")

    if "consequence_category" in df.columns:
        print("\n  Consequence categories:")
        for cat, n in df["consequence_category"].value_counts().items():
            print(f"    {cat:<20} : {n}")

    if "gene_symbol" in df.columns:
        genes = df["gene_symbol"].replace("", "intergenic").unique()
        print(f"\n  Genes/regions     : {', '.join(genes)}")

    print("\n  Output files:")
    for fpath in [OUT_CSV, OUT_PNG, OUT_HTM]:
        kb = os.path.getsize(fpath) // 1024 if os.path.exists(fpath) else 0
        print(f"    {fpath}  ({kb} KB)")
    print()


if __name__ == "__main__":
    main()
