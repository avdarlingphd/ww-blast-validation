#!/usr/bin/env python3
"""BLAST output of reads removed by the 60 bp rule and by the 10% rule (final rules: 2% margin, %id >= 90,
non-organism and plasmid hits ignored).

  python docs/make_rule_reads_figures.py DATA_DIR OUT_DIR

DATA_DIR  rr_<taxid>_{reads,examples,detections}.tsv (cluster job rule_reads.sbatch)
OUT_DIR   rule_reads_examples_<taxon>.png, ten_percent_rule_summary.tsv, ten_percent_rule_summary.png (+ captions)
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data_dir, out_dir = sys.argv[1:3]
os.makedirs(out_dir, exist_ok=True)
BLUE, ORANGE, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#b9b8b3", "#0b0b0b", "#52514e"
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#8a8984", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})

# ---------------- example tables ----------------
for f in sorted(glob.glob(f"{data_dir}/rr_*_examples.tsv")):
    ex = pd.read_csv(f, sep="\t")
    nm = ex["target"].iloc[0]
    RR = pd.read_csv(f[:-len("_examples.tsv")] + "_reads.tsv", sep="\t")
    n0, n60 = int((RR.lab0 == "CONFIRMED").sum()), int((RR.lab60 == "CONFIRMED").sum())
    blocks = [(q, e["why"].iloc[0], e) for q, e in ex.groupby("qseqid", sort=False)]
    nrows = sum(len(b[2]) + 2.6 for b in blocks)
    fig, ax = plt.subplots(figsize=(14.5, 0.29 * nrows + 1.4))
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, nrows + 2.4)
    ax.text(0, nrows + 1.8, f"{nm}: BLAST output of reads removed by a rule", fontsize=12.5, fontweight="bold", color=INK)
    ax.text(0, nrows + 1.0, f"CONFIRMED reads (all samples): {n0} with no length cutoff, {n60} at 60 bp.  "
            "Best hit per species, highest bitscore first; blue row = target; grey = ignored hit.",
            fontsize=8.6, color=INK2)
    ax.text(0, nrows + 0.4, "'vs top' = bitscore / top kept bitscore (CONFIRMED needs >= 98%, %id >= 90 and, at 60 bp, "
            "the target aligned over >= 60 of 150 bp).", fontsize=8.6, color=INK2)
    xs = [0.0, 0.33, 0.39, 0.46, 0.53, 0.60]
    y = nrows - 0.6
    for q, why, e in blocks:
        smp, tag = q.split(":", 1)
        ax.text(0, y, f"{why}:  sample {smp}, read {tag.split('|', 1)[1]}   "
                f"({int(e['n_hits'].iloc[0])} hits, {int(e['n_species'].iloc[0])} species in total)",
                fontsize=9.2, fontweight="bold", color=INK2)
        y -= 1
        for h, x in zip(["hit species (NCBI taxid)", "%id", "aligned bp", "bitscore", "vs top", "title of best hit"], xs):
            ax.text(x, y, h, fontsize=8.4, fontweight="bold", color=INK2)
        kept = e[e["ignored"].isna()]
        top = (kept if len(kept) else e)["bitscore"].max()
        for _, r in e.iterrows():
            y -= 1
            is_t = r["hit_species"] == nm
            ign = isinstance(r["ignored"], str) and r["ignored"] != ""
            if is_t:
                ax.add_patch(plt.Rectangle((-0.005, y - 0.38), 1.0, 0.8, color="#cde2fb", lw=0, zorder=0))
            vals = [f"{r['hit_species']} ({str(r['staxids']).split(';')[0]})", f"{r['pident']:.1f}", f"{int(r['length'])}",
                    f"{r['bitscore']:.0f}", (f"ignored ({r['ignored']})" if ign else f"{r['bitscore'] / top:.1%}"),
                    str(r["stitle"])[:68]]
            for v, x in zip(vals, xs):
                ax.text(x, y, v, fontsize=8.1, color=GREY if ign else INK,
                        fontweight="bold" if is_t and x == 0 else "normal")
        y -= 1.6
    fn = f"{out_dir}/rule_reads_examples_{nm.replace(' ', '_')}.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("saved", fn)

# ---------------- 10% rule: which reads does it remove? ----------------
rows = []
for f in sorted(glob.glob(f"{data_dir}/rr_*_reads.tsv")):
    R = pd.read_csv(f, sep="\t")
    tid = os.path.basename(f).split("_")[1]
    exf = f[:-len("_reads.tsv")] + "_examples.tsv"
    nm = pd.read_csv(exf, sep="\t")["target"].iloc[0] if os.path.exists(exf) else tid
    conf = R[R["lab0"] == "CONFIRMED"]
    gone = conf[conf["lab10"] == "EXCLUDED_NONSPECIFIC"]
    kept = conf[conf["lab10"] == "CONFIRMED"]
    s0 = R[R["lab0"] == "CONFIRMED"].groupby("sample").size()
    s10 = R[R["lab10"] == "CONFIRMED"].groupby("sample").size()
    fill = gone["top_filler"].value_counts()
    rows.append(dict(
        taxon=nm, confirmed_reads=len(conf), removed_by_10pct=len(gone),
        removed_share=round(len(gone) / len(conf), 3) if len(conf) else None,
        removed_target_sole_best=round((gone["n_sp_2pct"] == 1).mean(), 3) if len(gone) else None,
        removed_other_genus_tie=round(gone["other_genus_tie"].astype(bool).mean(), 3) if len(gone) else None,
        kept_other_genus_tie=round(kept["other_genus_tie"].astype(bool).mean(), 3) if len(kept) else None,
        removed_median_target_frac=round(gone["tfrac"].median(), 3) if len(gone) else None,
        removed_hits_below_2pct=round(gone["share_hits_below_2pct"].median(), 3) if len(gone) else None,
        removed_list_filled_by="; ".join(f"{s} {c / len(gone):.0%}" for s, c in fill.head(3).items()) if len(gone) else "",
        samples_with_confirmed_0=int((s0 > 0).sum()), samples_with_confirmed_10pct=int((s10 > 0).sum())))
S = pd.DataFrame(rows)
S.to_csv(f"{out_dir}/ten_percent_rule_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(S.to_string(index=False))
