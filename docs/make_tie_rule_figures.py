#!/usr/bin/env python3
"""10% rule (decision #20), applied only when another species ties: compare the tie trigger
'another species within 2% of the top' with 'another species at exactly the top bitscore'.
Read level only (no minimum-confirmed-reads rule; that is an option in the Time-Series R script).

  python docs/make_tie_rule_figures.py DATA_DIR OUT_DIR [TAXID_FOR_EXAMPLES ...]

DATA_DIR  tr_<taxid>_tie_reads.tsv and tr_<taxid>_tie_examples.tsv (cluster job tie_reads)
OUT_DIR   tie_rule_summary.png/.tsv and tie_rule_examples_<taxon>.png (+ _caption.txt)
"""
import os, sys, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data_dir, out_dir, *ex_taxids = sys.argv[1:]
os.makedirs(out_dir, exist_ok=True)
GREY, INK, INK2 = "#b9b8b3", "#0b0b0b", "#52514e"
C_NONE, C_2, C_EX = "#b9b8b3", "#eb6834", "#2a78d6"
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#8a8984", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})

# ---------------- summary over taxa ----------------
rows = []
for f in sorted(glob.glob(f"{data_dir}/tr_*_tie_reads.tsv")):  # summary only for the tie_reads job
    R = pd.read_csv(f, sep="\t")
    exf = f.replace("_tie_reads", "_tie_examples")
    nm = pd.read_csv(exf, sep="\t")["target"].iloc[0] if os.path.exists(exf) else os.path.basename(f).split("_")[1]
    c = R[R.lab_none == "CONFIRMED"]
    d2 = c[c.lab_mid2 == "EXCLUDED_NONSPECIFIC"]
    rows.append(dict(taxon=nm, confirmed=len(c),
                     kept_2pct=int((c.lab_mid2 == "CONFIRMED").sum()), kept_exact=int((c.lab_exact == "CONFIRMED").sum()),
                     samples_none=c["sample"].nunique(),
                     samples_2pct=c.loc[c.lab_mid2 == "CONFIRMED", "sample"].nunique(),
                     samples_exact=c.loc[c.lab_exact == "CONFIRMED", "sample"].nunique(),
                     dropped_2pct_only=int((d2.lab_exact == "CONFIRMED").sum()),
                     dropped_2pct_only_top_filler="; ".join(f"{s} {n / max(len(d2[d2.lab_exact == 'CONFIRMED']), 1):.0%}" for s, n in
                                                            d2[d2.lab_exact == "CONFIRMED"]["best_other"].value_counts().head(3).items())))
S = pd.DataFrame(rows)
if len(S): S = S.sort_values("confirmed", ascending=False)
if len(S):
    S.to_csv(f"{out_dir}/tie_rule_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); print(S.to_string(index=False))

    fig, ax = plt.subplots(figsize=(9.5, 0.55 * len(S) + 1.6))
    y = np.arange(len(S))[::-1]
    for off, col, lab, key in [(0.25, C_NONE, "no 10% rule", None), (0, C_EX, "10% rule, exact tie", "kept_exact"),
                               (-0.25, C_2, "10% rule, tie within 2%", "kept_2pct")]:
        v = np.ones(len(S)) if key is None else (S[key] / S["confirmed"]).to_numpy()
        ax.barh(y + off, v, height=0.24, color=col, label=lab)
        if key:
            for yy, vv in zip(y + off, v):
                ax.text(vv + 0.01, yy, f"{vv:.0%}", va="center", fontsize=7.8, color=INK2)
    ax.set_yticks(y); ax.set_yticklabels([f"{t}  ({n:,})" for t, n in zip(S.taxon, S.confirmed)], fontsize=9)
    ax.set_xlim(0, 1.12); ax.set_xlabel("share of CONFIRMED reads (no 10% rule) still CONFIRMED")
    ax.legend(fontsize=8.5, frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=3)
    ax.set_title("10% rule: which tie triggers it?", loc="left", fontsize=12, fontweight="bold", color=INK, pad=24)
    fn = f"{out_dir}/tie_rule_summary.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    open(fn.replace(".png", "_caption.txt"), "w").write(
        "Read-level BLAST validation, v4 final rules (NCBI nt 2026_09, taxid matching, 2% near-tie margin, %id >= 90, uncultured / "
        "non-organism / plasmid hits ignored, no alignment-length cutoff), all 250 samples, <= 500 reads per mate per detection. "
        "Bars: share of reads CONFIRMED without the 10% rule that stay CONFIRMED when the rule (target >= 10% of the read's kept "
        "named hits, up to 100) is applied only to reads where another species scores within 2% of the top bitscore (orange) or "
        "exactly equal to it (blue). Label: taxon (CONFIRMED reads without the rule). No minimum-read rule is applied.\n")
    print("saved", fn)

# ---------------- example BLAST tables ----------------
for t in ex_taxids:
    cand = glob.glob(f"{data_dir}/*_{t}_*examples.tsv")
    if not cand:
        continue
    exf = cand[0]
    ex = pd.read_csv(exf, sep="\t"); nm = ex["target"].iloc[0]
    blocks = [(q, e["why"].iloc[0], e) for q, e in ex.groupby("qseqid", sort=False)]
    nrows = sum(len(b[2]) + 2.6 for b in blocks)
    fig, ax = plt.subplots(figsize=(15, 0.29 * nrows + 1.4))
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, nrows + 2.4)
    x15 = "n_top_hits_species" in ex.columns
    ax.text(0, nrows + 1.8, f"{nm}: reads around the target-share cutoff" if x15 else f"{nm}: reads the target-share rule removes",
            fontsize=12.5, fontweight="bold", color=INK)
    ax.text(0, nrows + 1.0, "Best hit per species, highest bitscore first; blue row = target; grey = ignored hit. "
            "'hits' = that species' kept hits in the list (up to 100); 'top hits' = those within 2% of the top.", fontsize=8.6, color=INK2)
    ax.text(0, nrows + 0.4, "Rule: the target must be >= X% of the kept hits scoring within 2% of the top. "
            "'vs top' = bitscore / top kept bitscore.", fontsize=8.6, color=INK2)
    xs = [0.0, 0.28, 0.33, 0.385, 0.44, 0.50, 0.545, 0.605] if x15 else [0.0, 0.29, 0.345, 0.405, 0.465, 0.525, 0.58]
    heads = ["hit species (NCBI taxid)", "%id", "aligned", "bitscore", "vs top", "hits"] + (["top hits"] if x15 else []) + ["title of best hit"]
    yy = nrows - 0.6
    for q, why, e in blocks:
        smp, tag = q.split(":", 1)
        ax.text(0, yy, f"{why}:  sample {smp}, read {tag.split('|', 1)[1]}   " +
                (f"   target share = {int(e.loc[e.hit_species == nm, 'n_top_hits_species'].max() if (e.hit_species == nm).any() else 0)} of "
                 f"{int(e['n_top'].iloc[0])} top hits ({int(e['n_top_species'].iloc[0])} species) = {e['share'].iloc[0]:.0%}" if x15 else
                 f"({int(e['n_kept'].iloc[0])} kept hits of {int(e['n_hits'].iloc[0])})"), fontsize=9.2, fontweight="bold", color=INK2)
        yy -= 1
        for h, x in zip(heads, xs):
            ax.text(x, yy, h, fontsize=8.4, fontweight="bold", color=INK2)
        kept = e[e["ignored"].isna()]
        top = (kept if len(kept) else e)["bitscore"].max()
        for _, r in e.iterrows():
            yy -= 1
            is_t = r["hit_species"] == nm
            ign = isinstance(r["ignored"], str) and r["ignored"] != ""
            if is_t:
                ax.add_patch(plt.Rectangle((-0.005, yy - 0.38), 1.0, 0.8, color="#cde2fb", lw=0, zorder=0))
            vals = [f"{r['hit_species']} ({str(r['staxids']).split(';')[0]})", f"{r['pident']:.1f}", f"{int(r['length'])}",
                    f"{r['bitscore']:.0f}", (f"ignored ({r['ignored']})" if ign else f"{r['bitscore'] / top:.1%}"),
                    ("" if ign else str(int(r["n_hits_species"])))] + \
                   ([("" if ign else str(int(r["n_top_hits_species"])))] if x15 else []) + [str(r["stitle"])[:68]]
            for v, x in zip(vals, xs):
                ax.text(x, yy, v, fontsize=8.1, color=GREY if ign else INK, fontweight="bold" if is_t and x == 0 else "normal")
        yy -= 1.6
    fn = f"{out_dir}/{'x15' if x15 else 'tie_rule'}_examples_{nm.replace(' ', '_')}.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("saved", fn)
