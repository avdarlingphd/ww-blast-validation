#!/usr/bin/env python3
"""Read-level overview for one taxon under the final rules (2% margin, %id >= 90, non-organism and plasmid hits
ignored): what happens to every read, and what the reads removed by the 60 bp rule matched. Read level only:
no minimum-confirmed-reads rule (that is applied downstream in the Time-Series R script).

  python docs/make_failed_reads_overview.py DATA_DIR OUT_DIR TAXID [TAXID ...]

DATA_DIR  rr_<taxid>_reads.tsv (+ rr_<taxid>_examples.tsv for the name) from the cluster job rule_reads.sbatch
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data_dir, out_dir, *taxids = sys.argv[1:]
os.makedirs(out_dir, exist_ok=True)
BLUE, ORANGE, GREY, INK, INK2 = "#2a78d6", "#eb6834", "#b9b8b3", "#0b0b0b", "#52514e"
LABCOL = {"CONFIRMED": "#3a9d5d", "EXCLUDED_SHORT": ORANGE, "EXCLUDED_PLASMID": "#8e7cc3",
          "MISASSIGNED": "#c0392b", "EXCLUDED_COMPLEX": "#7f9cc6", "EXCLUDED": GREY}
plt.rcParams.update({"font.family": "DejaVu Sans", "axes.edgecolor": "#8a8984", "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})


def short(n):
    p = str(n).split()
    if len(p) >= 2 and p[0][:1].isupper() and p[0].isalpha() and p[1].islower() and p[1] != "sp." and p[1].isalpha():
        return f"{p[0][0]}. {' '.join(p[1:])}"
    return str(n)


for t in taxids:
    R = pd.read_csv(f"{data_dir}/rr_{t}_reads.tsv", sep="\t")
    nm = pd.read_csv(f"{data_dir}/rr_{t}_examples.tsv", sep="\t")["target"].iloc[0]
    conf0 = R[R["lab0"] == "CONFIRMED"].copy()
    conf0["removed"] = conf0["lab60"] == "EXCLUDED_SHORT"
    fig, ax = plt.subplots(1, 4, figsize=(21, 4.6), gridspec_kw={"width_ratios": [0.9, 1.1, 1.1, 1.2], "wspace": 0.45})

    # A: every read, label with no length cutoff vs at 60 bp
    a = ax[0]
    order = ["CONFIRMED", "EXCLUDED_SHORT", "EXCLUDED_COMPLEX", "EXCLUDED_PLASMID", "EXCLUDED", "MISASSIGNED"]
    for i, col in enumerate(["lab0", "lab60"]):
        c = R[col].value_counts()
        left = 0
        for lab in order:
            v = int(c.get(lab, 0))
            if v:
                a.barh(i, v, left=left, color=LABCOL[lab], edgecolor="white", lw=1.5, height=0.6)
                if v / len(R) > 0.06:
                    a.text(left + v / 2, i, str(v), ha="center", va="center", fontsize=8.5, color="white", fontweight="bold")
                left += v
    a.set_yticks([0, 1]); a.set_yticklabels(["no length\ncutoff", "60 bp"]); a.invert_yaxis()
    a.set_xlabel(f"reads ({len(R):,} with a BLAST hit, all samples)")
    a.set_title("what happens to every read", loc="left", fontsize=10, color=INK2)
    hs = [plt.Rectangle((0, 0), 1, 1, color=LABCOL[l]) for l in order if l in set(R["lab0"]) | set(R["lab60"])]
    a.legend(hs, [l for l in order if l in set(R["lab0"]) | set(R["lab60"])], fontsize=7.5, frameon=False,
             loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2)

    # B: CONFIRMED reads (no cutoff): aligned length vs species within 2%
    b = ax[1]
    jit = np.random.default_rng(1).uniform(-0.25, 0.25, len(conf0))
    for rem, col, lab in [(False, "#3a9d5d", "kept at 60 bp"), (True, ORANGE, "removed at 60 bp")]:
        d = conf0[conf0["removed"] == rem]
        b.scatter(d["tlen"], d["n_sp_2pct"] + jit[conf0["removed"].values == rem], s=22, color=col,
                  edgecolor="white", lw=0.5, label=f"{lab} ({len(d)})")
    b.axvline(60, color=INK2, ls=":", lw=0.9)
    b.set_xlabel("target alignment length (bp of 150)"); b.set_ylabel("species tied within 2% (incl. target)")
    b.set_xlim(40, 155); b.legend(fontsize=8, frameon=False, loc="upper right")
    b.set_title("reads CONFIRMED with no length cutoff", loc="left", fontsize=10, color=INK2)

    # C: which species tie on the removed reads
    c = ax[2]
    rem = conf0[conf0["removed"]]
    oth = rem["sp_2pct"].fillna("").str.split(";").apply(lambda L: [s for s in L if s and s != nm])
    cnt = pd.Series([s for L in oth for s in L]).value_counts().head(8)
    alone = (oth.apply(len) == 0).mean() if len(rem) else 0
    labels = ["target only"] + [short(s) for s in cnt.index]
    vals = [alone] + list(cnt.values / max(len(rem), 1))
    y = np.arange(len(labels))[::-1]
    c.barh(y, vals, color=[BLUE] + [ORANGE] * len(cnt), height=0.62)
    for yy, v in zip(y, vals):
        c.text(v + 0.01, yy, f"{v:.0%}", va="center", fontsize=8, color=INK2)
    c.set_yticks(y); c.set_yticklabels(labels, fontsize=8.2); c.set_xlim(0, 1.15)
    c.set_xlabel(f"share of the {len(rem)} reads removed at 60 bp")
    c.set_title("species tied within 2% on removed reads", loc="left", fontsize=10, color=INK2)

    # D: per sample, CONFIRMED reads with no cutoff vs at 60 bp (no minimum applied)
    d = ax[3]
    s0 = R[R["lab0"] == "CONFIRMED"].groupby("sample").size()
    s60 = R[R["lab60"] == "CONFIRMED"].groupby("sample").size().reindex(s0.index, fill_value=0)
    o = s0.sort_values(ascending=False).index
    x = np.arange(len(o))
    d.bar(x - 0.2, s0[o], width=0.4, color=GREY, label="no length cutoff")
    d.bar(x + 0.2, s60[o], width=0.4, color="#3a9d5d", label="60 bp")
    d.set_xticks(x); d.set_xticklabels([s.replace("GKWWBAA", "") for s in o], rotation=90, fontsize=7)
    d.set_ylabel("CONFIRMED reads"); d.set_xlabel("sample (GKWWBAA...)")
    d.legend(fontsize=8, frameon=False)
    d.set_title(f"per sample: {len(s0)} samples with a CONFIRMED read, {int((s60 > 0).sum())} at 60 bp",
                loc="left", fontsize=10, color=INK2)

    fig.suptitle(f"{nm}: every read, final rules, no length cutoff vs 60 bp", x=0.08, ha="left",
                 fontsize=13, fontweight="bold", color=INK, y=1.04)
    for j, L in enumerate("ABCD"):
        ax[j].text(-0.14, 1.1, L, transform=ax[j].transAxes, fontsize=14, fontweight="bold", color=INK)
    fn = f"{out_dir}/failed_reads_overview_{nm.replace(' ', '_')}.png"
    fig.savefig(fn, dpi=200, bbox_inches="tight", facecolor="white"); plt.close(fig)
    open(fn.replace(".png", "_caption.txt"), "w").write(
        f"{nm}, read-level BLAST validation (v4 final rules: NCBI nt 2026_09, taxid matching, near-tie margin 2%, "
        "%id >= 90, uncultured / non-organism / plasmid hits ignored), all 250 samples, <= 500 reads per mate per "
        "detection. (A) Label of every read with a BLAST hit, with no alignment-length cutoff and at 60 bp. "
        "(B) Reads CONFIRMED with no cutoff: target alignment length vs number of species whose best hit scores "
        "within 2% of the top kept hit; orange = removed at 60 bp (EXCLUDED_SHORT). (C) Species tied within 2% on "
        "the removed reads (a read can list several). (D) CONFIRMED reads per sample; no minimum-read rule is "
        "applied here (that is an option in the downstream R workflow).\n")
    print("saved", fn, "| reads", len(R), "| confirmed no cutoff", len(conf0), "| removed at 60", len(rem),
          "| confirmed at 60", int((R["lab60"] == "CONFIRMED").sum()))
