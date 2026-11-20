#!/usr/bin/env python3
"""Old (v2.1 report + v3 name-based read-level sweep) vs new (v4 taxid) per taxon, for review before
anything reaches the paper.

  python tools/compare_v2_v4.py OLD_REPORT OLD_SWEEP CLADE_COUNTS NEW_SWEEP OUT_PREFIX

Writes OUT_PREFIX_per_taxon.tsv (one row per taxon x margin) and OUT_PREFIX_per_detection.tsv
(detections whose supported call changed at any margin).
"""
import sys
import pandas as pd

old_rep, old_sw, cc, new_sw, out = sys.argv[1:6]
rep = pd.read_csv(old_rep, sep="\t", usecols=["sample", "taxid", "expected_organism", "classification"])
o = pd.read_csv(old_sw, sep="\t")
n = pd.read_csv(new_sw, sep="\t")
c = pd.read_csv(cc, sep="\t")
o = o.merge(c[["sample", "taxid", "reads_clade", "total_reads"]], on=["sample", "taxid"], how="left")
o["relative_abundance"] = o["validated_fraction"] * o["reads_clade"] / o["total_reads"]
n = n.merge(c[["sample", "taxid", "total_reads"]], on=["sample", "taxid"], how="left")
for d in (o, n):
    d["margin"] = d["margin"].round(2)
    d["supported"] = d["supported"].astype(str).str.lower().eq("true")

k = ["sample", "taxid", "margin"]
m = o[k + ["expected_organism", "misassigned", "supported", "validated_fraction", "relative_abundance"]].merge(
    n[k + ["misassigned", "excluded_complex", "supported", "validated_fraction", "relative_abundance"]],
    on=k, how="outer", suffixes=("_old", "_new"))
m = m.merge(rep.rename(columns={"classification": "v21_class", "expected_organism": "name_rep"}),
            on=["sample", "taxid"], how="left")
m["expected_organism"] = m["expected_organism"].fillna(m["name_rep"])
m["v21_TP"] = m["v21_class"].eq("TRUE_POSITIVE")
for c_ in ("supported_old", "supported_new"):
    m[c_] = m[c_].fillna(False).astype(bool)
nsamp = m["sample"].nunique()

g = m.groupby(["taxid", "expected_organism", "margin"])
t = pd.DataFrame({
    "detections": g.size(),
    "v21_true_positive": g["v21_TP"].sum(),
    "supported_old": g["supported_old"].sum(),
    "supported_new": g["supported_new"].sum(),
    "misassigned_reads_old": g["misassigned_old"].sum(),
    "misassigned_reads_new": g["misassigned_new"].sum(),
    "excluded_complex_reads_new": g["excluded_complex"].sum(),
    "mean_vf_old": g["validated_fraction_old"].mean().round(3),
    "mean_vf_new": g["validated_fraction_new"].mean().round(3),
    "mean_rel_abund_old_pct": (g["relative_abundance_old"].sum() / nsamp * 100),
    "mean_rel_abund_new_pct": (g["relative_abundance_new"].sum() / nsamp * 100),
}).reset_index()
t["supported_change"] = t["supported_new"] - t["supported_old"]
t["abund_ratio_new_old"] = (t["mean_rel_abund_new_pct"] / t["mean_rel_abund_old_pct"]).round(3)
t = t.sort_values(["margin", "supported_change", "misassigned_reads_old"], ascending=[True, False, False])
t.to_csv(f"{out}_per_taxon.tsv", sep="\t", index=False, float_format="%.6g")
ch = m[m["supported_old"] != m["supported_new"]].sort_values(["margin", "expected_organism", "sample"])
ch.to_csv(f"{out}_per_detection_changed.tsv", sep="\t", index=False, float_format="%.6g")

print(f"{nsamp} samples; detections per margin: {m.groupby('margin').size().to_dict()}")
for mg, s in m.groupby("margin"):
    print(f"margin {mg:.2f}: supported old {s.supported_old.sum():,} -> new {s.supported_new.sum():,}  "
          f"(gained {(~s.supported_old & s.supported_new).sum()}, lost {(s.supported_old & ~s.supported_new).sum()}); "
          f"misassigned reads {s.misassigned_old.sum():,.0f} -> {s.misassigned_new.sum():,.0f}; "
          f"v2.1 TRUE_POSITIVE {s.v21_TP.sum():,}")
