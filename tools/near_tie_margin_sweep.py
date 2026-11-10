#!/usr/bin/env python3
"""Sweep the read-level near-tie margin (2% / 5% / 10%) over the saved per-sample BLAST TSVs.

The margin only changes ONE per-read test -- whether the target's best NAMED bitscore is within
(1 - margin) of the read's top NAMED bitscore -- so we compute that ratio once per read and threshold
it at every margin in a single pass. Vectorized (no per-read Python apply), so it is far faster than
re-running false_positive_detection.py once per margin.

Output: one tidy row per (sample, taxid, margin) with confirmed / misassigned / excluded counts,
validated_fraction = CONFIRMED / (CONFIRMED + MISASSIGNED + EXCLUDED), and supported (CONFIRMED > 2).
Optionally joins a clade_counts TSV (sample, taxid, reads_clade, total_reads) to add corrected and
relative abundance per margin. Reuses load_sample / is_uncultured / organism_matches from the pipeline.

  python near_tie_margin_sweep.py --blast_dir DIR --taxid_list TSV --taxid_names CSV --output OUT.tsv
      [--margins 0.02,0.05,0.10] [--id_min 80] [--clade_counts CC.tsv]
"""
import os, sys, argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))
from false_positive_detection import (load_sample, is_uncultured, organism_matches,  # noqa: E402
                                      per_read_frame_tax, label_reads, summarize_labels, load_synonyms)


def per_read_frame(sub: pd.DataFrame, expected: str) -> pd.DataFrame:
    """One row per read (qseqid) with the margin-independent pieces:
       has_named, has_target, ratio = target_best_bitscore / top_named_bitscore, tpid = target %id."""
    reads = pd.Index(sub["qseqid"].unique(), name="qseqid")
    # cache organism predicates over the (few) unique organism strings in this detection
    uorg = sub["hit_organism"].dropna().unique()
    uncult = {o: is_uncultured(o) for o in uorg}
    tgtm = {o: organism_matches(o, expected) for o in uorg}
    s = sub.assign(
        _uncult=sub["hit_organism"].map(uncult).fillna(False),
        _target=sub["hit_organism"].map(tgtm).fillna(False),
    )
    named = s[(~s["_uncult"]) & s["bitscore"].notna()]
    R = pd.DataFrame(index=reads)
    if named.empty:
        R["has_named"] = False
        R["has_target"] = False
        R["ratio"] = np.nan
        R["tpid"] = np.nan
        return R
    top_named = named.groupby("qseqid")["bitscore"].max()
    tgt = named[named["_target"]].sort_values("bitscore", ascending=False, na_position="last")
    tgt_best = tgt.groupby("qseqid").first()  # bitscore + pident of the target's best named hit
    R["has_named"] = R.index.isin(top_named.index)
    R["has_target"] = R.index.isin(tgt_best.index)
    R["ratio"] = (tgt_best["bitscore"] / top_named).reindex(R.index)
    R["tpid"] = tgt_best["pident"].reindex(R.index)
    return R


def tally(R: pd.DataFrame, margin: float, id_min: float) -> tuple[int, int, int]:
    """(confirmed, misassigned, excluded) read counts at one margin."""
    n = len(R)
    if n == 0:
        return 0, 0, 0
    excluded = ~R["has_named"].to_numpy()
    ratio = R["ratio"].to_numpy(dtype=float)
    tpid = R["tpid"].to_numpy(dtype=float)
    ratio = np.where(np.isnan(ratio), -1.0, ratio)
    tpid = np.where(np.isnan(tpid), -1.0, tpid)
    confirmed = R["has_target"].to_numpy() & (ratio >= (1.0 - margin)) & (tpid >= id_min)
    nc = int(confirmed.sum())
    ne = int(excluded.sum())
    nm = n - nc - ne
    return nc, nm, ne


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast_dir", required=True)
    ap.add_argument("--taxid_list", required=True)
    ap.add_argument("--taxid_names", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--margins", default="0.02,0.05,0.10")
    ap.add_argument("--id_min", type=float, default=80.0)
    ap.add_argument("--min_confirmed", type=int, default=3, help="supported if CONFIRMED >= this (>2 => 3)")
    ap.add_argument("--clade_counts", default=None)
    ap.add_argument("--synonyms", default=None, help="resources/taxon_synonyms.csv (synonym-aware name matching)")
    ap.add_argument("--taxonomy", default=None,
                    help="v4 taxid mode (needs staxids in the BLAST TSVs); adds excluded_complex")
    a = ap.parse_args()
    if a.synonyms:
        print(f"synonyms: {load_synonyms(a.synonyms)} target names", flush=True)
    tax = None
    if a.taxonomy:
        from taxonomy import Taxonomy
        tax = Taxonomy(a.taxonomy)

    margins = [float(x) for x in a.margins.split(",")]
    taxid_to_name = {}
    if a.taxid_names and os.path.exists(a.taxid_names):
        tn = pd.read_csv(a.taxid_names)
        taxid_to_name = dict(zip(tn["taxid"].astype(str), tn["organism_name"]))
    clade = {}
    if a.clade_counts and os.path.exists(a.clade_counts):
        cc = pd.read_csv(a.clade_counts, sep="\t")
        cc.columns = [c.strip().lower() for c in cc.columns]
        for _, r in cc.iterrows():
            clade[(str(r["sample"]), str(r["taxid"]))] = (
                pd.to_numeric(r.get("reads_clade"), errors="coerce"),
                pd.to_numeric(r.get("total_reads"), errors="coerce"))
    pairs = pd.read_csv(a.taxid_list, sep="\t", header=None, names=["sample", "taxid"], dtype=str).dropna()
    pairs["taxid"] = pairs["taxid"].astype(str)

    rows = []
    n_samp = pairs["sample"].nunique()
    for si, (sample, grp) in enumerate(pairs.groupby("sample", sort=True), 1):
        df = load_sample(os.path.join(a.blast_dir, f"{sample}_blast_results.tsv"))
        if si % 25 == 0 or si == n_samp:
            print(f"[{si}/{n_samp}] {sample}", flush=True)
        for taxid in grp["taxid"]:
            expected = taxid_to_name.get(taxid, f"taxid_{taxid}")
            sub = df[df["taxid"] == taxid] if df is not None and not df.empty else None
            if tax is not None:
                R = per_read_frame_tax(sub, taxid, tax, expected) if sub is not None else None
            else:
                R = per_read_frame(sub, expected) if sub is not None else pd.DataFrame()
            cc = clade.get((sample, taxid))
            for m in margins:
                if tax is not None:
                    sm = (summarize_labels(label_reads(R, m, a.id_min), a.min_confirmed) if R is not None
                          else summarize_labels(pd.Series(dtype=str), a.min_confirmed))
                    nc, nm, ne, nx = (sm["rl_confirmed"], sm["rl_misassigned"], sm["rl_excluded"],
                                      sm["rl_excluded_complex"])
                else:
                    nc, nm, ne = tally(R, m, a.id_min); nx = 0
                denom = nc + nm + ne + nx
                frac = round(nc / denom, 4) if denom else None
                row = {"sample": sample, "taxid": taxid, "expected_organism": expected,
                       "margin": m, "confirmed": nc, "misassigned": nm, "excluded": ne,
                       "excluded_complex": nx, "reads_hit": denom, "validated_fraction": frac,
                       "validated_fraction_informative": round(nc / (nc + nm), 4) if (nc + nm) else None,
                       "supported": nc >= a.min_confirmed}
                if cc is not None:
                    reads_clade, total_reads = cc
                    if frac is not None and pd.notna(reads_clade):
                        corr = round(float(frac) * float(reads_clade), 1)
                        row["reads_clade"] = reads_clade
                        row["corrected_abundance"] = corr
                        row["relative_abundance"] = (round(corr / float(total_reads), 8)
                                                     if pd.notna(total_reads) and total_reads else None)
                rows.append(row)

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(a.output) or ".", exist_ok=True)
    out.to_csv(a.output, sep="\t", index=False)

    print("\nNEAR-TIE MARGIN SWEEP  (uncultured ignored, %%id >= %.0f)" % a.id_min)
    print(f"{'margin':>7} {'CONFIRMED':>11} {'MISASSIGNED':>12} {'EXCLUDED':>10} {'EXCL_CPLX':>10} "
          f"{'mean_vf':>8} {'supported':>10} {'detections':>11}")
    for m in margins:
        d = out[out["margin"] == m]
        vf = d["validated_fraction"].dropna()
        print(f"{m:>7.2f} {int(d['confirmed'].sum()):>11,} {int(d['misassigned'].sum()):>12,} "
              f"{int(d['excluded'].sum()):>10,} {int(d['excluded_complex'].sum()):>10,} {vf.mean():>8.3f} "
              f"{int(d['supported'].sum()):>10,} {len(d):>11,}")
    print(f"\nreport -> {a.output}")


if __name__ == "__main__":
    main()
