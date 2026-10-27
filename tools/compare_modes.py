#!/usr/bin/env python3
"""Compare matching modes on the SAME saved BLAST TSVs in a single pass.

Reads each per-sample <sample>_blast_results.tsv ONCE and scores every (sample, taxid) under both the
default (single best-hit) and tie-aware matching rules, then prints the before/after comparison and
writes a per-detection table. Free: no re-BLAST. Single-threaded (run with -c 1).

    python tools/compare_modes.py --blast_dir DIR --taxid_list TSV --taxid_names CSV --output OUT.tsv
"""
import os, argparse, importlib.util
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "fp", os.path.join(_HERE, "..", "workflow", "scripts", "false_positive_detection.py"))
fp = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(fp)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast_dir", required=True)
    ap.add_argument("--taxid_list", required=True)
    ap.add_argument("--taxid_names", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--threshold", type=float, default=80.0)
    ap.add_argument("--min_match_reads", type=int, default=3)
    a = ap.parse_args()

    names = pd.read_csv(a.taxid_names)
    t2n = dict(zip(names["taxid"].astype(str), names["organism_name"]))
    pairs = pd.read_csv(a.taxid_list, sep="\t", header=None, names=["sample", "taxid"], dtype=str).dropna()

    rows = []
    for i, (sample, grp) in enumerate(pairs.groupby("sample", sort=True), 1):
        df = fp.load_sample(os.path.join(a.blast_dir, f"{sample}_blast_results.tsv"))
        for taxid in grp["taxid"].astype(str):
            exp = t2n.get(taxid, f"taxid_{taxid}")
            if df is None:
                rows.append({"sample": sample, "taxid": taxid, "expected_organism": exp,
                             "class_default": "BLAST_NOT_RUN", "class_tie": "BLAST_NOT_RUN"})
                continue
            sub = df[df["taxid"] == taxid]
            r1, r2 = sub[sub["read"] == "R1"], sub[sub["read"] == "R2"]
            def label(**kw):
                return fp.classify(fp.analyze_reads(r1, exp, **kw), fp.analyze_reads(r2, exp, **kw),
                                   a.threshold, a.min_match_reads)[0]
            rows.append({"sample": sample, "taxid": taxid, "expected_organism": exp,
                         "class_default": label(), "class_tie": label(tie_aware=True)})
        if i % 25 == 0:
            print(f"  ...{i} samples", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(a.output, sep="\t", index=False)
    print("\n=== default ===\n" + out["class_default"].value_counts().to_string())
    print("\n=== tie_aware ===\n" + out["class_tie"].value_counts().to_string())
    mv = out[out["class_default"] != out["class_tie"]]
    print(f"\n=== changed: {len(mv)} detections ===")
    print(pd.crosstab(mv["class_default"], mv["class_tie"]).to_string())
    up = mv[(mv["class_default"].isin(["UNCERTAIN", "FALSE_POSITIVE"])) & (mv["class_tie"] == "TRUE_POSITIVE")]
    print("\n=== organisms rescued UNCERTAIN/FALSE -> TRUE_POSITIVE ===\n"
          + up["expected_organism"].value_counts().head(25).to_string())
    print(f"\nwrote {a.output}")

if __name__ == "__main__":
    main()
