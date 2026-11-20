#!/usr/bin/env python3
"""Who wins the MISASSIGNED reads? For each target taxid in --taxa, over every sample's saved BLAST TSV,
classify each read at --margin (same rule as classify_one_read) and, for MISASSIGNED reads, record the
top NAMED hit organism (the "winner"), whether its name came from sscinames (NCBI taxonomy) or the
stitle fallback (submitter free text), and the raw stitle. Also dumps every BLAST row for a few example
reads of --example_taxid (CONFIRMED and MISASSIGNED), to show what a read looks like.

  python tools/misassigned_winners.py --blast_dir DIR --taxid_list sample_taxid.tsv --taxa taxa.tsv \
      --out_prefix /n/netscratch/.../misassigned --example_taxid 54291
taxa.tsv: taxid<TAB>organism_name (no header).
Outputs: <prefix>_winners.tsv (per target x winner counts), <prefix>_example_reads.tsv.
"""
import os, sys, argparse
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))
from false_positive_detection import COLS, is_uncultured, organism_matches, organism_from_stitle  # noqa: E402


def name_of(r):
    s = r["sscinames"]
    if pd.notna(s) and str(s).strip() not in ("N/A", ""):
        return s, "sscinames"
    return organism_from_stitle(r["stitle"]), "stitle"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast_dir", required=True)
    ap.add_argument("--taxid_list", required=True)
    ap.add_argument("--taxa", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--margin", type=float, default=0.02)
    ap.add_argument("--id_min", type=float, default=80.0)
    ap.add_argument("--example_taxid", default="54291")
    ap.add_argument("--n_examples", type=int, default=6)
    a = ap.parse_args()

    taxa = pd.read_csv(a.taxa, sep="\t", header=None, names=["taxid", "name"], dtype=str)
    want = dict(zip(taxa.taxid, taxa.name))
    pairs = pd.read_csv(a.taxid_list, sep="\t", header=None, names=["sample", "taxid"], dtype=str).dropna()
    samples = sorted(pairs[pairs.taxid.isin(want)]["sample"].unique())

    win_rows, ex_rows, ex_seen = [], [], {"CONFIRMED": 0, "MISASSIGNED": 0}
    for i, smp in enumerate(samples, 1):
        p = os.path.join(a.blast_dir, f"{smp}_blast_results.tsv")
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            continue
        df = pd.read_csv(p, sep="\t", header=None, names=COLS, dtype=str, on_bad_lines="skip")
        df["taxid"] = df["qseqid"].str.split("|", n=1).str[0]
        df = df[df.taxid.isin(want)].copy()
        if df.empty:
            continue
        df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
        df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
        nm = df.apply(name_of, axis=1, result_type="expand")
        df["hit_organism"], df["name_source"] = nm[0], nm[1]
        for (tx, q), h in df.groupby(["taxid", "qseqid"], sort=False):
            exp = want[tx]
            named = h[(~h.hit_organism.map(is_uncultured)) & h.bitscore.notna()]
            if named.empty:
                continue  # EXCLUDED
            top = named.bitscore.max()
            tgt = named[named.hit_organism.map(lambda o: organism_matches(o, exp))]
            if tgt.empty:
                lab = "MISASSIGNED"
            else:
                r = tgt.loc[tgt.bitscore.idxmax()]
                ok = r.bitscore >= (1 - a.margin) * top and pd.notna(r.pident) and r.pident >= a.id_min
                lab = "CONFIRMED" if ok else "MISASSIGNED"
            if lab == "MISASSIGNED":
                w = named.loc[named.bitscore.idxmax()]
                win_rows.append({"taxid": tx, "expected_organism": exp, "winner": w.hit_organism,
                                 "winner_source": w.name_source, "winner_stitle": str(w.stitle)[:120],
                                 "target_present": not tgt.empty,
                                 "ratio": (tgt.bitscore.max() / top) if not tgt.empty else None})
            if tx == a.example_taxid and ex_seen[lab] < a.n_examples // 2:
                ex_seen[lab] += 1
                hh = h.sort_values("bitscore", ascending=False).assign(sample=smp, read_label=lab)
                ex_rows.append(hh[["sample", "qseqid", "read_label", "hit_organism", "name_source", "sscinames",
                                   "pident", "length", "evalue", "bitscore", "stitle"]])
        print(f"[{i}/{len(samples)}] {smp}", flush=True)

    W = pd.DataFrame(win_rows)
    W.to_csv(a.out_prefix + "_reads.tsv", sep="\t", index=False)
    agg = (W.groupby(["taxid", "expected_organism", "winner", "winner_source"])
             .agg(reads=("winner", "size"), target_present=("target_present", "mean"),
                  median_ratio=("ratio", "median"), example_stitle=("winner_stitle", "first"))
             .reset_index().sort_values(["expected_organism", "reads"], ascending=[True, False]))
    agg.to_csv(a.out_prefix + "_winners.tsv", sep="\t", index=False)
    if ex_rows:
        pd.concat(ex_rows).to_csv(a.out_prefix + "_example_reads.tsv", sep="\t", index=False)
    print("done ->", a.out_prefix + "_{reads,winners,example_reads}.tsv")


if __name__ == "__main__":
    main()
