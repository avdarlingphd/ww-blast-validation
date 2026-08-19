#!/usr/bin/env python3
"""
BLAST false-positive detection (v2) — classify each (sample, taxid) candidate detection as
TRUE_POSITIVE / FALSE_POSITIVE / UNCERTAIN / UNCULTURED_DOMINANT / NO_DATA.

Detection-source agnostic: works on read sets extracted from Kraken2 OR Bowtie2 candidate calls;
the only difference upstream is how the per-(sample, taxid) FASTAs were produced.

v2 refinements (vs v1):
  * TOP-HIT-BY-BITSCORE: each read is judged by its single best hit (max bitscore), not "any hit among
    the returned set". Requires bitscore in the BLAST outfmt.
  * RESCUE high-%match / low-read: a pair with high pct_match but few matching reads is now TRUE_POSITIVE
    (few reads that all hit the right organism is confident evidence), instead of UNCERTAIN.
  * Pairs with -max_target_seqs 100 (set in the blast rule) so the true organism is not crowded out.

Usage (called by the Snakemake classify rule):
    python false_positive_detection.py --blast_dir DIR --taxid_list TSV --taxid_names CSV
        --output report.tsv [--threshold 80] [--min_match_reads 1]
"""
import os, re, argparse
import pandas as pd
from collections import Counter

# BLAST outfmt 6 columns (MUST match the blast rule): bitscore added for top-hit selection
COLS = ["qseqid", "sscinames", "pident", "length", "evalue", "bitscore", "stitle"]

_STOP = re.compile(r"\b(strain|subsp\.?|var\.|pv\.|chromosome|plasmid|DNA|RNA|scaffold|contig|complete|"
                   r"partial|whole|genome|sequence|assembly|node|clone|isolate|sp\.\s|cf\.\s)\b", re.I)

def organism_from_stitle(stitle: str) -> str:
    if not stitle or str(stitle).strip() in ("", "N/A", "nan"):
        return "Unknown"
    s = str(stitle).strip()
    m = _STOP.search(s)
    org = s[:m.start()].strip().rstrip(",;") if m else s
    w = org.split()
    return " ".join(w[:3]) if len(w) > 3 else org

def organism_matches(hit: str, expected: str) -> bool:
    if not hit or hit.lower() in ("unknown", "n/a", ""):
        return False
    h, e = hit.lower(), expected.lower()
    if e in h or h in e:
        return True
    ew = e.split()
    if len(ew) == 1:                       # genus-only expected -> genus match allowed
        hg = h.split()[0] if h.split() else ""
        return bool(hg and hg == ew[0])
    return False                            # species expected -> require species-level match

def is_uncultured(org: str) -> bool:
    if not org:
        return False
    lo = org.lower()
    return any(k in lo for k in ("uncultured", "unclassified", "environmental sample",
                                 "metagenome", "synthetic construct"))

def analyze_tsv(path: str, expected: str) -> dict:
    if not path or not os.path.exists(path):
        return {"status": "missing"}
    if os.path.getsize(path) == 0:
        return {"status": "empty"}
    try:
        df = pd.read_csv(path, sep="\t", header=None, names=COLS, dtype=str, on_bad_lines="skip")
    except Exception as e:
        return {"status": f"read_error: {e}"}
    if df.empty:
        return {"status": "empty"}
    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["bitscore"] = pd.to_numeric(df["bitscore"], errors="coerce")
    df["hit_organism"] = df.apply(
        lambda r: r["sscinames"] if pd.notna(r["sscinames"]) and str(r["sscinames"]).strip() not in ("N/A", "")
        else organism_from_stitle(r["stitle"]), axis=1)

    n_reads = df["qseqid"].nunique()
    # TOP HIT PER READ = highest bitscore (ties -> first). This is the v2 change.
    top = (df.sort_values("bitscore", ascending=False, na_position="last")
             .groupby("qseqid", sort=False).first().reset_index())

    n_match = int(top["hit_organism"].apply(lambda o: organism_matches(o, expected)).sum())
    pct_match = round(n_match / n_reads * 100, 1) if n_reads else 0.0

    counts = Counter(top["hit_organism"]).most_common()   # rank by TOP-hit read counts
    expected_rank = next((i for i, (o, _) in enumerate(counts, 1) if organism_matches(o, expected)), None)

    n_uncult = int(top["hit_organism"].apply(is_uncultured).sum())
    pct_uncult = round(n_uncult / n_reads * 100, 1) if n_reads else 0.0
    mean_pident = round(df["pident"].mean(), 2) if not df["pident"].isna().all() else None
    top_orgs = "; ".join(f"{o}({c})" for o, c in counts[:5])
    return {"status": "complete", "n_reads": n_reads, "n_match": n_match, "pct_match": pct_match,
            "pct_uncultured": pct_uncult, "expected_rank": expected_rank,
            "mean_pident": mean_pident, "top_organisms": top_orgs}

def classify(r1, r2, threshold=80.0, min_match_reads=1):
    valid = [s for s in (r1, r2) if s and s.get("status") == "complete" and s.get("n_reads", 0) > 0]
    if not valid:
        return "NO_DATA", None, None, None
    avg_pct = round(sum(s["pct_match"] for s in valid) / len(valid), 1)
    avg_unc = round(sum(s.get("pct_uncultured", 0) for s in valid) / len(valid), 1)
    ranks = [s.get("expected_rank") for s in valid]
    n_matches = [s.get("n_match", 0) for s in valid]
    if avg_unc > 50.0:                                        # uncultured/environmental dominates
        return "UNCULTURED_DOMINANT", avg_pct, avg_unc, ranks
    if avg_pct >= threshold:
        rank_ok = all(r is not None and r <= 2 for r in ranks)
        # v2 RESCUE: drop the strict n_match>2 gate. A high-%match pair with few reads (all matching)
        # is confident TRUE. Require only >= min_match_reads (default 1) matching read in each set.
        reads_ok = all(n >= min_match_reads for n in n_matches)
        return ("TRUE_POSITIVE" if (rank_ok and reads_ok) else "UNCERTAIN", avg_pct, avg_unc, ranks)
    if avg_pct < 10.0:
        return "FALSE_POSITIVE", avg_pct, avg_unc, ranks
    return "UNCERTAIN", avg_pct, avg_unc, ranks

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast_dir", required=True)
    ap.add_argument("--taxid_list", required=True)
    ap.add_argument("--taxid_names", default=None)
    ap.add_argument("--output", required=True)
    ap.add_argument("--threshold", type=float, default=80.0)
    ap.add_argument("--min_match_reads", type=int, default=1)
    a = ap.parse_args()

    taxid_to_name = {}
    if a.taxid_names and os.path.exists(a.taxid_names):
        tn = pd.read_csv(a.taxid_names)
        taxid_to_name = dict(zip(tn["taxid"].astype(str), tn["organism_name"]))
    samples = pd.read_csv(a.taxid_list, sep="\t", header=None, names=["sample", "taxid"])
    samples["taxid"] = samples["taxid"].astype(str)

    rows = []
    for _, row in samples.iterrows():
        sample, taxid = row["sample"], row["taxid"]
        expected = taxid_to_name.get(taxid, f"taxid_{taxid}")
        d = os.path.join(a.blast_dir, f"blast_validation_{sample}_{taxid}")
        r1 = analyze_tsv(os.path.join(d, f"{sample}_taxid_{taxid}_R1_blast_results.tsv"), expected)
        r2 = analyze_tsv(os.path.join(d, f"{sample}_taxid_{taxid}_R2_blast_results.tsv"), expected)
        label, avg_pct, avg_unc, ranks = classify(r1, r2, a.threshold, a.min_match_reads)
        if not os.path.isdir(d):
            label = "BLAST_NOT_RUN"
        rows.append({"sample": sample, "taxid": taxid, "expected_organism": expected,
                     "classification": label, "avg_pct_match": avg_pct, "avg_pct_uncultured": avg_unc,
                     "R1_n_reads": r1.get("n_reads"), "R1_n_match": r1.get("n_match"),
                     "R1_pct_match": r1.get("pct_match"), "R1_expected_rank": r1.get("expected_rank"),
                     "R1_top_organisms": r1.get("top_organisms"),
                     "R2_n_reads": r2.get("n_reads"), "R2_n_match": r2.get("n_match"),
                     "R2_pct_match": r2.get("pct_match"), "R2_expected_rank": r2.get("expected_rank"),
                     "R2_top_organisms": r2.get("top_organisms")})
    out = pd.DataFrame(rows)
    print("CLASSIFICATION SUMMARY")
    for lab, c in out["classification"].value_counts().items():
        print(f"  {lab:<22} {c:>5} ({c/len(out)*100:.1f}%)")
    os.makedirs(os.path.dirname(a.output) or ".", exist_ok=True)
    out.to_csv(a.output, sep="\t", index=False)
    print(f"\nreport -> {a.output}")

if __name__ == "__main__":
    main()
