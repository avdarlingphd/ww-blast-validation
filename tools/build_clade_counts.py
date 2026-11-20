#!/usr/bin/env python3
"""Build a clade_counts TSV (sample, taxid, reads_clade, total_reads) from the upstream Kraken2 reports,
for every (sample, taxid) in the detection list. Feeds the abundance columns of the read-level validator
and the near-tie sweep.

  reads_clade  = kreport field 2 for that taxid (reads in the clade rooted at the taxon, incl. descendants)
  total_reads  = reads_clade(taxid 0, unclassified) + reads_clade(taxid 1, root)
                 = total QC'd + dehosted reads Kraken saw for that sample (the library-size denominator)

  python build_clade_counts.py KRAKEN_DIR TAXID_LIST OUT.tsv
"""
import os, sys
import pandas as pd

kraken_dir, taxid_list, out = sys.argv[1], sys.argv[2], sys.argv[3]
pairs = pd.read_csv(taxid_list, sep="\t", header=None, names=["sample", "taxid"], dtype=str).dropna()
pairs["taxid"] = pairs["taxid"].astype(str)

cache = {}  # sample -> (dict taxid->reads_clade, total_reads)
missing = []


def load(sample):
    f = os.path.join(kraken_dir, f"{sample}.kreport")
    d = {}
    if not os.path.exists(f):
        missing.append(sample)
        return d, 0
    with open(f) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            try:
                d[p[4].strip()] = int(p[1].strip())
            except ValueError:
                continue
    return d, d.get("0", 0) + d.get("1", 0)


rows = []
for sample, grp in pairs.groupby("sample", sort=True):
    if sample not in cache:
        cache[sample] = load(sample)
    d, tot = cache[sample]
    for taxid in grp["taxid"]:
        rows.append({"sample": sample, "taxid": taxid,
                     "reads_clade": d.get(taxid, 0), "total_reads": tot})

df = pd.DataFrame(rows)
df.to_csv(out, sep="\t", index=False)
n_zero = int((df["reads_clade"] == 0).sum())
print(f"wrote {out}: {len(df)} rows, {df['sample'].nunique()} samples")
print(f"  reads_clade == 0 (taxid not in kreport): {n_zero}")
if missing:
    print(f"  kreport missing for {len(missing)} sample(s), e.g. {missing[:3]}")
