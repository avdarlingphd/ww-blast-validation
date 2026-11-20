#!/usr/bin/env python3
"""Build the compact NCBI taxonomy table used for taxid-based read validation.

  python tools/build_taxonomy_table.py TAXDUMP_DIR OUT.tsv

TAXDUMP_DIR holds nodes.dmp, names.dmp, merged.dmp from an NCBI taxdump (pin a dated archive that matches
the nt snapshot, e.g. taxdmp_2026-09-01.zip for nt 2026_09). OUT.tsv columns:
  taxid  parent  rank  name            (one row per node; name = scientific name)
  plus rows with rank == "MERGED" mapping an old taxid -> its current taxid (in `parent`).
"""
import os, sys


def main(d, out):
    names = {}
    with open(os.path.join(d, "names.dmp")) as f:
        for line in f:
            p = line.split("\t|\t")
            if p[3].rstrip("\t|\n") == "scientific name":
                names[p[0]] = p[1]
    n = 0
    with open(out, "w") as o:
        o.write("taxid\tparent\trank\tname\n")
        with open(os.path.join(d, "nodes.dmp")) as f:
            for line in f:
                p = line.split("\t|\t")
                o.write(f"{p[0]}\t{p[1]}\t{p[2]}\t{names.get(p[0], '')}\n"); n += 1
        m = 0
        with open(os.path.join(d, "merged.dmp")) as f:
            for line in f:
                p = line.rstrip("\t|\n").split("\t|\t")
                o.write(f"{p[0]}\t{p[1]}\tMERGED\t\n"); m += 1
    print(f"{n:,} nodes + {m:,} merged -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
