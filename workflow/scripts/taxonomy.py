#!/usr/bin/env python3
"""NCBI taxonomy lookups for taxid-based read validation (v4).

Why: hit NAMES are unstable (NCBI renamed Raoultella -> Klebsiella and kept the taxids; nt titles carry
submitter names, "[Mycobacterium] chelonae", "MAG:" prefixes). Matching on taxids rolled up to species
is immune to all of that.

Table: tools/build_taxonomy_table.py builds it from a DATED NCBI taxdump that matches the nt snapshot
(nt 2026_09 -> taxdmp_2026-09-01). Columns taxid, parent, rank, name; rank == "MERGED" rows map an old
taxid to its current one.

Definitions:
  species_of(t)  the species-rank ancestor of t (t itself if it is a species); None above species.
  target_key(t)  the node a hit must descend from to count as the target: the species if t is at or
                 below species rank, else t itself (genus-level Kraken targets match at the genus).
  complex_of(t)  the BROADEST "species group"/"species subgroup" ancestor of t's species below the genus,
                 i.e. NCBI's own complex (e.g. 3390273 Klebsiella pneumoniae complex). None if not in one.
  is_uncultured(t)  name is uncultured/unclassified/metagenome/... or any ancestor is "environmental samples".
"""
import re
import pandas as pd

_COMPLEX_RANKS = {"species group", "species subgroup"}
_SPECIES_OR_BELOW = {"species", "subspecies", "strain", "serotype", "serogroup", "varietas", "forma",
                     "forma specialis", "isolate", "genotype", "biotype", "morph", "pathogroup"}
_UNCULT = re.compile(r"uncultured|unclassified|environmental sample|metagenome|synthetic construct", re.I)


class Taxonomy:
    def __init__(self, path: str):
        df = pd.read_csv(path, sep="\t", dtype={"taxid": "int64", "parent": "int64", "rank": str, "name": str},
                         keep_default_na=False)
        merged = df[df["rank"] == "MERGED"]
        self.merged = dict(zip(merged["taxid"], merged["parent"]))
        nodes = df[df["rank"] != "MERGED"]
        self.parent = dict(zip(nodes["taxid"], nodes["parent"]))
        self.rank = dict(zip(nodes["taxid"], nodes["rank"]))
        self.name = dict(zip(nodes["taxid"], nodes["name"]))
        self._lin, self._sp, self._cx, self._unc = {}, {}, {}, {}

    # ---- basics ----
    def resolve(self, t):
        """int taxid (merged -> current); None for missing/0/unknown."""
        try:
            t = int(t)
        except (TypeError, ValueError):
            return None
        t = self.merged.get(t, t)
        return t if t in self.parent and t != 0 else None

    def lineage(self, t):
        """tuple of taxids from t up to the root (inclusive)."""
        if t in self._lin:
            return self._lin[t]
        out, cur, seen = [], t, set()
        while cur is not None and cur not in seen:
            out.append(cur); seen.add(cur)
            p = self.parent.get(cur)
            if p is None or p == cur:
                break
            cur = p
        self._lin[t] = tuple(out)
        return self._lin[t]

    # ---- derived ----
    def species_of(self, t):
        if t in self._sp:
            return self._sp[t]
        sp = next((a for a in self.lineage(t) if self.rank.get(a) == "species"), None)
        self._sp[t] = sp
        return sp

    def target_key(self, t):
        t = self.resolve(t)
        if t is None:
            return None
        return self.species_of(t) if self.rank.get(t) in _SPECIES_OR_BELOW else t

    def complex_of(self, t):
        if t in self._cx:
            return self._cx[t]
        cx = None
        sp = self.species_of(t)
        if sp is not None:
            for a in self.lineage(sp)[1:]:      # walk up from the species, stop at the genus
                r = self.rank.get(a)
                if r in ("genus", "subgenus"):
                    break
                if r in _COMPLEX_RANKS:
                    cx = a                      # keep walking: broadest complex below the genus
        self._cx[t] = cx
        return cx

    def is_uncultured(self, t):
        if t in self._unc:
            return self._unc[t]
        u = bool(_UNCULT.search(self.name.get(t, ""))) or any(
            self.name.get(a, "") == "environmental samples" for a in self.lineage(t))
        self._unc[t] = u
        return u


def parse_staxids(v):
    """BLAST staxids cell ('54291', '54291;1234', 'N/A', '0') -> list of ints (may be empty)."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return []
    out = []
    for x in str(v).split(";"):
        x = x.strip()
        if x.isdigit() and x != "0":
            out.append(int(x))
    return out
