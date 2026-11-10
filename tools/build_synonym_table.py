#!/usr/bin/env python3
"""Build the taxon synonym table: every name a target taxon goes by, from NCBI + a manual override file.

  python tools/build_synonym_table.py TAXDUMP_DIR TARGETS.csv MANUAL.csv OUT.csv

TAXDUMP_DIR  nodes.dmp, names.dmp, merged.dmp of a DATED NCBI taxdump (taxdmp_2026-09-01 for nt 2026_09)
TARGETS.csv  taxid,organism_name (the Kraken target list, resources/taxid_to_name_all_pathogens.csv)
MANUAL.csv   taxid,synonym,note (resources/taxon_synonyms_manual.csv): names NCBI lacks, e.g. title-only
             spellings seen in nt
OUT.csv      one row per (taxid, synonym):
             taxid, kraken_name, ncbi_current_name, synonym, source
             source = kraken | ncbi_scientific | ncbi_synonym | ncbi_equivalent | ncbi_genus_combo | manual

Species synonyms come from NCBI name classes "scientific name", "synonym", "equivalent name" of the
target (merged taxids resolved; strain/subspecies targets also take their species' names).
"ncbi_genus_combo": the species epithet combined with every NCBI synonym of its GENUS (catches renames
where NCBI lists the genus synonym but not each species combination).
"""
import csv, os, re, sys

KEEP = {"scientific name": "ncbi_scientific", "synonym": "ncbi_synonym", "equivalent name": "ncbi_equivalent"}
BELOW_SPECIES = {"subspecies", "strain", "serotype", "serogroup", "varietas", "forma", "isolate", "genotype",
                 "biotype", "forma specialis", "morph", "pathogroup"}


# a synonym naming a SUBSPECIES / variety / pathovar of another species (e.g. "Aeromonas hydrophila subsp.
# anaerogenes" for A. caviae) must not be used: its first two words name a DIFFERENT species.
INFRA = re.compile(r"\b(subsp|ssp|var|pv|bv|serovar|biovar|forma|f\. sp|genomosp|genomovar)\b\.?", re.I)


def clean(n):
    n = re.sub(r"\s+", " ", re.sub(r"[\[\]\"']", "", n)).strip()
    return n


def key(n):
    """genus+epithet key, lowercase alphanumerics; None if not a binomial."""
    t = [w for w in (re.sub(r"[^a-z0-9]", "", x) for x in n.lower().split()) if w]
    return (t[0], t[1]) if len(t) >= 2 else None


def main(d, targets_csv, manual_csv, out_csv):
    targets = list(csv.DictReader(open(targets_csv)))
    want = set()
    merged = {}
    for line in open(os.path.join(d, "merged.dmp")):
        p = line.rstrip("\t|\n").split("\t|\t")
        merged[p[0]] = p[1]
    parent, rank = {}, {}
    for line in open(os.path.join(d, "nodes.dmp")):
        p = line.split("\t|\t")
        parent[p[0]], rank[p[0]] = p[1], p[2]

    def species_of(t):
        seen = set()
        while t in parent and t not in seen:
            if rank[t] == "species":
                return t
            seen.add(t); t = parent[t]
        return None

    def genus_of(t):
        seen = set()
        while t in parent and t not in seen:
            if rank[t] == "genus":
                return t
            seen.add(t); t = parent[t]
        return None

    _org = {}

    def is_organism(t):
        """True if t descends from cellular organisms (131567). False for vectors / "other sequences" and for
        viruses named after their host ("Aspergillus fumigatus partitivirus"), which NCBI also files at
        species rank and whose first two words would otherwise collide with the host species."""
        path, cur = [], t
        while cur not in _org:
            path.append(cur)
            p = parent.get(cur)
            if cur == "131567":
                _org[cur] = True; break
            if p is None or p == cur:
                _org[cur] = False; break
            cur = p
        v = _org[cur]
        for x in path:
            _org[x] = v
        return v

    resolved = {}
    for r in targets:
        t = merged.get(r["taxid"], r["taxid"])
        sp = species_of(t) if rank.get(t) in BELOW_SPECIES or rank.get(t) == "species" else None
        g = genus_of(t)
        resolved[r["taxid"]] = (t, sp, g)
        want.update(x for x in (t, sp, g) if x)

    # current scientific name of EVERY species: a synonym that is now another species' valid name is dropped
    # (NCBI keeps historical synonyms such as "Klebsiella aerogenes" under K. pneumoniae).
    names, species_key = {}, {}
    for line in open(os.path.join(d, "names.dmp")):
        p = line.rstrip("\t|\n").split("\t|\t")
        if p[3] == "scientific name" and rank.get(p[0]) == "species" and is_organism(p[0]):
            k = key(clean(p[1]))
            if k:
                species_key.setdefault(k, p[0])
        if p[0] in want and p[3] in KEEP and not (p[3] != "scientific name" and INFRA.search(p[1])):
            names.setdefault(p[0], []).append((clean(p[1]), KEEP[p[3]]))

    rows, seen, dropped = [], set(), []
    own_species = {}

    def add(taxid, kname, cur, syn, src):
        k = key(syn)
        other = species_key.get(k) if k else None
        if src not in ("kraken", "ncbi_scientific", "manual") and other and other != own_species.get(taxid):
            dropped.append({"taxid": taxid, "kraken_name": kname, "synonym": syn, "source": src,
                            "reason": f"current name of another species (taxid {other})"})
            return
        dk = (taxid, syn.lower())
        if syn and dk not in seen:
            seen.add(dk)
            rows.append({"taxid": taxid, "kraken_name": kname, "ncbi_current_name": cur, "synonym": syn, "source": src})

    for r in targets:
        tid, kname = r["taxid"], r["organism_name"]
        t, sp, g = resolved[tid]
        own_species[tid] = sp or t
        cur = next((n for n, s in names.get(t, []) if s == "ncbi_scientific"), "")
        add(tid, kname, cur, clean(kname), "kraken")
        for node in (t, sp):
            for n, s in names.get(node, []) if node else []:
                add(tid, kname, cur, n, s)
        # epithet x genus synonyms (e.g. Raoultella + ornithinolytica)
        sp_names = [n for n, s in names.get(sp or t, [])] + [clean(kname)]
        epithets = {n.split()[1] for n in sp_names if len(n.split()) >= 2}
        for gn, s in names.get(g, []) if g else []:
            gword = gn.split()[0]
            for ep in epithets:
                add(tid, kname, cur, f"{gword} {ep}", "ncbi_genus_combo")

    if manual_csv and os.path.exists(manual_csv):
        for m in csv.DictReader(open(manual_csv)):
            r = next((x for x in rows if x["taxid"] == m["taxid"]), None)
            add(m["taxid"], r["kraken_name"] if r else "", r["ncbi_current_name"] if r else "",
                clean(m["synonym"]), "manual")

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["taxid", "kraken_name", "ncbi_current_name", "synonym", "source"])
        w.writeheader(); w.writerows(rows)
    dpath = out_csv.replace(".csv", "_dropped.csv")
    with open(dpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["taxid", "kraken_name", "synonym", "source", "reason"])
        w.writeheader(); w.writerows(dropped)
    print(f"dropped {len(dropped)} synonyms that are now other species' names -> {dpath}")
    renamed = sum(1 for r in targets if (lambda c: c and c != clean(r["organism_name"]))(
        next((x["ncbi_current_name"] for x in rows if x["taxid"] == r["taxid"]), "")))
    print(f"{len(rows)} synonym rows for {len(targets)} targets; {renamed} targets renamed by NCBI -> {out_csv}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
