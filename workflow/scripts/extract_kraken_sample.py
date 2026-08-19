#!/usr/bin/env python3
"""
Batched, single-pass read extraction for ONE sample (kraken variant).

Reads the sample's ~5 GB Kraken2 output (.kraken) ONCE and its two FASTQs ONCE each, and writes the
per-(sample, taxid) R1/R2 FASTAs for EVERY taxid that sample was called for. This replaces looping
KrakenTools extract_kraken_reads.py per taxid (which reloaded the 5 GB .kraken once per taxid).

Behaviour is faithful to KrakenTools extract_kraken_reads.py with --include-children:
  * The .kreport is parsed into a taxonomy tree (2 spaces = 1 level); for each target taxid we collect
    its whole subtree, so reads classified to any strain/subspecies under the taxid are captured.
  * A read whose taxid falls under several targets (nested targets, e.g. genus + species both in the
    list) is written to EACH matching target — exactly as running the tool once per taxid would.
  * FASTQ read-IDs are matched on the first whitespace-delimited token, with a trailing /1 or /2
    stripped (same normalisation as KrakenTools).

Zero-read handling: every target taxid gets BOTH output FASTAs created, EMPTY if no reads matched, so
the DAG never breaks and the classifier reads them as NO_DATA. This is the case that produced the
original 1,415 NO_DATA calls.

Output layout (identical to the old per-pair rule):
  {out_dir}/blast_validation_{sample}_{taxid}/{sample}_taxid_{taxid}_R{1,2}.fasta

Usage:
  extract_kraken_sample.py --sample S --kraken K.kraken --kreport K.kreport
      --r1 S_R1.fastq.gz --r2 S_R2.fastq.gz --taxid_list list.tsv --out_dir FASTA_TREE
"""
import os, sys, gzip, argparse


# ---- KrakenTools-faithful report/kraken parsing (adapted from extract_kraken_reads.py) -----------
class Tree:
    def __init__(self, taxid, level_num, parent=None):
        self.taxid = taxid
        self.level_num = level_num
        self.parent = parent
        self.children = []


def process_kraken_report(report_line):
    l_vals = report_line.strip().split('\t')
    if len(l_vals) < 5:
        return []
    try:
        int(l_vals[1])
    except ValueError:
        return []
    try:
        taxid = int(l_vals[-3]); _ = l_vals[-2]        # kuniq-style layout
    except ValueError:
        taxid = int(l_vals[-2])                        # standard kraken2 .kreport
    spaces = 0
    for ch in l_vals[-1]:
        if ch == ' ':
            spaces += 1
        else:
            break
    return [taxid, int(spaces / 2)]


def process_kraken_output(kraken_line):
    l_vals = kraken_line.split('\t')
    if len(l_vals) < 5:
        return [-1, '']
    if "taxid" in l_vals[2]:
        tax_id = l_vals[2].split("taxid ")[-1][:-1]
    else:
        tax_id = l_vals[2]
    read_id = l_vals[1]
    tax_id = 81077 if tax_id == 'A' else int(tax_id)
    return [tax_id, read_id]


def subtree_taxids(kreport_path, targets):
    """Return {target_taxid: set(all taxids in its subtree, incl. itself)} using the kreport tree."""
    target_set = set(targets)
    base_nodes = {}
    prev = None
    with open(kreport_path) as fh:
        for line in fh:
            vals = process_kraken_report(line)
            if not vals:
                continue
            taxid, level_num = vals
            if taxid == 0:
                continue
            if taxid == 1:
                prev = Tree(taxid, level_num)
                if taxid in target_set:
                    base_nodes[taxid] = prev
                continue
            if prev is None:
                continue
            while level_num != (prev.level_num + 1) and prev.parent is not None:
                prev = prev.parent
            node = Tree(taxid, level_num, prev)
            prev.children.append(node)
            prev = node
            if taxid in target_set:
                base_nodes[taxid] = node
    # every target gets at least itself; add full descendant set where the node exists in the report
    result = {t: {t} for t in targets}
    for t, node in base_nodes.items():
        stack = list(node.children)
        while stack:
            n = stack.pop()
            result[t].add(n.taxid)
            stack.extend(n.children)
    return result


# ---- FASTQ -> FASTA streaming (manual; no BioPython) ---------------------------------------------
def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')


def norm_id(header_token):
    rid = header_token
    if rid.endswith('/1') or rid.endswith('/2'):
        rid = rid[:-2]
    return rid


def write_pass(fastq_path, read_to_targets, handles):
    """Stream one FASTQ, appending matching reads (as FASTA) to each target's handle."""
    n_out = 0
    with _open(fastq_path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            seq = fh.readline()
            fh.readline()   # +
            fh.readline()   # qual
            if not h.startswith('@'):
                continue
            tok = h[1:].split(None, 1)[0]
            targets = read_to_targets.get(norm_id(tok))
            if not targets:
                continue
            rec = ">" + tok + "\n" + seq.strip() + "\n"
            for t in targets:
                handles[t].write(rec)
            n_out += 1
    return n_out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--kraken", required=True)
    ap.add_argument("--kreport", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--taxid_list", required=True)
    ap.add_argument("--out_dir", required=True)
    a = ap.parse_args()

    # taxids this sample was called for
    targets = []
    with open(a.taxid_list) as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 2 and p[0] == a.sample:
                try:
                    targets.append(int(p[1]))
                except ValueError:
                    pass
    targets = sorted(set(targets))
    if not targets:
        sys.stderr.write(f"WARNING: no taxids for sample {a.sample} in {a.taxid_list}\n")

    # output paths + pre-create empty FASTAs so every (sample,taxid) always has R1/R2 (zero-read safe)
    def fa(taxid, read):
        d = os.path.join(a.out_dir, f"blast_validation_{a.sample}_{taxid}")
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, f"{a.sample}_taxid_{taxid}_{read}.fasta")
    r1_paths = {t: fa(t, "R1") for t in targets}
    r2_paths = {t: fa(t, "R2") for t in targets}
    for t in targets:
        open(r1_paths[t], 'w').close()
        open(r2_paths[t], 'w').close()

    if not targets:
        print(f"[{a.sample}] no targets; nothing to extract")
        return

    # STEP 0: kreport -> per-target subtree taxid sets; invert to taxid -> [targets]
    subtrees = subtree_taxids(a.kreport, targets)
    taxid_to_targets = {}
    for t, tids in subtrees.items():
        for d in tids:
            taxid_to_targets.setdefault(d, []).append(t)

    # STEP 1: single pass over the .kraken -> read_id -> [targets]
    read_to_targets = {}
    n_lines = 0
    with open(a.kraken) as fh:
        for line in fh:
            n_lines += 1
            tax_id, read_id = process_kraken_output(line)
            if tax_id == -1:
                continue
            tgts = taxid_to_targets.get(tax_id)
            if tgts:
                read_to_targets[read_id] = tgts
    print(f"[{a.sample}] kraken lines={n_lines:,}  targets={len(targets)}  "
          f"reads_of_interest={len(read_to_targets):,}")

    # STEP 2: single pass over each FASTQ -> per-target FASTAs
    h1 = {t: open(r1_paths[t], 'a') for t in targets}
    n1 = write_pass(a.r1, read_to_targets, h1)
    for f in h1.values():
        f.close()
    h2 = {t: open(r2_paths[t], 'a') for t in targets}
    n2 = write_pass(a.r2, read_to_targets, h2)
    for f in h2.values():
        f.close()
    print(f"[{a.sample}] wrote R1 reads={n1:,}  R2 reads={n2:,} across {len(targets)} taxids")


if __name__ == "__main__":
    main()
