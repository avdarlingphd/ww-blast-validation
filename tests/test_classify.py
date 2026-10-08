#!/usr/bin/env python3
"""
Unit tests for the BLAST false-positive classification logic
(workflow/scripts/false_positive_detection.py).

These lock down the thresholding DECISIONS independently of BLAST:
  - species-level organism matching (genus + epithet must agree),
  - the TRUE / FALSE / UNCERTAIN / UNCULTURED_DOMINANT labels,
  - the ">2 matching reads" requirement and the 2-read floor.

Run either way:
    pytest tests/                   # if pytest is installed
    python tests/test_classify.py   # no pytest needed (self-running)
"""
import importlib.util
import os
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.join(_HERE, "..", "workflow", "scripts", "false_positive_detection.py")
_spec = importlib.util.spec_from_file_location("fp", _SRC)
fp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fp)


def _set(n_reads, n_match, pct_match, rank, pct_uncultured=0.0):
    """One read-set summary, shaped like analyze_tsv's return value."""
    return {"status": "complete", "n_reads": n_reads, "n_match": n_match,
            "pct_match": pct_match, "expected_rank": rank,
            "pct_uncultured": pct_uncultured}


# ---- organism_matches: species-level matching ------------------------------
def test_exact_species_matches():
    assert fp.organism_matches("Klebsiella pneumoniae", "Klebsiella pneumoniae")

def test_subspecies_hit_matches_species():
    assert fp.organism_matches("Klebsiella pneumoniae subsp ozaenae", "Klebsiella pneumoniae")

def test_strain_suffix_matches_species():
    assert fp.organism_matches("Escherichia coli O157:H7", "Escherichia coli")

def test_sibling_species_does_not_match():
    assert not fp.organism_matches("Klebsiella variicola", "Klebsiella pneumoniae")

def test_genus_only_hit_does_not_confirm_species():
    # the v2.1 tightening: a bare genus hit must NOT confirm a species-level candidate
    assert not fp.organism_matches("Klebsiella", "Klebsiella pneumoniae")

def test_genus_level_candidate_matches_at_genus():
    assert fp.organism_matches("Acinetobacter baumannii", "Acinetobacter")

def test_wrong_genus_does_not_match():
    assert not fp.organism_matches("Pseudomonas aeruginosa", "Acinetobacter")


# ---- classify: label thresholds --------------------------------------------
def test_true_positive():
    label, *_ = fp.classify(_set(10, 10, 100.0, 1), _set(10, 10, 100.0, 1))
    assert label == "TRUE_POSITIVE"

def test_false_positive():
    label, *_ = fp.classify(_set(50, 0, 0.0, None), _set(50, 0, 0.0, None))
    assert label == "FALSE_POSITIVE"

def test_uncultured_dominant():
    label, *_ = fp.classify(_set(20, 18, 90.0, 1, pct_uncultured=80.0), None)
    assert label == "UNCULTURED_DOMINANT"

def test_mid_range_is_uncertain():
    # 50% match is between 10 and 80 -> UNCERTAIN
    label, *_ = fp.classify(_set(20, 10, 50.0, 1), None)
    assert label == "UNCERTAIN"


# ---- classify: >2 matching reads + 2-read floor ----------------------------
def test_two_matching_reads_is_uncertain():
    # high %match but only 2 matching reads -> not enough (min_match_reads=3)
    label, *_ = fp.classify(_set(2, 2, 100.0, 1), None)
    assert label == "UNCERTAIN"

def test_three_matching_reads_is_true():
    label, *_ = fp.classify(_set(3, 3, 100.0, 1), None)
    assert label == "TRUE_POSITIVE"

def test_single_read_is_uncertain():
    # a lone read is trivially 100% match but fails both the read floor and min_match_reads
    label, *_ = fp.classify(_set(1, 1, 100.0, 1), None)
    assert label == "UNCERTAIN"

def test_expected_rank_worse_than_two_is_uncertain():
    label, *_ = fp.classify(_set(10, 10, 100.0, 3), None)
    assert label == "UNCERTAIN"

def test_no_data_when_no_reads():
    label, *_ = fp.classify(None, None)
    assert label == "NO_DATA"


# ---- complex-aware scoring --------------------------------------------------
def _crosswalk(members_by_complex):
    """Build (member_to_complex, complex_members) like load_complex_crosswalk would, in memory."""
    m2c, members = {}, {}
    for cx, orgs in members_by_complex.items():
        for org in orgs:
            m2c[fp._name_key(org)] = cx
            members.setdefault(cx, []).append(fp._tokens(org))
    return m2c, members

def _reads_df(pairs):
    """A per-read-set DataFrame (one row per read = its top hit), shaped like load_sample output."""
    return pd.DataFrame([
        {"qseqid": q, "sscinames": "N/A", "pident": 100.0, "length": 150,
         "evalue": 0.0, "bitscore": 200.0, "stitle": org, "hit_organism": org}
        for q, org in pairs])

def test_candidate_complex_lookup():
    m2c, _ = _crosswalk({"ACB": ["Acinetobacter baumannii", "Acinetobacter nosocomialis"]})
    assert fp.candidate_complex("Acinetobacter baumannii", m2c) == "ACB"
    assert fp.candidate_complex("Acinetobacter nosocomialis", m2c) == "ACB"
    assert fp.candidate_complex("Escherichia coli", m2c) is None

def test_complex_matcher_matches_any_member():
    _, members = _crosswalk({"ACB": ["Acinetobacter baumannii", "Acinetobacter nosocomialis"]})
    mf = fp.make_complex_matcher(members["ACB"])
    assert mf("Acinetobacter baumannii")
    assert mf("Acinetobacter nosocomialis strain X")   # sibling member -> complex match
    assert not mf("Klebsiella pneumoniae")             # outside the complex

def test_complex_rescues_species_split():
    # 6 reads: 3 hit A. baumannii, 3 hit A. nosocomialis (same ACB complex)
    df = _reads_df([("r1", "Acinetobacter baumannii"), ("r2", "Acinetobacter baumannii"),
                    ("r3", "Acinetobacter baumannii"), ("r4", "Acinetobacter nosocomialis"),
                    ("r5", "Acinetobacter nosocomialis"), ("r6", "Acinetobacter nosocomialis")])
    _, members = _crosswalk({"ACB": ["Acinetobacter baumannii", "Acinetobacter nosocomialis"]})
    # species level: only 3/6 match A. baumannii -> UNCERTAIN
    sp = fp.analyze_reads(df, "Acinetobacter baumannii")
    assert sp["pct_match"] == 50.0
    assert fp.classify(sp, None)[0] == "UNCERTAIN"
    # complex level: all 6 match the ACB complex -> TRUE_POSITIVE
    mf = fp.make_complex_matcher(members["ACB"])
    cx = fp.analyze_reads(df, "Acinetobacter baumannii", match_fn=mf)
    assert cx["pct_match"] == 100.0
    assert fp.classify(cx, None)[0] == "TRUE_POSITIVE"


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)
