#!/usr/bin/env python3
"""Synonym-aware name matching (organism_matches + load_synonyms, resources/taxon_synonyms.csv).

    python tests/test_synonyms.py      (or pytest tests/)
"""
import os, sys
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "workflow", "scripts"))
import false_positive_detection as fpd  # noqa: E402

SYN = os.path.join(_HERE, "..", "resources", "taxon_synonyms.csv")


def test_without_synonyms_rename_does_not_match():
    fpd.load_synonyms(None)
    assert not fpd.organism_matches("Klebsiella ornithinolytica", "Raoultella ornithinolytica")


def test_rename_matches_with_synonyms():
    assert fpd.load_synonyms(SYN) > 100
    assert fpd.organism_matches("Klebsiella ornithinolytica", "Raoultella ornithinolytica")
    assert fpd.organism_matches("Klebsiella planticola", "Raoultella planticola")
    assert fpd.organism_matches("[Mycobacterium] chelonae", "Mycobacteroides chelonae")


def test_title_prefixes_are_stripped():
    fpd.load_synonyms(SYN)
    assert fpd.organism_matches("MAG: Rothia mucilaginosa", "Rothia mucilaginosa")
    assert fpd.organism_matches("MAG TPA_asm: Delftia acidovorans", "Delftia acidovorans")
    assert not fpd.organism_matches("MAG: Naegleria", "Naegleria fowleri"), "genus-only stays a non-match"


def test_synonyms_do_not_leak_to_other_species():
    fpd.load_synonyms(SYN)
    assert not fpd.organism_matches("Klebsiella pneumoniae", "Raoultella ornithinolytica")
    assert not fpd.organism_matches("Klebsiella ornithinolytica", "Klebsiella pneumoniae")
    assert not fpd.organism_matches("Mycobacteroides abscessus", "Mycobacteroides chelonae")
    # NCBI synonym "Aeromonas hydrophila subsp. anaerogenes" of A. caviae must NOT make A. hydrophila match
    assert not fpd.organism_matches("Aeromonas hydrophila", "Aeromonas caviae")
    assert fpd.organism_matches("Aeromonas punctata", "Aeromonas caviae")
    assert fpd.organism_matches("Peptoclostridium difficile", "Clostridioides difficile")
    assert fpd.organism_matches("Enterobacter aerogenes", "Klebsiella aerogenes")
    assert fpd.organism_matches("Clostridium difficile", "Clostridioides difficile"), "vector names must not block"
    assert not fpd.organism_matches("Klebsiella aerogenes", "Klebsiella pneumoniae"), "historic synonym now its own species"
    assert not fpd.organism_matches("Mycobacterium tuberculosis", "Mycobacterium avium")


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t(); print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1; print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)
