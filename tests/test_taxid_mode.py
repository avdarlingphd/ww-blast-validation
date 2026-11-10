#!/usr/bin/env python3
"""
Unit tests for v4 taxid-mode read validation (per_read_frame_tax / label_reads / summarize_labels)
and the Taxonomy helper. Uses a tiny synthetic NCBI-style taxonomy, so no taxdump is needed.

Locks down:
  - renamed genus (Raoultella -> Klebsiella, same taxid) CONFIRMS, whatever the title says,
  - "MAG:" titles and strain-level / merged taxids still match the target species,
  - same-complex sibling winning -> EXCLUDED_COMPLEX; wider margin can turn it CONFIRMED,
  - other-genus winner -> MISASSIGNED; tie between a sibling and a non-complex organism -> MISASSIGNED,
  - uncultured-only hits -> EXCLUDED.

    python tests/test_taxid_mode.py      (or pytest tests/)
"""
import os, sys, tempfile
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "workflow", "scripts"))
from taxonomy import Taxonomy                                                    # noqa: E402
from false_positive_detection import per_read_frame_tax, label_reads, summarize_labels  # noqa: E402

NODES = [  # taxid, parent, rank, name
    (1, 1, "no rank", "root"),
    (543, 1, "family", "Enterobacteriaceae"),
    (570, 543, "genus", "Klebsiella"),
    (3390273, 570, "species group", "Klebsiella pneumoniae complex"),
    (573, 3390273, "species", "Klebsiella pneumoniae"),
    (244366, 3390273, "species", "Klebsiella variicola"),
    (72407, 573, "subspecies", "Klebsiella pneumoniae subsp. pneumoniae"),
    (54291, 570, "species", "Klebsiella ornithinolytica"),
    (575, 570, "species", "Klebsiella planticola"),
    (561, 543, "genus", "Escherichia"),
    (562, 561, "species", "Escherichia coli"),
    (48479, 1, "no rank", "environmental samples"),
    (77133, 48479, "species", "uncultured bacterium"),
]
MERGED = [(99999, 54291)]  # an old taxid now merged into 54291


def _tax():
    f = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    f.write("taxid\tparent\trank\tname\n")
    for t, p, r, n in NODES:
        f.write(f"{t}\t{p}\t{r}\t{n}\n")
    for o, n in MERGED:
        f.write(f"{o}\t{n}\tMERGED\t\n")
    f.close()
    return Taxonomy(f.name)


TAX = _tax()


def _hits(rows):
    """rows: (read, staxids, bitscore, pident, title)"""
    df = pd.DataFrame(rows, columns=["read", "staxids", "bitscore", "pident", "stitle"])
    df["qseqid"] = "T|R1|" + df["read"]
    df["hit_organism"] = df["stitle"]
    return df


def _labels(rows, target, margin=0.02, id_min=80.0):
    return label_reads(per_read_frame_tax(_hits(rows), target, TAX), margin, id_min).to_dict()


def test_taxonomy_basics():
    assert TAX.target_key(72407) == 573, "subspecies target rolls up to species"
    assert TAX.target_key(570) == 570, "genus target stays a genus"
    assert TAX.complex_of(573) == 3390273 and TAX.complex_of(244366) == 3390273
    assert TAX.complex_of(54291) is None
    assert TAX.resolve(99999) == 54291, "merged taxid resolves"
    assert TAX.is_uncultured(77133)


def test_renamed_genus_confirms():
    # The Raoultella bug: NCBI calls 54291 Klebsiella ornithinolytica; title says Klebsiella.
    lab = _labels([("r1", "54291", 279, 100.0, "Klebsiella ornithinolytica strain WM1 chromosome"),
                   ("r1", "54291", 274, 99.3, "Raoultella ornithinolytica strain X")], 54291)
    assert lab["T|R1|r1"] == "CONFIRMED", lab


def test_mag_title_strain_and_merged_taxids_confirm():
    lab = _labels([("a", "54291", 279, 100.0, "MAG: Raoultella ornithinolytica isolate 7"),
                   ("b", "99999", 279, 100.0, "old record"),
                   ("c", "72407", 279, 100.0, "K. pneumoniae subsp. pneumoniae")], 54291)
    assert lab["T|R1|a"] == "CONFIRMED" and lab["T|R1|b"] == "CONFIRMED", lab
    lab = _labels([("c", "72407", 279, 100.0, "K. pneumoniae subsp. pneumoniae")], 573)
    assert lab["T|R1|c"] == "CONFIRMED", "strain/subspecies hit counts for the species target"


def test_complex_sibling_excluded_then_confirmed_at_wider_margin():
    rows = [("r", "244366", 279, 100.0, "Klebsiella variicola"),
            ("r", "573", 260, 97.0, "Klebsiella pneumoniae")]       # ratio 0.932
    assert _labels(rows, 573, margin=0.02)["T|R1|r"] == "EXCLUDED_COMPLEX"
    assert _labels(rows, 573, margin=0.10)["T|R1|r"] == "CONFIRMED"
    # sibling wins and the target is absent entirely -> still EXCLUDED_COMPLEX
    assert _labels([("s", "244366", 279, 100.0, "K. variicola")], 573)["T|R1|s"] == "EXCLUDED_COMPLEX"


def test_other_genus_is_misassigned_and_mixed_tie_is_misassigned():
    assert _labels([("r", "562", 279, 100.0, "E. coli"),
                    ("r", "573", 200, 90.0, "K. pneumoniae")], 573)["T|R1|r"] == "MISASSIGNED"
    # top tie between a complex sibling and E. coli -> not ALL top hits are siblings -> MISASSIGNED
    assert _labels([("t", "244366", 279, 100.0, "K. variicola"),
                    ("t", "562", 279, 100.0, "E. coli")], 573)["T|R1|t"] == "MISASSIGNED"
    # Raoultella has no NCBI complex: a Klebsiella sibling winning is MISASSIGNED, not EXCLUDED_COMPLEX
    assert _labels([("u", "575", 279, 100.0, "K. planticola")], 54291)["T|R1|u"] == "MISASSIGNED"


def test_uncultured_only_is_excluded_and_ignored_for_top():
    lab = _labels([("x", "77133", 300, 100.0, "uncultured bacterium clone 5"),
                   ("y", "77133", 300, 100.0, "uncultured bacterium"),
                   ("y", "573", 279, 100.0, "K. pneumoniae")], 573)
    assert lab["T|R1|x"] == "EXCLUDED", lab
    assert lab["T|R1|y"] == "CONFIRMED", "uncultured top hit is ignored"


def test_low_identity_target_is_misassigned():
    assert _labels([("r", "573", 279, 79.0, "K. pneumoniae")], 573)["T|R1|r"] == "MISASSIGNED"


def test_summary_fractions():
    lab = pd.Series(["CONFIRMED"] * 6 + ["MISASSIGNED"] * 2 + ["EXCLUDED_COMPLEX"] * 2)
    s = summarize_labels(lab)
    assert s["validated_fraction"] == 0.6 and s["validated_fraction_informative"] == 0.75
    assert s["supported"] and s["rl_excluded_complex"] == 2


if __name__ == "__main__":
    tests = [v for k, v in dict(globals()).items() if k.startswith("test_")]
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
