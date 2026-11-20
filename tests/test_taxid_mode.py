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
    (131567, 1, "no rank", "cellular organisms"),
    (2, 131567, "superkingdom", "Bacteria"),
    (2323, 2, "no rank", "unclassified Bacteria"),
    (1869227, 2323, "species", "bacterium"),
    (28384, 1, "no rank", "other sequences"),
    (29278, 28384, "no rank", "vectors"),
    (2911, 29278, "species", "Shuttle vector pTEX"),
    (10239, 1, "superkingdom", "Viruses"),
    (2731619, 10239, "class", "Caudoviricetes"),
    (1980, 2731619, "species", "Caudoviricetes sp."),
    (543, 2, "family", "Enterobacteriaceae"),
    (570, 543, "genus", "Klebsiella"),
    (3390273, 570, "species group", "Klebsiella pneumoniae complex"),
    (573, 3390273, "species", "Klebsiella pneumoniae"),
    (244366, 3390273, "species", "Klebsiella variicola"),
    (72407, 573, "subspecies", "Klebsiella pneumoniae subsp. pneumoniae"),
    (54291, 570, "species", "Klebsiella ornithinolytica"),
    (575, 570, "species", "Klebsiella planticola"),
    (561, 543, "genus", "Escherichia"),
    (562, 561, "species", "Escherichia coli"),
    (48479, 2, "no rank", "environmental samples"),
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
    """rows: (read, staxids, bitscore, pident, title[, length])"""
    rows = [r if len(r) == 6 else (*r, 150) for r in rows]
    df = pd.DataFrame(rows, columns=["read", "staxids", "bitscore", "pident", "stitle", "length"])
    df["qseqid"] = "T|R1|" + df["read"]
    df["hit_organism"] = df["stitle"]
    return df


def _labels(rows, target, margin=0.02, id_min=80.0, min_len=0, noncell=False, plasmid=False, frac=0):
    return label_reads(per_read_frame_tax(_hits(rows), target, TAX, ignore_noncellular=noncell,
                                          ignore_plasmid=plasmid), margin, id_min, min_len, frac).to_dict()


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


def test_min_alignment_length_16():
    rows = [("r", "573", 95, 100.0, "K. pneumoniae", 51)]
    assert _labels(rows, 573)["T|R1|r"] == "CONFIRMED", "no length rule by default"
    assert _labels(rows, 573, min_len=100)["T|R1|r"] == "EXCLUDED_SHORT"
    assert _labels([("s", "573", 270, 100.0, "K. pneumoniae", 150)], 573, min_len=100)["T|R1|s"] == "CONFIRMED"
    # a short alignment never rescues a read the target did not win: still MISASSIGNED
    assert _labels([("t", "562", 279, 100.0, "E. coli", 150), ("t", "573", 95, 100.0, "K. pn", 51)], 573,
                   min_len=100)["T|R1|t"] == "MISASSIGNED"


def test_noncellular_hits_ignored_17():
    assert TAX.is_noncellular(2911) and TAX.is_noncellular(1980) and TAX.is_noncellular(1869227)
    assert not TAX.is_noncellular(573)
    # vector / phage / "bacterium" outscore the target: MISASSIGNED unless #17 ignores them
    rows = [("r", "2911", 279, 100.0, "Shuttle vector pTEX"), ("r", "1980", 279, 100.0, "Caudoviricetes sp."),
            ("r", "1869227", 279, 100.0, "bacterium"), ("r", "573", 240, 99.0, "K. pneumoniae")]
    assert _labels(rows, 573)["T|R1|r"] == "MISASSIGNED"
    assert _labels(rows, 573, noncell=True)["T|R1|r"] == "CONFIRMED"
    # only non-cellular hits -> EXCLUDED (nothing named left), not MISASSIGNED
    assert _labels([("v", "2911", 279, 100.0, "Shuttle vector")], 573, noncell=True)["T|R1|v"] == "EXCLUDED"


def test_plasmid_hits_ignored_19():
    # best hit is a C.-difficile-filed plasmid: counts as target without #19, is dropped with it
    rows = [("r", "573", 182, 99.0, "Klebsiella pneumoniae strain X plasmid pKP1, complete sequence"),
            ("r", "562", 150, 99.0, "Escherichia coli chromosome")]
    assert _labels(rows, 573)["T|R1|r"] == "CONFIRMED"
    # read on a mobile element: neither for nor against the target
    assert _labels(rows, 573, plasmid=True)["T|R1|r"] == "EXCLUDED_PLASMID"
    assert _labels([("p", "573", 182, 99.0, "K. pneumoniae plasmid pKP1")], 573, plasmid=True)["T|R1|p"] == "EXCLUDED_PLASMID"
    # target chromosome scores with the plasmid: still CONFIRMED on the non-plasmid hit
    assert _labels(rows + [("r", "573", 181, 99.0, "Klebsiella pneumoniae chromosome")], 573,
                   plasmid=True)["T|R1|r"] == "CONFIRMED"
    # plasmid hit well below a non-target chromosome: an ordinary MISASSIGNED read
    rows2 = [("m", "562", 279, 100.0, "Escherichia coli chromosome"),
             ("m", "573", 200, 99.0, "Klebsiella pneumoniae plasmid pKP1")]
    assert _labels(rows2, 573, plasmid=True)["T|R1|m"] == "MISASSIGNED"
    assert summarize_labels(pd.Series(["CONFIRMED", "EXCLUDED_PLASMID"]))["rl_excluded_plasmid"] == 1


def test_identity_floor_90():
    rows = [("r", "573", 182, 87.8, "K. pneumoniae chromosome")]
    assert _labels(rows, 573, id_min=80)["T|R1|r"] == "CONFIRMED"
    assert _labels(rows, 573, id_min=90)["T|R1|r"] == "MISASSIGNED"


def test_target_fraction_rule_option():
    # #20: share of the target among kept named hits within 2% of the top; uncultured and lower-scoring hits do not count
    rows = [("r", "573", 279, 100.0, "K. pneumoniae")] + [("r", "562", 279, 100.0, "E. coli")] * 9 \
        + [("r", "77133", 279, 100.0, "uncultured bacterium")] * 20
    assert _labels(rows, 573)["T|R1|r"] == "CONFIRMED", "rule off by default"
    assert _labels(rows, 573, frac=0.10)["T|R1|r"] == "CONFIRMED", "1 of 10 at the top = 10% passes; uncultured not counted"
    assert _labels(rows + [("r", "562", 279, 100.0, "E. coli")], 573, frac=0.10)["T|R1|r"] == "EXCLUDED_NONSPECIFIC"
    # target alone at the top: lower-scoring hits filling the list never count
    solo = [("r", "573", 279, 100.0, "K. pneumoniae")] + [("r", "1280", 250, 95.0, "S. aureus")] * 20
    assert _labels(solo, 573, frac=0.50)["T|R1|r"] == "CONFIRMED", "sole best target = 100% of the top hits"
    # a near tie (within 2%) counts in the share: 1 target + 3 E. coli at 275/279 -> 25%
    near = solo + [("r", "562", 275, 99.6, "E. coli")] * 3
    assert _labels(near, 573, frac=0.25)["T|R1|r"] == "CONFIRMED"
    assert _labels(near, 573, frac=0.30)["T|R1|r"] == "EXCLUDED_NONSPECIFIC"
    # a hit just outside 2% (273/279 = 97.8%) does not count
    assert _labels(solo + [("r", "562", 273, 99.0, "E. coli")] * 3, 573, frac=0.50)["T|R1|r"] == "CONFIRMED"


def test_summary_fractions():
    lab = pd.Series(["CONFIRMED"] * 6 + ["MISASSIGNED"] * 2 + ["EXCLUDED_COMPLEX"] * 2)
    s = summarize_labels(lab)
    assert s["validated_fraction"] == 0.6 and s["validated_fraction_informative"] == 0.75
    s2 = summarize_labels(pd.Series(["CONFIRMED"] * 6 + ["MISASSIGNED"] * 2 + ["EXCLUDED_SHORT"] * 2))
    assert s2["rl_excluded_short"] == 2 and s2["validated_fraction"] == 0.6 and s2["validated_fraction_informative"] == 0.75
    assert "supported" not in s and s["rl_excluded_complex"] == 2


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
