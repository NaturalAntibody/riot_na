"""Regression: our translate() must match BioPython Seq.translate(gap='.')."""

from itertools import product

import pytest
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from riot_na.alignment.alignment_utils import translate

# Matches the previous riot_na call: Seq(...).translate(gap=".")
_GAP = "."


def _biopython_translate(sequence: str, coding_frame: int = 0) -> str:
    sequence = sequence[coding_frame:]
    sequence = sequence[: len(sequence) - len(sequence) % 3]
    return str(Seq(sequence).translate(gap=_GAP))


@pytest.mark.parametrize(
    ("sequence", "coding_frame"),
    [
        ("ATG", 0),
        ("ATGTAA", 0),
        ("ATG...", 0),
        ("...", 0),
        ("NNN", 0),
        ("ATGNNN", 0),
        ("ACN", 0),
        ("CTN", 0),
        ("TAR", 0),
        ("TRA", 0),
        ("TAN", 0),
        ("AAR", 0),
        ("AAY", 0),
        ("RAT", 0),
        ("WTA", 0),
        ("atg", 0),
        ("XATGCCCTAG", 1),
        ("XXATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 2),
        ("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 0),
        ("ATGCC", 0),  # partial trailing codon dropped
    ],
)
def test_translate_matches_biopython(sequence: str, coding_frame: int):
    assert translate(sequence, coding_frame) == _biopython_translate(sequence, coding_frame)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGT", repeat=3)])
def test_translate_all_unambiguous_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGTRYWSMKHBVDN", repeat=3)])
def test_translate_all_iupac_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


def test_translate_rejects_partial_gap_codon():
    with pytest.raises(ValueError, match="invalid"):
        translate("AT.", 0)

    with pytest.raises(TranslationError):
        _biopython_translate("AT.")
