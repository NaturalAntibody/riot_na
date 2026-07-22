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
        # Leading junk skipped by frame (X is not a valid codon base in BioPython)
        ("XATGCCCTAG", 1),
        ("XXATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 2),
        ("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG", 0),
        ("ATGCC", 0),  # partial trailing codon dropped
        # RNA / uracil (BioPython Standard table is DNA/RNA-generic)
        ("AUG", 0),
        ("UUU", 0),
        ("UAA", 0),
        ("UAG", 0),
        ("UGA", 0),
        ("UAR", 0),
        ("URA", 0),
        ("auguuu", 0),
        ("AUG...UUU", 0),
        ("AUGCCCUAG", 0),
        # Mixed T/U in one sequence
        ("AUT", 0),
        ("AUGTTT", 0),
        ("ATGUUU", 0),
        ("YUA", 0),
        ("WUA", 0),
    ],
)
def test_translate_matches_biopython(sequence: str, coding_frame: int):
    assert translate(sequence, coding_frame) == _biopython_translate(sequence, coding_frame)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGT", repeat=3)])
def test_translate_all_unambiguous_dna_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGU", repeat=3)])
def test_translate_all_unambiguous_rna_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGTU", repeat=3)])
def test_translate_all_mixed_tu_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGTRYWSMKHBVDN", repeat=3)])
def test_translate_all_iupac_dna_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize("codon", ["".join(bases) for bases in product("ACGURYWSKMBDHVN", repeat=3)])
def test_translate_all_iupac_rna_codons(codon: str):
    assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize(
    ("codon", "expect_error"),
    [
        ("XXX", True),
        ("TAX", True),
        ("XAA", True),
        ("ATX", False),  # resolves to amino-acid X (I/M), not an error
        ("GGX", False),  # all expansions are G
    ],
)
def test_translate_nucleotide_x_matches_biopython(codon: str, expect_error: bool):
    if expect_error:
        with pytest.raises(ValueError, match="invalid"):
            translate(codon, 0)
        with pytest.raises(TranslationError):
            _biopython_translate(codon)
    else:
        assert translate(codon, 0) == _biopython_translate(codon)


@pytest.mark.parametrize(
    ("rna", "dna", "expected"),
    [
        ("AUG", "ATG", "M"),
        ("UUU", "TTT", "F"),
        ("UUC", "TTC", "F"),
        ("UUA", "TTA", "L"),
        ("UUG", "TTG", "L"),
        ("UCU", "TCT", "S"),
        ("UAU", "TAT", "Y"),
        ("UGU", "TGT", "C"),
        ("UGG", "TGG", "W"),
        ("UAA", "TAA", "*"),
        ("UAG", "TAG", "*"),
        ("UGA", "TGA", "*"),
        ("UAR", "TAR", "*"),
        ("URA", "TRA", "*"),
        ("AUGCCCUAG", "ATGCCCTAG", "MP*"),
        ("auguuu", "atgttt", "MF"),
    ],
)
def test_translate_rna_matches_dna_and_expected(rna: str, dna: str, expected: str):
    assert translate(rna, 0) == expected
    assert translate(dna, 0) == expected
    assert translate(rna, 0) == _biopython_translate(rna)
    assert translate(dna, 0) == _biopython_translate(dna)


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("AUT", "I"),  # A + U/T mix → Ile
        ("AUGTTT", "MF"),
        ("ATGUUU", "MF"),
        ("TTU", "F"),
        ("UUT", "F"),
        ("AUG...UUU", "M.F"),
        ("YUA", "L"),
        ("WUA", "J"),  # A/U + A → Ile/Leu
    ],
)
def test_translate_mixed_dna_rna_matches_expected(sequence: str, expected: str):
    assert translate(sequence, 0) == expected
    assert translate(sequence, 0) == _biopython_translate(sequence)


@pytest.mark.parametrize(
    ("codon", "expected"),
    [
        # B = Asp/Asn (D/N)
        ("RAC", "B"),
        ("RAT", "B"),
        ("RAY", "B"),
        ("RAU", "B"),  # RNA
        # Z = Glu/Gln (E/Q)
        ("SAA", "Z"),
        ("SAG", "Z"),
        ("SAR", "Z"),
        # J = Ile/Leu (I/L)
        ("WTA", "J"),
        ("WUA", "J"),  # RNA
        ("MTA", "J"),
        ("MTC", "J"),
        ("MTT", "J"),
        ("MTU", "J"),  # RNA
        ("MTY", "J"),
        ("MTW", "J"),
        ("MTM", "J"),
        ("MTH", "J"),
        ("HTA", "J"),
        ("MUA", "J"),  # RNA
    ],
)
def test_translate_degenerate_amino_acids(codon: str, expected: str):
    assert translate(codon, 0) == expected
    assert _biopython_translate(codon) == expected


def test_translate_rejects_partial_gap_codon():
    with pytest.raises(ValueError, match="invalid"):
        translate("AT.", 0)

    with pytest.raises(TranslationError):
        _biopython_translate("AT.")
