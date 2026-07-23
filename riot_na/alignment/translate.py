"""Translate DNA/RNA sequences with the standard genetic code."""

# Standard genetic code (NCBI table 1).
_CODONS = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}
_IUPAC = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
}
# Letters accepted when a codon cannot be resolved to a single residue / stop
# (ambiguous DNA+RNA alphabets; U is normalized to T first).
_VALID_NUCLEOTIDES = frozenset("ACGTRYWSMKHBVDN")
# Degenerate amino-acid codes when an ambiguous codon encodes exactly these sets.
_DEGENERATE_AA = {
    frozenset({"D", "N"}): "B",
    frozenset({"E", "Q"}): "Z",
    frozenset({"I", "L"}): "J",
}
_GAP = "."


def _ambiguous_placeholder(codon: str) -> str:
    if _VALID_NUCLEOTIDES.issuperset(codon):
        return "X"
    raise ValueError(f"Codon {codon!r} is invalid")


def _translate_codon(codon: str) -> str:
    if codon == _GAP * 3:
        return _GAP
    if _GAP in codon:
        raise ValueError(f"Codon {codon!r} is invalid")

    known = _CODONS.get(codon)
    if known is not None:
        return known

    try:
        expansions = [a + b + c for a in _IUPAC[codon[0]] for b in _IUPAC[codon[1]] for c in _IUPAC[codon[2]]]
    except KeyError as exc:
        raise ValueError(f"Codon {codon!r} is invalid") from exc

    amino_acids: set[str] = set()
    has_stop = False
    for exp in expansions:
        aa = _CODONS[exp]
        if aa == "*":
            has_stop = True
        else:
            amino_acids.add(aa)

    if has_stop and not amino_acids:
        return "*"
    if has_stop:
        # Mixed stop + amino acid (e.g. TAN, NNN): X only for valid IUPAC letters.
        return _ambiguous_placeholder(codon)
    if len(amino_acids) == 1:
        return next(iter(amino_acids))
    degenerate = _DEGENERATE_AA.get(frozenset(amino_acids))
    if degenerate is not None:
        return degenerate
    # Multiple amino acids with no tighter degenerate code (e.g. ATX → I/M) → X.
    return "X"


def translate(query_sequence: str, coding_frame: int) -> str:
    """Translate DNA/RNA in ``coding_frame`` with the standard genetic code.

    Accepts DNA (``T``) and RNA (``U``), including mixed sequences. Uracil is
    treated as thymine. Stop → ``*``, unresolved / mixed ambiguous codons →
    ``X``, and gap codon ``...`` → ``.``.
    """
    assert coding_frame in [0, 1, 2]

    # Standard table is DNA/RNA-generic; U is equivalent to T.
    query_sequence = query_sequence[coding_frame:].upper().replace("U", "T")
    n = len(query_sequence) - len(query_sequence) % 3
    return "".join(_translate_codon(query_sequence[i : i + 3]) for i in range(0, n, 3))
