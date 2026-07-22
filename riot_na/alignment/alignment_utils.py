from collections.abc import Iterator
from dataclasses import replace
from itertools import groupby

from riot_na.data.model import AlignmentEntry, AlignmentString, Cigar

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
# Degenerate amino-acid codes when an ambiguous codon encodes exactly these sets.
_DEGENERATE_AA = {
    frozenset({"D", "N"}): "B",
    frozenset({"E", "Q"}): "Z",
    frozenset({"I", "L"}): "J",
}
_GAP = "."


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
    if not has_stop and len(amino_acids) == 1:
        return next(iter(amino_acids))
    if not has_stop:
        degenerate = _DEGENERATE_AA.get(frozenset(amino_acids))
        if degenerate is not None:
            return degenerate
    return "X"


def translate(query_sequence: str, coding_frame: int) -> str:
    """Translate DNA in ``coding_frame`` with the standard genetic code.

    Behaviour matches common translate defaults: stop → ``*``, unresolved /
    mixed ambiguous codons → ``X``, and gap codon ``...`` → ``.``.
    """
    assert coding_frame in [0, 1, 2]

    query_sequence = query_sequence[coding_frame:].upper()
    n = len(query_sequence) - len(query_sequence) % 3
    return "".join(_translate_codon(query_sequence[i : i + 3]) for i in range(0, n, 3))


def get_cigar_op_groups(cigar: Cigar) -> Iterator[tuple[int, str]]:
    cigar_groups = (list(grouper) for _, grouper in groupby(cigar, lambda character: character.isdigit()))
    for group_size_grouper, op_grouper in zip(cigar_groups, cigar_groups):
        yield int("".join(group_size_grouper)), op_grouper[0]


def unfold_cigar(cigar: Cigar) -> AlignmentString:
    result = []
    for group_size, operation in get_cigar_op_groups(cigar):
        for _ in range(group_size):
            result.append(operation)
    return AlignmentString("".join(result))


def fold_cigar(alignment_str: AlignmentString) -> Cigar:
    return Cigar("".join(f"{sum(1 for _ in op_group)}{op}" for op, op_group in groupby(alignment_str)))


def has_frameshift(cigar: Cigar) -> bool:
    counter = 0

    for size, operation in get_cigar_op_groups(cigar):
        if operation == "M":
            if counter % 3 != 0:
                return True
            counter = 0

        if operation == "I":
            counter -= size
        elif operation == "D":
            counter += size

    return False


def infer_reading_frame(t_start: int, t_frame: int) -> int:
    #             AT GAC TATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCAG
    #                 MM MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    #  T TAC TAC TAC TAC TACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG target reading frame is 1
    #    |            |
    #  t_frame     t_start
    # t_start = 11
    # t_frame = 1
    # frame_offset = ((t_start -t_frame) %3) = 10 % 3 = 1
    # q_frame = (3 - frame_offset) %3 = (3 - 1) % 3 = 2

    # offset from the reading frame of the first aligned position
    frame_offset = (t_start - t_frame) % 3
    q_frame = (3 - frame_offset) % 3
    return q_frame


def offset_alignments(offset: int, aln: AlignmentEntry) -> AlignmentEntry:
    return replace(aln, q_start=aln.q_start + offset, q_end=aln.q_end + offset)


def align_sequences(query: str, target: str, alignment: AlignmentString) -> tuple[str, str]:
    query_aln = []
    target_aln = []

    for op in alignment:
        if op == "M":
            query_aln.append(query[0])
            target_aln.append(target[0])
            query = query[1:]
            target = target[1:]
        elif op == "I":
            query_aln.append(query[0])
            target_aln.append("-")
            query = query[1:]
        elif op == "D":
            query_aln.append("-")
            target_aln.append(target[0])
            target = target[1:]

    return "".join(query_aln), "".join(target_aln)
