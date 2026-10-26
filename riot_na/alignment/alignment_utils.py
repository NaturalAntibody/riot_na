from collections.abc import Iterator, Sequence
from dataclasses import replace
from itertools import groupby
from typing import Optional, Protocol

from riot_na.data.model import AlignmentEntry, AlignmentString, Cigar, Locus


class QueryAlignment(Protocol):
    @property
    def q_end(self) -> int: ...


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


def v_domain_end(v_alignment: QueryAlignment, j_alignment: Optional[QueryAlignment], j_query_start: int = 0) -> int:
    """
    End of a V domain in the query of its V alignment: end of its J alignment or, if J was not found, of its V
    alignment. C region follows it. j_query_start is the start of the query of the J alignment, if J was aligned
    to the query left after the V alignment.
    """
    return j_query_start + j_alignment.q_end if j_alignment else v_alignment.q_end


def select_c_alignment(
    c_alignments: Sequence[Optional[AlignmentEntry]], v_locus: Optional[Locus]
) -> Optional[AlignmentEntry]:
    """
    Select the best of C alignments to C genes of different loci. C region can follow V region of another locus
    in engineered antibodies (e.g. VL-CH1 in CrossMabs), so C genes of all loci are aligned.
    Alignments are compared by score, as e-values depend on the size of the gene database of each locus.
    C genes of the V locus win ties.
    """
    return max(
        (alignment for alignment in c_alignments if alignment is not None),
        key=lambda alignment: (alignment.alignment_score or 0, alignment.locus == v_locus),
        default=None,
    )


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
