# pylint: skip-file
# flake8: noqa
import blosum
import pytest
from skbio.alignment import StripedSmithWaterman

from riot_na.alignment.skbio_alignment import align_aa
from riot_na.data.model import AlignmentEntryAA, AlignmentString, Cigar
from riot_na.schemes.scheme_alignment import (
    _merge_cigars,
    fix_imgt_cdrs_numbering,
    force_c_terminus_matches,
    force_n_terminus_matches,
    infer_last_scheme_position_aligned,
    merge_v_j,
)
from riot_na.schemes.smooth_alignment import (
    reorder_buffer_imgt,
    reorder_buffer_other_schemes,
)

AA_ALIGNER_PARAMS = {
    "gap_open_penalty": 11,
    "gap_extend_penalty": 1,
    "protein": True,
    "substitution_matrix": blosum.BLOSUM(62),
}


def test_align_aa():
    # given
    query = "AAAA"
    target = "BAAABBB"
    expected_cigar = Cigar("3M")

    # when
    aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)
    result = align_aa(aligner, query, "target_id", target, db_length=None, calculate_score=False)

    # then
    assert result.cigar == expected_cigar
    assert result.q_start == 0
    assert result.q_end == 3
    assert result.t_start == 1
    assert result.t_end == 4


def test_force_n_terminus_alignments():
    # given alignment of query and target
    alignment = AlignmentEntryAA(
        target_id="target_id",
        alignment_score=None,
        seq_identity=None,
        e_value=None,
        q_start=3,
        q_end=5,
        t_start=2,
        t_end=4,
        cigar=Cigar("3M"),
        species=None,
        locus=None,
        q_seq=None,
        t_seq=None,
    )

    # when
    result = force_n_terminus_matches(alignment)

    # then
    # verify 1S was appended to cigra (q_start - t_start == 1)
    # verify that no N was appended
    # verify that min(q_start, t_start) matches was appended to the cigar
    assert result.cigar == Cigar("2M3M")
    assert result.q_start == 1


def test_fix_c_terminus_alignment():
    # given alignment of query and target
    alignment = AlignmentEntryAA(
        target_id="target_id",
        alignment_score=None,
        seq_identity=None,
        e_value=None,
        q_start=1,
        q_end=4,
        t_start=1,
        t_end=4,
        cigar="3M",
        species=None,
        locus=None,
        q_seq=None,
        t_seq=None,
    )
    query_sequence = "BAAAA"
    target_sequence = "CAAABBB"
    # alignment:
    # q:    B-AAAA
    # t:    -CAAABBB
    # aln:  SNMMM...

    # when
    result = force_c_terminus_matches(query_sequence, target_sequence, alignment)

    # then
    # verify N nor S was not appended to cigar
    # verify that min(q_start, t_start) matches was appended to the cigar
    assert result.cigar == Cigar("3M1M")


def test_merge_cigars():
    # I M -> I: it1++ -ok
    # I D -> M:it1++, it2++ -ok
    # I I -> I: it2++ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # M M ->M:it1++, it2++ -ok
    # M D -> D: it2++ -ok
    # M I -> I: it1++, it2++ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # D D -> D: it1++ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # D M -> : D : it1++, it2++
    # D I -> X: it1++, it2++

    # M M
    assert _merge_cigars("MMMM", "MMMM") == "MMMM"
    # M D
    assert _merge_cigars("MMM", "MMDM") == "MMDM"
    assert _merge_cigars("MMM", "MMDDM") == "MMDDM"
    # M I
    assert _merge_cigars("MMMM", "MMIM") == "MMIM"
    assert _merge_cigars("MMMM", "MMIIM") == "MMII"

    # I M
    assert _merge_cigars("MMIM", "MMMM") == "MMIM"
    assert _merge_cigars("MMIIM", "MMMMD") == "MMIIM"
    # I D
    assert _merge_cigars("MMIM", "MMDM") == "MMMM"
    assert _merge_cigars("MMIIM", "MMDDM") == "MMMMM"
    # I I
    assert _merge_cigars("MMIM", "MMIM") == "MMII"
    assert _merge_cigars("MMIIM", "MMIM") == "MMIII"
    assert _merge_cigars("MMIM", "MMIIM") == "MMII"
    assert _merge_cigars("MMIIM", "MMIIM") == "MMIII"

    # D M
    assert _merge_cigars("MMDM", "MMMM") == "MMDM"
    assert _merge_cigars("MMDDM", "MMMMM") == "MMDDM"
    # D D
    with pytest.raises(AssertionError, match="Query gene alignment is longer than gene scheme alignment"):
        _merge_cigars("MMDM", "MMDM") == "MMDDM"

    assert _merge_cigars("MMDM", "MMDDMM") == "MMDDDM"
    assert _merge_cigars("MMDM", "MMDMDM") == "MMDDDM"


def test_merge_long_cigars():
    q_g = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMMMMM"
    # QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYIS SSGSTIYYADSVKGRFTISRDNAKNSLYLQM   NSLRAEDTAVYYCAR
    g_sch = "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIMMMMMMMMMMMM"

    #        NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMMMMMMMMMMMIIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMIIIIIIMMMMMMMMMMMMMMMMMMMIIIIIIIIMMMMMMMMMMMM

    res = _merge_cigars(q_g, g_sch)
    print(res)


def test_infer_last_scheme_position_aligned():
    # given
    query = "CAAAAAAAAAAAA"
    gene = "BAAAAAAABC"
    aln = "SNMMMMMMM"
    g_sch = "MMMDMMMMMM"
    q_sch = "SNMMDMMMMM"

    # query =  "C-AA-AAAAAAAAAA"
    # gene =   "-BAA-AAAAABC"
    # aln =    "SNMM-MMMMM"
    # g_sch =   "MMMDMMMMMM"
    # q_sch2 = "IDMMDMMMMM"
    # last_pos="         |"  # 9 (0-based)

    # when
    last_pos = infer_last_scheme_position_aligned(q_sch)

    # then
    assert last_pos == 9


def test_merge_v_j():
    # given
    query_imgt_cigar_v = AlignmentString(
        "MMMMMMMMMDMMMMMMMMMMMMMMMMMMMMDDDDMMMMMMMMMMMMMMMMMMMMMMMMMDDDMMMMMMMMMMDMMMMMMMMMMMMMMMMMMDMMMMMMMMMMMMMM"
    )  # pylint: disable=line-too-long
    unaligned_middle_len = 18
    query_imgt_cigar_j = AlignmentString(unaligned_middle_len * "S" + "DMMMMMMMMMMM")
    j_gene_start_on_scheme = 117
    # when
    result = merge_v_j(query_imgt_cigar_v, query_imgt_cigar_j, j_gene_start_on_scheme)

    # then
    assert result == AlignmentString(
        "MMMMMMMMMDMMMMMMMMMMMMMMMMMMMMDDDDMMMMMMMMMMMMMMMMMMMMMMMMMDDDMMMMMMMMMMDMMMMMMMMMMMMMMMMMMDMMMMMMMMMMMMMMMMMMMMIIIIIIIMMMMMMMMMMMMMMMM"
    )


def test_merge_v_j_where_j_starts_on_deletion():
    # given
    # query:
    # genes
    #
    # gene_imgt

    # query_imgt

    query = "AAAAAAAAAABBB DDDDDDDDB"  # query is masked
    gene = "CCC DDDDDDDD"
    aln = "             NNN MMMMMMMM"
    g_sch = "MMMDMMMMMMMM"  # NNNNNNNNNNNNNNN -15N -> gene_start = 15
    q_sch = "DDDDMMMMMM"
    gene_start_on_scheme = 15

    query_imgt_cigar_v = AlignmentString(
        "MMMMMMMMMDMMMMMMMMMMMMMMMMMMMMDDDDMMMMMMMMMMMMMMMMMMMMMMMMMDDDMMMMMMMMMMDMMMMMMMMMMMMMMMMMMDMMMMMMMMMMMMMM"
    )  # pylint: disable=line-too-long
    unaligned_middle_len = 18
    query_imgt_cigar_j = AlignmentString(unaligned_middle_len * "S" + "DMMMMMMMMMMMMMM")
    j_gene_start_on_scheme = 114

    # when
    result = merge_v_j(query_imgt_cigar_v, query_imgt_cigar_j, j_gene_start_on_scheme)

    # then
    assert result == AlignmentString(
        "MMMMMMMMMDMMMMMMMMMMMMMMMMMMMMDDDDMMMMMMMMMMMMMMMMMMMMMMMMMDDDMMMMMMMMMMDMMMMMMMMMMMMMMMMMMDMMMMMMMMMMMMMMMMMMIIIIIIIIIIMMMMMMMMMMMMMMMMMM"
    )


def test_reorder_buffer_matches():
    # given buffer = cdr + 2 neigbour residues on both sides
    buffer = AlignmentString("MMMMMMMM")

    # when
    result = reorder_buffer_other_schemes(buffer, relative_indel_position=5)

    # then
    assert result == buffer


def test_reorder_buffer_insertions():
    # given buffer = cdr + 2 neigbour residues on both sides
    buffer = AlignmentString("MIMMMMMMMI")

    # when
    result = reorder_buffer_other_schemes(buffer, relative_indel_position=5)

    # then
    assert result == AlignmentString("MMMMMIIMMM")


def test_reorder_buffer_deletions():
    # given buffer = cdr + 2 neigbour residues on both sides
    buffer = AlignmentString("MDMMMMMMDM")

    # when
    result = reorder_buffer_other_schemes(buffer, relative_indel_position=5)

    # then
    assert result == AlignmentString("MMMDDMMMMM")


def test_reorder_buffer_imgt_matches():
    # given even len(buffer) - cdr1 and cdr2
    buffer = AlignmentString("MMMMMMMM")
    assert len(buffer) % 2 == 0

    # when
    result = reorder_buffer_imgt(buffer, relative_indel_position=5)

    # then
    assert result == buffer

    # given odd len(buffer) - cdr3
    buffer = AlignmentString("MMMMMMMMM")
    assert len(buffer) % 2 == 1

    # when
    result = reorder_buffer_imgt(buffer, relative_indel_position=5)

    # then
    assert result == buffer


def test_reorder_buffer_imgt_insertions():
    # given even len(buffer) - cdr1 and cdr2
    buffer = AlignmentString("MIMMMMMMMI")
    assert len(buffer) % 2 == 0

    # when
    result = reorder_buffer_imgt(buffer, relative_indel_position=5)

    # then
    assert result == AlignmentString("MMMMMIIMMM")

    # given odd len(buffer) - cdr3
    buffer = AlignmentString("MIMIIIMMMMMMI")
    assert len(buffer) % 2 == 1

    # when
    result = reorder_buffer_imgt(buffer, relative_indel_position=5)

    # then
    assert result == AlignmentString("MMMMMIIIIIMMM")


def test_reorder_buffer_imgt_deletions():
    # given even len(buffer) - cdr1 and cdr2
    buffer = AlignmentString("MDMMMMMMDM")
    assert len(buffer) % 2 == 0

    # when
    result = reorder_buffer_imgt(buffer, relative_indel_position=5)

    # then
    assert result == AlignmentString("MMMMDDMMMM")

    # given odd len(buffer) - cdr3
    buffer = AlignmentString("MDMMMMMMDMM")
    assert len(buffer) % 2 == 1

    # when
    result = reorder_buffer_imgt(buffer, relative_indel_position=5)

    # then
    assert result == AlignmentString("MMMMDDMMMMM")


def test_fix_imgt_cdrs_numbering():
    # CDR1
    # No change
    numbering = {32: [("A", 0), ("B", 1)], 33: [("D", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {32: [("A", 0), ("B", 1)], 33: [("D", 0)]}

    # Odd insertion number
    numbering = {32: [("A", 0), ("B", 1), ("C", 2), ("D", 3)], 33: [("E", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {32: [("A", 0), ("B", 1), ("C", 2)], 33: [("E", 0), ("D", 1)]}

    # Even insertion number
    numbering = {32: [("A", 0), ("B", 1), ("C", 2), ("D", 3), ("E", 4)], 33: [("F", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {
        32: [("A", 0), ("B", 1), ("C", 2)],
        33: [("F", 0), ("E", 1), ("D", 2)],
    }

    # CDR2
    # No change
    numbering = {60: [("A", 0), ("B", 1)], 61: [("D", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {60: [("A", 0), ("B", 1)], 61: [("D", 0)]}

    # Odd insertion number
    numbering = {60: [("A", 0), ("B", 1), ("C", 2), ("D", 3)], 61: [("E", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {60: [("A", 0), ("B", 1), ("C", 2)], 61: [("E", 0), ("D", 1)]}

    # Even insertion number
    numbering = {60: [("A", 0), ("B", 1), ("C", 2), ("D", 3), ("E", 4)], 61: [("F", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {
        60: [("A", 0), ("B", 1), ("C", 2)],
        61: [("F", 0), ("E", 1), ("D", 2)],
    }

    # CDR3 - swap priority
    # No change
    numbering = {111: [("A", 0)], 112: [("B", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {111: [("A", 0)], 112: [("B", 0)]}

    # Odd insertion number
    numbering = {111: [("A", 0), ("B", 1)], 112: [("C", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {111: [("A", 0)], 112: [("C", 0), ("B", 1)]}

    # Even insertion number
    numbering = {111: [("A", 0), ("B", 1), ("C", 2)], 112: [("D", 0)]}
    assert fix_imgt_cdrs_numbering(numbering) == {111: [("A", 0), ("B", 1)], 112: [("D", 0), ("C", 1)]}


if __name__ == "__main__":
    test_align_aa()
    test_force_n_terminus_alignments()
    test_fix_c_terminus_alignment()
    test_merge_cigars()

    test_infer_last_scheme_position_aligned()

    test_merge_v_j()
    test_merge_v_j_where_j_starts_on_deletion()

    test_reorder_buffer_matches()
    test_reorder_buffer_insertions()
    test_reorder_buffer_deletions()

    test_merge_long_cigars()
