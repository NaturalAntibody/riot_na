from riot_na.alignment.alignment_utils import fold_cigar, get_cigar_op_groups, unfold_cigar


def test_get_cigar_op_groups_parses_multi_digit_groups():
    assert list(get_cigar_op_groups(cigar="10M2I3D")) == [(10, "M"), (2, "I"), (3, "D")]


def test_unfold_cigar_expands_groups():
    assert unfold_cigar(cigar="3M1I2D") == "MMMIDD"


def test_fold_cigar_collapses_alignment_string():
    assert fold_cigar(alignment_str="MMMIDD") == "3M1I2D"


def test_fold_and_unfold_cigar_round_trip():
    cigar = "12M1I4M2D"
    assert fold_cigar(alignment_str=unfold_cigar(cigar=cigar)) == cigar


def test_fold_cigar_handles_empty_alignment_string():
    assert fold_cigar(alignment_str="") == ""
