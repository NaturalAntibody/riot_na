import math

from riot_na.alignment.alignment_metrics import (
    GumbellParams,
    calculate_seq_identity,
    compute_bit_score,
    compute_evalue,
    compute_raw_score_aa,
    compute_raw_score_from_bit_score,
)
from riot_na.data.constants import BLOSUM_62


def test_blosum_62_is_available():
    assert BLOSUM_62["A"]["A"] == 4.0
    assert BLOSUM_62["W"]["W"] == 11.0
    assert BLOSUM_62["A"]["R"] == -1.0


def test_calculate_seq_identity():
    assert calculate_seq_identity(cigar="4M", query="ACGT", target="ACGA") == 0.75


def test_calculate_seq_identity_with_offsets():
    assert calculate_seq_identity(cigar="4M", query="XXACGT", target="YYACGA", query_start=2, target_start=2) == 0.75


def test_calculate_seq_identity_with_indels():
    assert calculate_seq_identity(cigar="2M1I1D1M", query="ACGT", target="ACTT") == 0.6


def test_compute_raw_score_aa_for_match_and_gap():
    assert compute_raw_score_aa(query="AR", target="AR", cigar="2M") == 9.0
    assert compute_raw_score_aa(query="ARN", target="AR", cigar="2M1I") == -2.0


def test_compute_raw_score_aa_for_gap_extension_and_deletion():
    assert compute_raw_score_aa(query="ARND", target="AR", cigar="2M2I") == -3.0
    assert compute_raw_score_aa(query="AR", target="ARN", cigar="2M1D") == -2.0


def test_compute_bit_score_round_trip():
    raw_score = 42
    bit_score = compute_bit_score(raw_score=raw_score, gumbell_params=GumbellParams.IGBLAST)
    reconstructed = compute_raw_score_from_bit_score(bit_score=bit_score, gumbell_params=GumbellParams.IGBLAST)

    assert math.isclose(reconstructed, raw_score)


def test_compute_bit_score_round_trip_for_aa_params():
    raw_score = 87.5
    bit_score = compute_bit_score(raw_score=raw_score, gumbell_params=GumbellParams.AA)
    reconstructed = compute_raw_score_from_bit_score(bit_score=bit_score, gumbell_params=GumbellParams.AA)

    assert math.isclose(reconstructed, raw_score)


def test_compute_evalue():
    assert compute_evalue(query_length=10, db_length=100, bit_score=5) == 31.25


def test_compute_evalue_decreases_with_higher_bit_score():
    low_score = compute_evalue(query_length=100, db_length=1000, bit_score=10)
    high_score = compute_evalue(query_length=100, db_length=1000, bit_score=20)

    assert high_score < low_score
