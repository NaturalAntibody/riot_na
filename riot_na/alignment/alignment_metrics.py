import math
from dataclasses import dataclass
from typing import Optional

from riot_na.alignment.alignment_utils import get_cigar_op_groups, unfold_cigar
from riot_na.data.constants import BLOSUM_62
from riot_na.data.model import Cigar


@dataclass(frozen=True)
class GumbelParams:
    """
    Lambda and K of the Gumbel distribution of local alignment scores. They depend on the scoring system (matrix or
    match / mismatch scores and gap penalties), see riot_na.alignment.gumbel for how they are obtained (ALP,
    Ascending Ladder Program - https://doi.org/10.1093%2Fbioinformatics%2Fbtv575).
    """

    lambda_: float
    K: float  # pylint: disable=invalid-name


def compute_raw_score_aa(
    query: str,
    target: str,
    cigar: Cigar,
    gap_open_penalty=11,
    gap_extend_penalty=1,
    substitution_matrix=BLOSUM_62,
) -> float:
    raw_score = 0
    query_pos = 0
    target_pos = 0
    prev_match = True
    for op in unfold_cigar(cigar):
        if op == "M":
            raw_score += substitution_matrix[query[query_pos]][target[target_pos]]

            query_pos += 1
            target_pos += 1
            prev_match = True
        elif op == "I":
            if prev_match:
                raw_score -= gap_open_penalty
                prev_match = False
            else:
                raw_score -= gap_extend_penalty
            query_pos += 1
        elif op == "D":
            if prev_match:
                raw_score -= gap_open_penalty
                prev_match = False
            else:
                raw_score -= gap_extend_penalty
            target_pos += 1
    return raw_score


def compute_bit_score(raw_score: int, gumbel_params: GumbelParams) -> float:
    return (gumbel_params.lambda_ * raw_score - math.log(gumbel_params.K)) / math.log(2)


def compute_raw_score_from_bit_score(bit_score: float, gumbel_params: GumbelParams) -> float:
    return (math.log(gumbel_params.K) + bit_score * math.log(2)) / gumbel_params.lambda_


def log10_evalue_threshold(e_value_threshold: Optional[float]) -> Optional[float]:
    if e_value_threshold is None:
        return None
    assert e_value_threshold > 0, "E-value threshold must be greater than 0"
    return math.log10(e_value_threshold)


def compute_log10_evalue(query_length: int, db_length: int, bit_score: float) -> float:
    # log10(E) = log10(query) + log10(db) - bit_score * log10(2)
    return math.log10(query_length) + math.log10(db_length) - (bit_score * math.log10(2))


def calculate_seq_identity(cigar: Cigar, query: str, target: str, query_start: int = 0, target_start: int = 0) -> float:
    cigar_items = get_cigar_op_groups(cigar)
    query_pos = query_start
    target_pos = target_start
    match_cnt = 0
    total = 0
    for op_cnt, op in cigar_items:
        match op:
            case "M":
                for i in range(op_cnt):
                    if query[query_pos + i] == target[target_pos + i]:
                        match_cnt += 1
                query_pos += op_cnt
                target_pos += op_cnt
            case "I":
                query_pos += op_cnt
            case "D":
                target_pos += op_cnt
        total += op_cnt
    return match_cnt / total
