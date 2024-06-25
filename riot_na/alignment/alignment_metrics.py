import math
from enum import Enum

import blosum  # type: ignore

from riot_na.alignment.alignment_utils import get_cigar_op_groups, unfold_cigar
from riot_na.data.model import Cigar

# L (Lambda) and K constants are depended on scoring matrix and gap penalties
# and are calculated by ALP (Ascending Ladder Program) - https://doi.org/10.1093%2Fbioinformatics%2Fbtv575
# The present values were calculated for the following parameters (default SSW):
# Match score: 2
# Mismatch score:  -3
# Gap open penalty: 5
# Gap extend penalty: 2
# Background probabilities:  A, T, G, C 0.2499975, N 0.00001


class GumbellParams(Enum):
    IGBLAST = {"L": 1.08, "K": 0.28}
    AA = {"L": 0.26453605633241922, "K": 0.043186874595437463}


def compute_raw_score_aa(
    query: str,
    target: str,
    cigar: Cigar,
    gap_open_penalty=11,
    gap_extend_penalty=1,
    substitution_matrix=blosum.BLOSUM(62),
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


def compute_bit_score(raw_score: int, gumbell_params: GumbellParams = GumbellParams.IGBLAST) -> float:
    return (gumbell_params.value["L"] * raw_score - math.log(gumbell_params.value["K"])) / math.log(2)


def compute_raw_score_from_bit_score(bit_score: float, gumbell_params: GumbellParams = GumbellParams.IGBLAST) -> float:
    return (math.log(gumbell_params.value["K"]) + bit_score * math.log(2)) / gumbell_params.value["L"]


def compute_evalue(query_length: int, db_length: int, bit_score: float) -> float:
    return query_length * db_length * 2 ** (-bit_score)


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
