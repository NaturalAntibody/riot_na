from typing import Optional

from skbio.alignment import StripedSmithWaterman  # type: ignore

from riot_na.alignment.alignment_metrics import (
    GumbellParams,
    calculate_seq_identity,
    compute_bit_score,
    compute_evalue,
)
from riot_na.data.model import Cigar, InternalAlignmentEntry, InternalAlignmentEntryAA


# Extracted aligner call for profiling purposes.
def _align(aligner, target):
    res = aligner(target)
    return res


def align(
    aligner: StripedSmithWaterman,
    target_id: str,
    target: str,
    db_length: int,
    query: str,
    rev_comp: bool = False,
) -> InternalAlignmentEntry:
    res = _align(aligner, target)
    bit_score = compute_bit_score(res["optimal_alignment_score"])
    e_value = compute_evalue(len(query), db_length, res["optimal_alignment_score"])
    # Need to add 1 because of: https://github.com/biocore/scikit-bio/issues/1340

    q_start = res["query_begin"]
    q_end = res["query_end"] + 1
    t_start = res["target_begin"]
    t_end = res["target_end_optimal"] + 1

    # seq_identity = calculate_seq_identity(res["cigar"], query[q_start:q_end], target[t_start:t_end])
    return InternalAlignmentEntry(
        target_id=target_id,
        alignment_score=bit_score,
        seq_identity=calculate_seq_identity(res["cigar"], query, target, q_start, t_start),
        e_value=e_value,
        q_start=q_start,
        q_end=q_end,
        q_len=q_end - q_start,
        t_start=t_start,
        t_end=t_end,
        t_len=t_end - t_start,
        cigar=Cigar(res["cigar"]),
        query=query,
        rev_comp=rev_comp,
    )


def align_aa(
    aligner: StripedSmithWaterman,
    query: str,
    target_id: str,
    target: str,
    db_length: Optional[int] = None,
    calculate_score: bool = True,
) -> InternalAlignmentEntryAA:
    res = aligner(target)

    q_start = res["query_begin"]
    q_end = res["query_end"] + 1
    t_start = res["target_begin"]
    t_end = res["target_end_optimal"] + 1

    if calculate_score:
        assert db_length is not None
        bit_score = compute_bit_score(res["optimal_alignment_score"], gumbell_params=GumbellParams.AA)
        e_value = compute_evalue(len(query), db_length, res["optimal_alignment_score"])
        seq_identity = calculate_seq_identity(res["cigar"], query, target, q_start, t_start) if query else 0
    else:
        bit_score = None
        e_value = None
        seq_identity = None

    return InternalAlignmentEntryAA(
        target_id=target_id,
        alignment_score=bit_score,
        seq_identity=seq_identity,
        e_value=e_value,
        q_start=q_start,
        q_end=q_end,
        t_start=t_start,
        t_end=t_end,
        cigar=Cigar(res["cigar"]),
    )
