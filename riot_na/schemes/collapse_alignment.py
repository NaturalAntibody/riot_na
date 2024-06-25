from itertools import groupby
from math import ceil
from typing import Sequence

from riot_na.data.model import AlignmentString


def collapse_ins_del(ops: Sequence[str]) -> AlignmentString:
    match_count = ops.count("M")
    deletions = ops.count("D")
    insertions = ops.count("I")
    extra_matches = min(deletions, insertions)
    match_count = match_count + extra_matches
    left_match_count = ceil(match_count / 2)
    right_match_count = match_count - left_match_count
    if deletions > insertions:
        return AlignmentString(f"{left_match_count*'M'}{(deletions-insertions) *'D'}{right_match_count*'M'}")
    if deletions < insertions:
        return AlignmentString(f"{left_match_count*'M'}{(insertions-deletions) *'I'}{right_match_count*'M'}")
    return AlignmentString(f"{match_count*'M'}")


def _collapse_ins_del_ordered(ops: Sequence[str]) -> AlignmentString:
    deletions_indices: list[int] = []
    insertion_indices: list[int] = []
    deletion_indices_to_change_to_match = set()
    insertion_indices_to_remove = set()
    last_op = ""

    for i, op in enumerate(ops):
        match op:
            case "D":
                if last_op == "I":
                    deletions_indices = []
                if insertion_indices:
                    deletion_indices_to_change_to_match.add(i)
                    insertion_indices_to_remove.add(insertion_indices.pop())
                else:
                    deletions_indices.append(i)
                last_op = op
            case "I":
                if last_op == "D":
                    insertion_indices = []
                if deletions_indices:
                    deletion_indices_to_change_to_match.add(deletions_indices.pop())
                    insertion_indices_to_remove.add(i)
                else:
                    insertion_indices.append(i)
                last_op = op

    ops_without_extra_insertions = [
        op if i not in deletion_indices_to_change_to_match else "M"
        for i, op in enumerate(ops)
        if i not in insertion_indices_to_remove
    ]
    return AlignmentString("".join(ops_without_extra_insertions))


def collapse_alignment_str(alignment_str: AlignmentString, ordered: bool = False) -> AlignmentString:
    collapse_ins_del_fn = _collapse_ins_del_ordered if ordered else collapse_ins_del
    res = []
    for key, group in groupby(alignment_str, lambda x: x == "M"):
        ops = list(group)
        if key:
            res.extend(ops)
        else:
            res.append(collapse_ins_del_fn(ops))
    return AlignmentString("".join(res))
