from functools import partial
from math import ceil, floor
from typing import Callable, Final, Literal, cast

from riot_na.data.model import AlignmentString, ChainType, Scheme
from riot_na.data.scheme_definitions import get_indels_positions
from riot_na.data.scheme_regions import get_cdr_ranges
from riot_na.schemes.collapse_alignment import collapse_ins_del

Op = Literal["M", "D", "N", "I", "S"]
M: Final[Op] = "M"
D: Final[Op] = "D"
N: Final[Op] = "N"
I: Final[Op] = "I"
S: Final[Op] = "S"


def _reorder_cdr_deletions(relative_indel_position: int, matches: int, remaining_deletions: int) -> str:
    # deletions should be put to the left of specified position
    # indel_position = matches_before + deletions => matches_before = indel_position - deletions
    matches_before = relative_indel_position - remaining_deletions
    if matches_before >= 0:
        assert matches_before <= matches
        return f"{matches_before * M}{remaining_deletions * D}{(matches - matches_before) * M}"

    # if there are too many deletions, just put deletions first, then matches
    return f"{remaining_deletions * D}{matches  * M}"


def _reorder_cdr_deletions_imgt(relative_indel_position: int, matches: int, remaining_deletions: int) -> str:
    # deletions should be put in the middle of cdr, starting from indel_position + 1
    # for cdr1 and cdr2, and indel_position for cdr3. We can detect cdr number from
    # alignment string length, cdr1 and cdr2 have even number of positions, while cdr3 odd.
    legal_positions = matches + remaining_deletions
    round_fn = ceil if legal_positions % 2 else floor
    matches_before = relative_indel_position - round_fn(remaining_deletions / 2)
    matches_before = max(matches_before, 0)
    assert matches >= matches_before
    return f"{matches_before * M}{remaining_deletions * D}{(matches - matches_before) * M}"


def _reorder_cdr(
    alignment_str: AlignmentString,
    relative_indel_position: int,
    reorder_buffer_deletions: Callable[[int, int, int], str],
) -> AlignmentString:
    matches = alignment_str.count(M)
    dels = alignment_str.count(D)
    ins = alignment_str.count(I)

    assert matches + dels + ins == len(alignment_str)

    extra_matches = min(ins, dels)
    matches = matches + extra_matches

    # collapse insertions and deletions over the alignment
    if ins > dels:
        remaining_insertions = ins - dels
        # insertions should be put after specified position or in the middle
        assert relative_indel_position <= matches
        result = f"{relative_indel_position * M}{remaining_insertions * I}{(matches - relative_indel_position) * M}"

    elif dels > ins:
        remaining_deletions = dels - ins
        result = reorder_buffer_deletions(relative_indel_position, matches, remaining_deletions)
    else:
        result = matches * M

    return AlignmentString(result)


reorder_buffer_other_schemes = cast(
    Callable[[AlignmentString, int], AlignmentString],
    partial(_reorder_cdr, reorder_buffer_deletions=_reorder_cdr_deletions),
)
reorder_buffer_imgt = cast(
    Callable[[AlignmentString, int], AlignmentString],
    partial(_reorder_cdr, reorder_buffer_deletions=_reorder_cdr_deletions_imgt),
)


def smooth_cdr_junctions(
    alignment_str: AlignmentString, chain_type: ChainType, scheme: Scheme = Scheme.IMGT
) -> AlignmentString:
    """
    Function for reordering function for reordering cdr insertions and deletions
    to match numbering scheme standard.

    Each numbering scheme (IMGT, KABAT, CHOTHIA, MARTIN) specifies legal positions
    in CDR regions, where insertions or deletions should be put. This function takes
    each CDR and reorders alignment string to match those specifications.

    IMGT standard: https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
    KABAT, CHOTHIA, MARTIN: http://www.bioinf.org.uk/abs/info.html
    """

    reorder_buffer_fn: Callable[[AlignmentString, int], AlignmentString] = (
        reorder_buffer_imgt if scheme == Scheme.IMGT else reorder_buffer_other_schemes
    )
    scheme_position = 0
    smoothed = ""
    buffer = ""
    regions_it = iter(get_cdr_ranges(scheme, chain_type))
    indel_positions_it = iter(get_indels_positions(scheme, chain_type))

    region_start, region_end = next(regions_it, (None, None))
    indel_position = next(indel_positions_it, None)

    for op in alignment_str:
        if op in (M, D, N):
            scheme_position = scheme_position + 1

        # Append the rest of fwr4
        if not (region_start and region_end and indel_position):
            smoothed = smoothed + op
            continue

        # Take only cdr + insertions right before cdr start
        if (region_start - 1 == scheme_position and op == I) or region_start <= scheme_position <= region_end:
            buffer = buffer + op
        else:
            if buffer:
                # For cases when alignment starts on framework, do nothing.
                if set(buffer) == {N}:
                    smoothed = smoothed + buffer
                else:
                    # +1 because of 1-based indexing
                    relative_indel_position = indel_position - region_start + 1
                    fixed_buffer = reorder_buffer_fn(
                        AlignmentString(buffer.replace(N, D)),
                        relative_indel_position,
                    )
                    smoothed = smoothed + fixed_buffer

                buffer = ""
                region_start, region_end = next(regions_it, (None, None))
                indel_position = next(indel_positions_it, None)

            smoothed = smoothed + op

    # For cases when alignment ends on cdr
    if buffer:
        assert N not in buffer, "No framework aligned"
        assert region_start and region_end and indel_position

        region_length = region_end - region_start + 1
        collapsed_buffer = collapse_ins_del(buffer)
        missing_deletions_count = region_length - (collapsed_buffer.count(M) + collapsed_buffer.count(D))
        filled_buffer = AlignmentString(f"{collapsed_buffer}{missing_deletions_count * D}")
        relative_indel_position = indel_position - region_start + 1
        fixed_buffer = reorder_buffer_fn(filled_buffer, relative_indel_position)
        smoothed = smoothed + fixed_buffer

    return AlignmentString(smoothed)
