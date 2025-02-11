from typing import Optional

from riot_na.data.model import (
    ChainType,
    RegionOffsetsAA,
    RegionOffsetsNT,
    Scheme,
    SchemeAlignment,
    ShortRegion,
)
from riot_na.data.scheme_regions import get_region


def infer_region_offsets(
    scheme_alignment: SchemeAlignment, scheme: Scheme, chain_type: ChainType
) -> dict[ShortRegion, list[int]]:
    result: dict[ShortRegion, list[int]] = {}

    scheme_position = 0
    query_position = scheme_alignment.q_start or 0
    insertion_counter = 0
    alignment_str = scheme_alignment.alignment_str

    region_id = ShortRegion.FW1

    for op in alignment_str:
        if op == "M":
            insertion_counter = 0
            query_position = query_position + 1
            scheme_position = scheme_position + 1
            region_id = get_region(scheme_position, scheme, chain_type)

            residues = result.get(region_id, [])
            residues.append(query_position)
            result[region_id] = residues

        elif op in {"D", "N"}:
            insertion_counter = 0
            scheme_position = scheme_position + 1
            region_id = get_region(scheme_position, scheme, chain_type)
        elif op == "I":
            insertion_counter = insertion_counter + 1
            query_position = query_position + 1
            residues = result.get(region_id, [])
            residues.append(query_position)
            result[region_id] = residues

    return result


def get_first_aa(region_positions: Optional[list[int]]) -> int:
    if not region_positions or len(region_positions) == 0:
        return -1

    return region_positions[0]


def get_last_aa(region_positions: Optional[list[int]]) -> int:
    if not region_positions or len(region_positions) == 0:
        return -1

    return region_positions[-1]


def infer_aa_region_offsets(aa_regions: dict[ShortRegion, list[int]]) -> RegionOffsetsAA:
    return RegionOffsetsAA(
        fwr1_start_aa=get_first_aa(aa_regions.get(ShortRegion.FW1)),
        fwr1_end_aa=get_last_aa(aa_regions.get(ShortRegion.FW1)),
        cdr1_start_aa=get_first_aa(aa_regions.get(ShortRegion.CDR1)),
        cdr1_end_aa=get_last_aa(aa_regions.get(ShortRegion.CDR1)),
        fwr2_start_aa=get_first_aa(aa_regions.get(ShortRegion.FW2)),
        fwr2_end_aa=get_last_aa(aa_regions.get(ShortRegion.FW2)),
        cdr2_start_aa=get_first_aa(aa_regions.get(ShortRegion.CDR2)),
        cdr2_end_aa=get_last_aa(aa_regions.get(ShortRegion.CDR2)),
        fwr3_start_aa=get_first_aa(aa_regions.get(ShortRegion.FW3)),
        fwr3_end_aa=get_last_aa(aa_regions.get(ShortRegion.FW3)),
        cdr3_start_aa=get_first_aa(aa_regions.get(ShortRegion.CDR3)),
        cdr3_end_aa=get_last_aa(aa_regions.get(ShortRegion.CDR3)),
        fwr4_start_aa=get_first_aa(aa_regions.get(ShortRegion.FW4)),
        fwr4_end_aa=get_last_aa(aa_regions.get(ShortRegion.FW4)),
    )


def get_first_nt(region: Optional[list[int]], offset: int) -> int:
    if region is None:
        return -1

    return region[0] * 3 - 2 + offset


def get_last_nt(region: Optional[list[int]], offset: int) -> int:
    if region is None:
        return -1

    return region[-1] * 3 + offset


def infer_nt_region_offsets(
    aa_regions: dict[ShortRegion, list[int]], nt_alignment_start: int = 0, reading_frame: int = 0
) -> RegionOffsetsNT:
    # nt     1 2 3 4 5 6 7
    # pt       1     2
    # start  |
    # end              |
    offset = nt_alignment_start + reading_frame

    return RegionOffsetsNT(
        fwr1_start=get_first_nt(aa_regions.get(ShortRegion.FW1), offset),
        fwr1_end=get_last_nt(aa_regions.get(ShortRegion.FW1), offset),
        cdr1_start=get_first_nt(aa_regions.get(ShortRegion.CDR1), offset),
        cdr1_end=get_last_nt(aa_regions.get(ShortRegion.CDR1), offset),
        fwr2_start=get_first_nt(aa_regions.get(ShortRegion.FW2), offset),
        fwr2_end=get_last_nt(aa_regions.get(ShortRegion.FW2), offset),
        cdr2_start=get_first_nt(aa_regions.get(ShortRegion.CDR2), offset),
        cdr2_end=get_last_nt(aa_regions.get(ShortRegion.CDR2), offset),
        fwr3_start=get_first_nt(aa_regions.get(ShortRegion.FW3), offset),
        fwr3_end=get_last_nt(aa_regions.get(ShortRegion.FW3), offset),
        cdr3_start=get_first_nt(aa_regions.get(ShortRegion.CDR3), offset),
        cdr3_end=get_last_nt(aa_regions.get(ShortRegion.CDR3), offset),
        fwr4_start=get_first_nt(aa_regions.get(ShortRegion.FW4), offset),
        fwr4_end=get_last_nt(aa_regions.get(ShortRegion.FW4), offset),
    )
