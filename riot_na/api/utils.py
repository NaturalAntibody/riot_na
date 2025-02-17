from functools import cache

from riot_na.data.model import AirrRearrangementEntryAA, ShortRegion


def get_region_position_indices(airr: AirrRearrangementEntryAA, region: ShortRegion) -> list[int]:
    match region:
        case ShortRegion.CDR1:
            assert airr.cdr1_start_aa and airr.cdr1_end_aa
            return list(range(airr.cdr1_start_aa - 1, airr.cdr1_end_aa))
        case ShortRegion.CDR2:
            assert airr.cdr2_start_aa and airr.cdr2_end_aa
            return list(range(airr.cdr2_start_aa - 1, airr.cdr2_end_aa))
        case ShortRegion.CDR3:
            assert airr.cdr3_start_aa and airr.cdr3_end_aa
            return list(range(airr.cdr3_start_aa - 1, airr.cdr3_end_aa))
        case ShortRegion.FW1:
            assert airr.fwr1_start_aa and airr.fwr1_end_aa
            return list(range(airr.fwr1_start_aa - 1, airr.fwr1_end_aa))
        case ShortRegion.FW2:
            assert airr.fwr2_start_aa and airr.fwr2_end_aa
            return list(range(airr.fwr2_start_aa - 1, airr.fwr2_end_aa))
        case ShortRegion.FW3:
            assert airr.fwr3_start_aa and airr.fwr3_end_aa
            return list(range(airr.fwr3_start_aa - 1, airr.fwr3_end_aa))
        case ShortRegion.FW4:
            assert airr.fwr4_start_aa and airr.fwr4_end_aa
            return list(range(airr.fwr4_start_aa - 1, airr.fwr4_end_aa))


def get_regions_position_indices(airr: AirrRearrangementEntryAA, regions: list[ShortRegion]) -> list[int]:
    indices = set()
    for region in regions:
        indices.update(get_region_position_indices(airr, region))
    return sorted(list(indices))


def scheme_positions_to_index(airr: AirrRearrangementEntryAA, scheme_positions: list[str]) -> list[int]:
    assert airr.positional_scheme_mapping
    scheme_positional_mapping = {v: k for k, v in airr.positional_scheme_mapping.items()}
    return [scheme_positional_mapping[pos] for pos in scheme_positions]


def get_primary_seq(airr: AirrRearrangementEntryAA) -> str:
    assert airr.scheme_residue_mapping is not None
    return "".join(airr.scheme_residue_mapping.values())


@cache
def int_to_str_insertion(n: int) -> str:
    """
    Converts an integer (1-based) to an IMGT-style insertion letter.
    Example:
      1 -> 'A', 26 -> 'Z', 27 -> 'AA', 28 -> 'AB', ..., 52 -> 'AZ', 53 -> 'BA'
    """
    if n < 1:
        raise ValueError("Input must be a positive integer.")

    result = ""
    while n > 0:
        n -= 1  # Adjust for 1-based indexing
        result = chr(ord("A") + (n % 26)) + result
        n //= 26
    return result


def map_insertion_number_to_letter(position: str) -> str:
    # Translates 111.1, 111.2 to 111A, 111B ETC
    if "." not in position:
        return position
    position_number, insertion_number = position.split(".")
    return f"{position_number}{int_to_str_insertion(int(insertion_number))}"
