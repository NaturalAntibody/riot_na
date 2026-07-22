from functools import cache

from riot_na import ChainType, Locus, Scheme
from riot_na.data.model import AirrRearrangementEntryAA, ShortRegion
from riot_na.data.scheme_regions import get_region


def get_region_position_indices(airr: AirrRearrangementEntryAA, region: ShortRegion) -> list[int]:
    """
    Return zero-based amino acid indices for a single region in an AIRR entry.

    Parameters:
    -----------
    - airr : AIRR rearrangement entry with region boundary coordinates.
    - region : CDR or framework region to resolve.

    Returns:
    -----------
    List of zero-based indices covering the requested region.
    """
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
    """
    Return zero-based amino acid indices for multiple regions in an AIRR entry.

    Parameters:
    -----------
    - airr : AIRR rearrangement entry with region boundary coordinates.
    - regions : CDR and/or framework regions to resolve.

    Returns:
    -----------
    Sorted list of unique zero-based indices across all requested regions.
    """
    indices = set()
    for region in regions:
        indices.update(get_region_position_indices(airr, region))
    return sorted(list(indices))


def scheme_positions_to_index(airr: AirrRearrangementEntryAA, scheme_positions: list[str]) -> list[int]:
    """
    Map scheme position labels to zero-based sequence indices.

    Parameters:
    -----------
    - airr : AIRR rearrangement entry with positional_scheme_mapping populated.
    - scheme_positions : Scheme position labels to look up (e.g. "111A").

    Returns:
    -----------
    Zero-based indices corresponding to the given scheme positions.
    """
    assert airr.positional_scheme_mapping
    scheme_positional_mapping = {v: k for k, v in airr.positional_scheme_mapping.items()}
    return [scheme_positional_mapping[pos] for pos in scheme_positions]


def get_primary_seq(airr: AirrRearrangementEntryAA) -> str:
    """
    Reconstruct the numbered amino acid sequence from scheme residue mapping.

    Parameters:
    -----------
    - airr : AIRR rearrangement entry with scheme_residue_mapping populated.

    Returns:
    -----------
    Amino acid sequence formed by joining mapped residues in scheme order.
    """
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
    """
    Convert numeric insertion notation to IMGT-style position labels.

    Parameters:
    -----------
    - position : Scheme position, optionally with a numeric insertion suffix
      (e.g. "111.1", "111.2").

    Returns:
    -----------
    Position label with insertion letters (e.g. "111A", "111B"), or the
    original value when no insertion suffix is present.
    """
    if "." not in position:
        return position
    position_number, insertion_number = position.split(".")
    return f"{position_number}{int_to_str_insertion(int(insertion_number))}"


def str_to_int_insertion(position: str) -> int:
    """
    Extract the integer scheme position from a position label.

    Parameters:
    -----------
    - position : Scheme position label, with or without an insertion suffix
      (e.g. "111", "111.1", "112.1", "112.2", etc.).

    Returns:
    -----------
    Integer scheme position, ignoring any insertion suffix.
    """
    return int(position.split(".", maxsplit=1)[0])


def get_scheme_residue_mapping_insertions_as_letters(
    airr: AirrRearrangementEntryAA,
) -> dict[str, str]:
    """
    Convert scheme residue mapping insertions to IMGT-style position labels (e.g. "111.1" -> "111A", "111.2" -> "111B", etc.).

    Parameters:
    -----------
    - airr : AIRR rearrangement entry with scheme_residue_mapping populated.

    Returns:
    -----------
    Scheme residue mapping with insertions converted to letters (e.g. {"111A": "G", "111B": "K"}).
    """
    assert airr.scheme_residue_mapping is not None
    return {
        map_insertion_number_to_letter(scheme_position): residue
        for scheme_position, residue in airr.scheme_residue_mapping.items()
    }


def get_scheme_residue_mapping_by_region(
    airr: AirrRearrangementEntryAA,
) -> dict[ShortRegion, dict[str, str]]:
    """
    Group scheme residue mapping entries by CDR/framework region.

    Parameters:
    -----------
    - airr : AIRR rearrangement entry with scheme_residue_mapping, locus,
      and numbering_scheme populated.

    Returns:
    -----------
    Mapping from each ShortRegion to its scheme-position-to-residue entries.
    """
    assert airr.scheme_residue_mapping is not None
    assert airr.locus is not None
    assert airr.numbering_scheme is not None

    chain_type = ChainType.from_locus(Locus(airr.locus))
    region_scheme_residue_mapping: dict[ShortRegion, dict[str, str]] = {}
    for scheme_position, residue in airr.scheme_residue_mapping.items():
        region = get_region(
            str_to_int_insertion(scheme_position),
            Scheme(airr.numbering_scheme),
            chain_type,
        )
        region_scheme_residue_mapping.setdefault(region, {})[scheme_position] = residue
    return region_scheme_residue_mapping
