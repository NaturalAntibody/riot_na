from functools import cache
from typing import Optional

from riot_na.common.assert_never import assert_never
from riot_na.data.model import ChainType, Scheme, ShortRegion
from riot_na.data.scheme_definitions import (
    CHOTHIA_REGIONS,
    IMGT_REGIONS,
    KABAT_REGIONS,
    MARTIN_REGIONS,
    ChainRegions,
    Regions,
)


class UnknownPosition(Exception):
    pass


@cache
def get_regions_by_scheme_and_chain(scheme: Scheme, chain: ChainType) -> ChainRegions:
    """
    Get the region of a given position in a specified scheme entry.
    """
    match scheme:
        case Scheme.KABAT:
            match chain:
                case ChainType.HEAVY:
                    return KABAT_REGIONS[chain]
                case ChainType.LIGHT:
                    return KABAT_REGIONS[chain]
            assert_never(chain)
        case Scheme.IMGT:
            match chain:
                case ChainType.HEAVY:
                    return IMGT_REGIONS[chain]
                case ChainType.LIGHT:
                    return IMGT_REGIONS[chain]
            assert_never(chain)
        case Scheme.CHOTHIA:
            match chain:
                case ChainType.HEAVY:
                    return CHOTHIA_REGIONS[chain]
                case ChainType.LIGHT:
                    return CHOTHIA_REGIONS[chain]
            assert_never(chain)
        case Scheme.MARTIN:
            match chain:
                case ChainType.HEAVY:
                    return MARTIN_REGIONS[chain]
                case ChainType.LIGHT:
                    return MARTIN_REGIONS[chain]
            assert_never(chain)
    assert_never(scheme)


def position_to_region(regions: ChainRegions, position: int) -> Optional[ShortRegion]:
    if regions[ShortRegion.FW1]["min"] <= position <= regions[ShortRegion.FW1]["max"]:
        return ShortRegion.FW1
    if regions[ShortRegion.CDR1]["min"] <= position <= regions[ShortRegion.CDR1]["max"]:
        return ShortRegion.CDR1
    if regions[ShortRegion.FW2]["min"] <= position <= regions[ShortRegion.FW2]["max"]:
        return ShortRegion.FW2
    if regions[ShortRegion.CDR2]["min"] <= position <= regions[ShortRegion.CDR2]["max"]:
        return ShortRegion.CDR2
    if regions[ShortRegion.FW3]["min"] <= position <= regions[ShortRegion.FW3]["max"]:
        return ShortRegion.FW3
    if regions[ShortRegion.CDR3]["min"] <= position <= regions[ShortRegion.CDR3]["max"]:
        return ShortRegion.CDR3
    if regions[ShortRegion.FW4]["min"] <= position <= regions[ShortRegion.FW4]["max"]:
        return ShortRegion.FW4

    return None


def get_region(scheme_position: int, scheme: Scheme, chain_type: ChainType) -> ShortRegion:
    """
    Based on scheme position, get the region
    """
    region_mapping = get_regions_by_scheme_and_chain(scheme, chain_type)

    region = position_to_region(region_mapping, scheme_position)

    if region:
        return region

    raise UnknownPosition(f"Unknown position {scheme_position} for scheme {scheme} chain {chain_type}")


def get_regions_definitions() -> dict[Scheme, Regions]:
    return {
        Scheme.IMGT: IMGT_REGIONS,
        Scheme.KABAT: KABAT_REGIONS,
        Scheme.CHOTHIA: CHOTHIA_REGIONS,
        Scheme.MARTIN: CHOTHIA_REGIONS,
    }


@cache
def get_cdr_ranges(scheme: Scheme, chain_type: ChainType) -> list[tuple[int, int]]:
    def _get_cdrs(regions: ChainRegions):
        return [
            (regions[region]["min"], regions[region]["max"])
            for region in [ShortRegion.CDR1, ShortRegion.CDR2, ShortRegion.CDR3]
        ]

    match (scheme, chain_type):
        case Scheme.IMGT, ChainType.HEAVY:
            return _get_cdrs(IMGT_REGIONS[chain_type])
        case Scheme.IMGT, ChainType.LIGHT:
            return _get_cdrs(IMGT_REGIONS[chain_type])
        case Scheme.KABAT, ChainType.HEAVY:
            return _get_cdrs(KABAT_REGIONS[chain_type])
        case Scheme.KABAT, ChainType.LIGHT:
            return _get_cdrs(KABAT_REGIONS[chain_type])
        case Scheme.CHOTHIA, ChainType.HEAVY:
            return _get_cdrs(CHOTHIA_REGIONS[chain_type])
        case Scheme.CHOTHIA, ChainType.LIGHT:
            return _get_cdrs(CHOTHIA_REGIONS[chain_type])
        case Scheme.MARTIN, ChainType.HEAVY:
            return _get_cdrs(MARTIN_REGIONS[chain_type])
        case Scheme.MARTIN, ChainType.LIGHT:
            return _get_cdrs(MARTIN_REGIONS[chain_type])

    raise ValueError(f"Unknown scheme {scheme} or chain {chain_type}")


if __name__ == "__main__":
    print(list(get_cdr_ranges(Scheme.IMGT, ChainType.HEAVY)))
