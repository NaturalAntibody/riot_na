from functools import cache
from typing import Final, TypedDict

from riot_na.data.model import ChainType, Scheme, ShortRegion

# Number of scheme positions for each chain type
IMGT_POSITIONS: Final[int] = 128
KABAT_POSITIONS_LIGHT: Final[int] = 107
KABAT_POSITIONS_HEAVY: Final[int] = 113
CHOTHIA_POSITIONS_LIGHT: Final[int] = 107
CHOTHIA_POSITIONS_HEAVY: Final[int] = 113
MARTIN_POSITIONS_LIGHT: Final[int] = 107
MARTIN_POSITIONS_HEAVY: Final[int] = 113

# CDR insertions / deletions positions
IMGT_INDEL_POSITIONS: Final[list[int]] = [32, 60, 111]
KABAT_INDEL_POSITIONS_HEAVY: Final[list[int]] = [35, 52, 100]
KABAT_INDEL_POSITIONS_LIGHT: Final[list[int]] = [27, 52, 95]
CHOTHIA_INDEL_POSITIONS_HEAVY: Final[list[int]] = [31, 52, 100]
CHOTHIA_INDEL_POSITIONS_LIGHT: Final[list[int]] = [30, 52, 95]
MARTIN_INDEL_POSITIONS_HEAVY: Final[list[int]] = [31, 52, 100]
MARTIN_INDEL_POSITIONS_LIGHT: Final[list[int]] = [30, 52, 95]


# Region position ranges
class MinMaxRange(TypedDict):
    min: int
    max: int


ChainRegions = dict[ShortRegion, MinMaxRange]
Regions = dict[ChainType, ChainRegions]

IMGT_REGIONS: Final[Regions] = {
    ChainType.HEAVY: {
        ShortRegion.FW1: {"min": 1, "max": 26},
        ShortRegion.CDR1: {"min": 27, "max": 38},
        ShortRegion.FW2: {"min": 39, "max": 55},
        ShortRegion.CDR2: {"min": 56, "max": 65},
        ShortRegion.FW3: {"min": 66, "max": 104},
        ShortRegion.CDR3: {"min": 105, "max": 117},
        ShortRegion.FW4: {"min": 118, "max": 128},
    },
    ChainType.LIGHT: {
        ShortRegion.FW1: {"min": 1, "max": 26},
        ShortRegion.CDR1: {"min": 27, "max": 38},
        ShortRegion.FW2: {"min": 39, "max": 55},
        ShortRegion.CDR2: {"min": 56, "max": 65},
        ShortRegion.FW3: {"min": 66, "max": 104},
        ShortRegion.CDR3: {"min": 105, "max": 117},
        ShortRegion.FW4: {"min": 118, "max": 128},
    },
}

KABAT_REGIONS: Final[Regions] = {
    ChainType.HEAVY: {
        ShortRegion.FW1: {"min": 1, "max": 30},
        ShortRegion.CDR1: {"min": 31, "max": 35},
        ShortRegion.FW2: {"min": 36, "max": 49},
        ShortRegion.CDR2: {"min": 50, "max": 65},
        ShortRegion.FW3: {"min": 66, "max": 94},
        ShortRegion.CDR3: {"min": 95, "max": 102},
        ShortRegion.FW4: {"min": 103, "max": 113},
    },
    ChainType.LIGHT: {
        ShortRegion.FW1: {"min": 1, "max": 23},
        ShortRegion.CDR1: {"min": 24, "max": 34},
        ShortRegion.FW2: {"min": 35, "max": 49},
        ShortRegion.CDR2: {"min": 50, "max": 56},
        ShortRegion.FW3: {"min": 57, "max": 88},
        ShortRegion.CDR3: {"min": 89, "max": 97},
        ShortRegion.FW4: {"min": 98, "max": 107},
    },
}


CHOTHIA_REGIONS: Final[Regions] = {
    ChainType.HEAVY: {
        ShortRegion.FW1: {"min": 1, "max": 25},
        ShortRegion.CDR1: {"min": 26, "max": 32},
        ShortRegion.FW2: {"min": 33, "max": 51},
        ShortRegion.CDR2: {"min": 52, "max": 56},
        ShortRegion.FW3: {"min": 57, "max": 95},
        ShortRegion.CDR3: {"min": 96, "max": 101},
        ShortRegion.FW4: {"min": 102, "max": 113},
    },
    ChainType.LIGHT: {
        ShortRegion.FW1: {"min": 1, "max": 25},
        ShortRegion.CDR1: {"min": 26, "max": 32},
        ShortRegion.FW2: {"min": 33, "max": 49},
        ShortRegion.CDR2: {"min": 50, "max": 52},
        ShortRegion.FW3: {"min": 53, "max": 90},
        ShortRegion.CDR3: {"min": 91, "max": 96},
        ShortRegion.FW4: {"min": 97, "max": 107},
    },
}

MARTIN_REGIONS: Final[Regions] = {
    ChainType.HEAVY: {
        ShortRegion.FW1: {"min": 1, "max": 25},
        ShortRegion.CDR1: {"min": 26, "max": 32},
        ShortRegion.FW2: {"min": 33, "max": 51},
        ShortRegion.CDR2: {"min": 52, "max": 56},
        ShortRegion.FW3: {"min": 57, "max": 95},
        ShortRegion.CDR3: {"min": 96, "max": 101},
        ShortRegion.FW4: {"min": 102, "max": 113},
    },
    ChainType.LIGHT: {
        ShortRegion.FW1: {"min": 1, "max": 25},
        ShortRegion.CDR1: {"min": 26, "max": 32},
        ShortRegion.FW2: {"min": 33, "max": 49},
        ShortRegion.CDR2: {"min": 50, "max": 52},
        ShortRegion.FW3: {"min": 53, "max": 90},
        ShortRegion.CDR3: {"min": 91, "max": 96},
        ShortRegion.FW4: {"min": 97, "max": 107},
    },
}


@cache
def get_legal_positions(chain_type: ChainType, scheme: Scheme) -> int:
    match scheme, chain_type:
        case Scheme.IMGT, ChainType.HEAVY:
            return IMGT_POSITIONS
        case Scheme.IMGT, ChainType.LIGHT:
            return IMGT_POSITIONS
        case Scheme.KABAT, ChainType.HEAVY:
            return KABAT_POSITIONS_HEAVY
        case Scheme.KABAT, ChainType.LIGHT:
            return KABAT_POSITIONS_LIGHT
        case Scheme.CHOTHIA, ChainType.HEAVY:
            return CHOTHIA_POSITIONS_HEAVY
        case Scheme.CHOTHIA, ChainType.LIGHT:
            return CHOTHIA_POSITIONS_LIGHT
        case Scheme.MARTIN, ChainType.HEAVY:
            return MARTIN_POSITIONS_HEAVY
        case Scheme.MARTIN, ChainType.LIGHT:
            return MARTIN_POSITIONS_LIGHT
        case _:
            raise ValueError(f"Unknown scheme {scheme} or chain type {chain_type}")


@cache
def get_indels_positions(scheme: Scheme, chain_type: ChainType) -> list[int]:
    # this is scheme specific
    match (scheme, chain_type):
        case Scheme.IMGT, _:
            return IMGT_INDEL_POSITIONS
        case Scheme.KABAT, ChainType.HEAVY:
            return KABAT_INDEL_POSITIONS_HEAVY
        case Scheme.KABAT, ChainType.LIGHT:
            return KABAT_INDEL_POSITIONS_LIGHT
        case Scheme.CHOTHIA, ChainType.HEAVY:
            return CHOTHIA_INDEL_POSITIONS_HEAVY
        case Scheme.CHOTHIA, ChainType.LIGHT:
            return CHOTHIA_INDEL_POSITIONS_LIGHT
        case Scheme.MARTIN, ChainType.HEAVY:
            return MARTIN_INDEL_POSITIONS_HEAVY
        case Scheme.MARTIN, ChainType.LIGHT:
            return MARTIN_INDEL_POSITIONS_LIGHT

    raise ValueError(f"Unknown scheme {scheme}")
