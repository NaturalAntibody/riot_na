from functools import cache
from typing import Final, TypedDict

from riot_na.data.model import ChainType, CScheme, Scheme, ShortRegion

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


# Constant region scheme positions, in sequence order. The C scheme mappings of germlines
# (scheme_mappings/<species>/<c_scheme>/c_scheme_mapping.csv) are cigars over these positions: a residue on a position
# is a match, a missing position is a deletion and a residue between positions is an insertion.
EU_C_POSITIONS_HEAVY: Final[tuple[str, ...]] = tuple(str(position) for position in range(118, 448))
# EU numbering is defined for kappa light chains only
EU_C_POSITIONS_LIGHT: Final[tuple[str, ...]] = tuple(str(position) for position in range(108, 215))
KABAT_C_POSITIONS_HEAVY: Final[tuple[str, ...]] = tuple(str(position) for position in range(114, 479))
# lambda C starts with an insertion after 107 (Kabat 107A), so 107 is the first position
KABAT_C_POSITIONS_LIGHT: Final[tuple[str, ...]] = tuple(str(position) for position in range(107, 217))

# IMGT unique numbering for C-DOMAIN. Positions with a dot (1.8-1.1, 15.1-15.3, 45.1-45.7, 84.1-84.7, 85.7-85.1,
# 96.1-96.2) are insertions: 1.x precede position 1, 85.x precede 85, the others follow their position.
IMGT_C_DOMAIN_POSITIONS: Final[tuple[str, ...]] = (
    *(str(i) for i in range(1, 32)),
    *(str(i) for i in range(34, 46)),
    *(str(i) for i in range(77, 129)),
)
IMGT_C_DOMAINS: Final[tuple[str, ...]] = ("CH1", "CH2", "CH3", "CL")
# hinge positions are IMGT exon numbering; IGHG3 hinge is encoded by 4 exons (H1-H4), other IGHG by one (H);
# positions 129 and 130 of CH3 are CHS
IMGT_C_POSITIONS_HEAVY: Final[tuple[str, ...]] = (
    *(f"CH1_{position}" for position in IMGT_C_DOMAIN_POSITIONS),
    *(f"H_{i}" for i in range(1, 16)),
    *(f"H{exon}_{i}" for exon, exon_length in ((1, 17), (2, 15), (3, 15), (4, 15)) for i in range(1, exon_length + 1)),
    *(f"CH2_{position}" for position in IMGT_C_DOMAIN_POSITIONS),
    *(f"CH3_{position}" for position in IMGT_C_DOMAIN_POSITIONS),
    "CH3_129",
    "CH3_130",
)
IMGT_C_POSITIONS_LIGHT: Final[tuple[str, ...]] = tuple(f"CL_{position}" for position in IMGT_C_DOMAIN_POSITIONS)


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


def get_c_scheme_positions(c_scheme: CScheme, chain_type: ChainType) -> tuple[str, ...]:
    match c_scheme, chain_type:
        case CScheme.EU, ChainType.HEAVY:
            return EU_C_POSITIONS_HEAVY
        case CScheme.EU, ChainType.LIGHT:
            return EU_C_POSITIONS_LIGHT
        case CScheme.KABAT, ChainType.HEAVY:
            return KABAT_C_POSITIONS_HEAVY
        case CScheme.KABAT, ChainType.LIGHT:
            return KABAT_C_POSITIONS_LIGHT
        case CScheme.IMGT, ChainType.HEAVY:
            return IMGT_C_POSITIONS_HEAVY
        case CScheme.IMGT, ChainType.LIGHT:
            return IMGT_C_POSITIONS_LIGHT

    raise ValueError(f"Unknown C scheme {c_scheme} or chain type {chain_type}")
