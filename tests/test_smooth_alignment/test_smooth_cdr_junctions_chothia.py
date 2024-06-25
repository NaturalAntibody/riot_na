from riot_na.data.model import AlignmentString, ChainType, Scheme
from riot_na.data.scheme_definitions import (
    CHOTHIA_POSITIONS_HEAVY,
    CHOTHIA_POSITIONS_LIGHT,
)
from riot_na.schemes.smooth_alignment import smooth_cdr_junctions

SCHEME = Scheme.CHOTHIA


def test_smooth_cdr_junctions_chothia_insertions_and_deletions():
    regions = {
        "fwr1": "MMDMDMIIMMDMDDMMDMMDMMIIMMMMIMII",
        "cdr1": "DMDMIIMDM",
        "fwr2": "IMIMMDMDMMMDMMDDDMMMIMI",
        "cdr2": "MDDDM",
        "fwr3": "MDMMMMMDMMMMMMMMDMMMIIMMMMDMMIIMMMMIMMDMDMMIMIII",
        "cdr3": "MDMMMM",
        "fwr4": "IIMMMMDMMMDMMIIM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == CHOTHIA_POSITIONS_HEAVY

    regions = {
        "fwr1": "MMDMDMIIMMDMDDMMDMMDMMIIMMMMIM",
        "cdr1": "MMMMMMIIM",
        "fwr2": "MIMMDMDMMMDMMDDDMMMIM",
        "cdr2": "DDMMM",
        "fwr3": "MDMMMMMDMMMMMMMMDMMMIIMMMMDMMIIMMMMIMMDMDMMIM",
        "cdr3": "MMMMMIIIIM",
        "fwr4": "MMMMDMMMDMMIIM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == CHOTHIA_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=ChainType.HEAVY, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_chothia_insertions_and_deletions_light_chain():
    regions = {
        "fwr1": "MMMMIMMMMDMMMMIIMMDDMMMMMMDMI",
        "cdr1": "MDMMIMMIMI",
        "fwr2": "MDMDMMMMIIMMDMMMMMD",
        "cdr2": "MDM",
        "fwr3": "IMIMMDDMMMMMIIMMMDMMMDMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MDMDDD",
        "fwr4": "MIMMMMIIMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == CHOTHIA_POSITIONS_LIGHT

    regions = {
        "fwr1": "MMMMIMMMMDMMMMIIMMDDMMMMMMDM",
        "cdr1": "MMMMMIIIMM",
        "fwr2": "MDMDMMMMIIMMDMMMMMD",
        "cdr2": "MMM",
        "fwr3": "MIMMDDMMMMMIIMMMDMMMDMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMDDM",
        "fwr4": "MIMMMMIIMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == CHOTHIA_POSITIONS_LIGHT

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=ChainType.LIGHT, scheme=SCHEME)
        == then_alignment_str
    )


if __name__ == "__main__":
    test_smooth_cdr_junctions_chothia_insertions_and_deletions()
    test_smooth_cdr_junctions_chothia_insertions_and_deletions_light_chain()
