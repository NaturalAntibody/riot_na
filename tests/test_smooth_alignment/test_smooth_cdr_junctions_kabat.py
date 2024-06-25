from riot_na.data.model import AlignmentString, ChainType, Scheme
from riot_na.data.scheme_definitions import KABAT_POSITIONS_HEAVY, KABAT_POSITIONS_LIGHT
from riot_na.schemes.smooth_alignment import (
    _reorder_cdr_deletions,
    smooth_cdr_junctions,
)

SCHEME = Scheme.KABAT
CHAIN_TYPE = ChainType.HEAVY


def test_smooth_cdr_junctions_kabat_no_indels():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == given_alignment_str
    )


def test_smooth_cdr_junctions_kabat_deletions_in_frameworks():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMD",
        "cdr1": "MMMMM",
        "fwr2": "DMMMMMMMMMMMMD",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MDMMMMMMMMMMMMMMMMMMMMMMMMMMD",
        "cdr3": "MMMMMMMM",
        "fwr4": "DMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    then_alignment_str = given_alignment_str

    temp = smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
    print(temp)
    assert temp == then_alignment_str


def test_smooth_cdr_junctions_kabat_deletions_in_cdrs():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMD",
        "cdr1": "DDMMM",
        "fwr2": "DMMMMMMMMMMMMD",
        "cdr2": "MMMMMMMDMMMDMMMM",
        "fwr3": "MDMMMMMMMMMMMMMMMMMMMMMMMMMMD",
        "cdr3": "DMMMMMDM",
        "fwr4": "DMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMD",
        "cdr1": "MMMDD",
        "fwr2": "DMMMMMMMMMMMMD",
        "cdr2": "MDDMMMMMMMMMMMMM",
        "fwr3": "MDMMMMMMMMMMMMMMMMMMMMMMMMMMD",
        "cdr3": "MMMMDDMM",
        "fwr4": "DMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    temp = smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
    assert temp == then_alignment_str


def test_smooth_cdr_junctions_kabat_too_much_deletions():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "DDMDD",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMDDMDMMD",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DMDDDMDD",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MDDDD",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "DDDDMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DDDDDDMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_insertions_in_frameworks():
    regions = {
        "fwr1": "MMMMMIIMMMMIMMMMMMMIMMMMMMMMMMMMMIM",
        "cdr1": "MMMMM",
        "fwr2": "MIMMMIMMMIMMMMIMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MIMMIIMMMMMMMMMMMMMMMMMIIMMMMMMMMIM",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMMIMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    then_alignment_str = given_alignment_str
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_insertions_in_cdrs():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMII",
        "cdr1": "MIMMIMMII",
        "fwr2": "IIMMMMMMMMMMMMMMII",
        "cdr2": "MIIMMMMMMIIMMMMMMMMIMI",
        "fwr3": "IMMMMMMMMMMMMMMMMMMMMMMMMMMMMMI",
        "cdr3": "IMMMMIMMMIIM",
        "fwr4": "IMMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMIIIIIIII",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMIIIIIIIIIMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMIIIIIIMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_insertions_and_deletions():
    regions = {
        "fwr1": "MMMMIMMMMDMMMMIIMMDDMMMMMMMMMMMDMI",
        "cdr1": "MDMIMIMI",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMDMMMIMMIMDMMMMIMM",
        "fwr3": "IMIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MDMMMMMM",
        "fwr4": "MIMMMMIIMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "MMMMIMMMMDMMMMIIMMDDMMMMMMMMMMMDM",
        "cdr1": "MMMMMIII",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMMIIMMMMMMMMMMMMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIMM",
        "fwr4": "MIMMMMIIMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_start_on_last_fwr1():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNM",
        "cdr1": "MDMIMIMI",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMDMMMIMMIMDMMMMIMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNM",
        "cdr1": "MMMMMII",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMMIMMMMMMMMMMMMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIIMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_start_on_cdr1_start():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MDMIMIMI",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMDMMMIMMIMDMMMMIMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MMMMMII",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMMIMMMMMMMMMMMMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIIMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_start_on_cdr1():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNIMIM",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMDMMMIMMIMDMMMMIMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MMMMD",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMMIMMMMMMMMMMMMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIIMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_start_on_cdr1_end():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNM",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMDMMMIMMIMDMMMMIMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MDDDD",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMMIMMMMMMMMMMMMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIIMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_start_on_fwr2():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNN",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMDMMMIMMIMDMMMMIMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNN",
        "fwr2": "MDMDMMMIIMMDMMMD",
        "cdr2": "MMMIMMMMMMMMMMMMM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIIMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_start_on_cdr2_end():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNN",
        "fwr2": "NNNNNNNNNNNNNN",
        "cdr2": "NNNNNNNNNNNNNNNM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMMII",
        "cdr3": "MMMMMMMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNN",
        "fwr2": "NNNNNNNNNNNNNN",
        "cdr2": "DDDDDDDDDDDDDDDM",
        "fwr3": "MIMMMMIIMMMMMMMMMIIIMMMMDMMMMDMMMMM",
        "cdr3": "MMMMMMIIMM",
        "fwr4": "MIMMMIIMMMDMDM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_HEAVY

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_short_fwr4():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMM",
        "fwr4": "MMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 107

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMM",
        "fwr4": "MMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 107

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_fwr4_junction():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMM",
        "fwr4": "DM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 104

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMM",
        "fwr4": "DM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 104

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_in_fwr4_junction():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MDMMIMMMM",
        "fwr4": "M",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 103

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMM",
        "fwr4": "M",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 103

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_no_fwr4():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DMMMMMMM",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 102

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMDMM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 102

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_cdr3():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DMMMMM",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 100

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMDDDMM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 102

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_cdr3_before_indel_pos():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DMM",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 97

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DDDDDDMM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 102

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_single_match_on_cdr3():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "M",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 95

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "DDDDDDDM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 102

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_fwr3_junction():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 94

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 94

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_fwr3_junction_last_missing():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 93

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 93

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_fwr3():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 87

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 87

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_kabat_end_on_cdr2():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMMMMMDM",
        "fwr3": "",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 64

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMM",
        "fwr2": "MMMMMMMMMMMMMM",
        "cdr2": "MDDMMMMMMMMMMMMM",
        "fwr3": "",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 65

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_reorder_buffer_deletions_kabat_no_matches():
    # even buffer len
    relative_indel_position = 6
    buffer = 14 * "D"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 15 * "D"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer


def test_reorder_buffer_deletions_kabat_single_match_at_buffer_end():
    # even buffer len
    relative_indel_position = 6
    buffer = 13 * "D" + "M"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 14 * "D" + "M"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer


def test_reorder_buffer_deletions_kabat_match_on_each_buffer_end():
    # even buffer len
    relative_indel_position = 6
    buffer = "M" + 13 * "D" + "M"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == 13 * "D" + "MM"

    # odd buffer len
    relative_indel_position = 7
    buffer = "M" + 14 * "D" + "M"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == 14 * "D" + "MM"


def test_reorder_buffer_deletions_kabat_two_matches_on_buffer_right_end():
    # even buffer len
    relative_indel_position = 6
    buffer = 12 * "D" + "MM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 13 * "D" + "MM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer


def test_reorder_buffer_deletions_kabat_three_matches_on_buffer_right_end():
    # even buffer len
    relative_indel_position = 6
    buffer = 11 * "D" + "MMM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 12 * "D" + "MMM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions(relative_indel_position, matches, deletions) == buffer


def test_smooth_cdr_junctions_kabat_insertions_and_deletions_light_chain():
    regions = {
        "fwr1": "MMMIMMMMDMMMMIMMMMIMMMMMMMII",
        "cdr1": "MDMMMIMMIMMDM",
        "fwr2": "DMMIMMMMDMMMMMMMI",
        "cdr2": "MMIMMIIMMD",
        "fwr3": "IMIMMMMMMMMDMMMMDMMMMMMMMMIMMMMMMMIMI",
        "cdr3": "IMDDMDMDMD",
        "fwr4": "IMIMMMMIMMDMMI",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_LIGHT

    regions = {
        "fwr1": "MMMIMMMMDMMMMIMMMMIMMMMMMM",
        "cdr1": "MMMMIIMMMMMMM",
        "fwr2": "DMMIMMMMDMMMMMMM",
        "cdr2": "MMMIIIIMMMM",
        "fwr3": "MIMMMMMMMMDMMMMDMMMMMMMMMIMMMMMMMIM",
        "cdr3": "MMMMMDDMM",
        "fwr4": "MIMMMMIMMDMMI",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == KABAT_POSITIONS_LIGHT

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=ChainType.LIGHT, scheme=SCHEME)
        == then_alignment_str
    )


if __name__ == "__main__":
    test_smooth_cdr_junctions_kabat_no_indels()
    test_smooth_cdr_junctions_kabat_deletions_in_frameworks()
    test_smooth_cdr_junctions_kabat_deletions_in_cdrs()
    test_smooth_cdr_junctions_kabat_too_much_deletions()
    test_smooth_cdr_junctions_kabat_insertions_in_frameworks()
    test_smooth_cdr_junctions_kabat_insertions_in_cdrs()
    test_smooth_cdr_junctions_kabat_insertions_and_deletions()
    test_smooth_cdr_junctions_kabat_start_on_last_fwr1()
    test_smooth_cdr_junctions_kabat_start_on_cdr1_start()
    test_smooth_cdr_junctions_kabat_start_on_cdr1()
    test_smooth_cdr_junctions_kabat_start_on_cdr1_end()
    test_smooth_cdr_junctions_kabat_start_on_fwr2()
    test_smooth_cdr_junctions_kabat_start_on_cdr2_end()
    test_smooth_cdr_junctions_kabat_short_fwr4()
    test_smooth_cdr_junctions_kabat_end_on_fwr4_junction()
    test_smooth_cdr_junctions_kabat_end_in_fwr4_junction()
    test_smooth_cdr_junctions_kabat_no_fwr4()
    test_smooth_cdr_junctions_kabat_end_on_cdr3()
    test_smooth_cdr_junctions_kabat_end_on_cdr3_before_indel_pos()
    test_smooth_cdr_junctions_kabat_single_match_on_cdr3()
    test_smooth_cdr_junctions_kabat_end_on_fwr3_junction()
    test_smooth_cdr_junctions_kabat_end_on_fwr3_junction_last_missing()
    test_smooth_cdr_junctions_kabat_end_on_fwr3()
    test_smooth_cdr_junctions_kabat_end_on_cdr2()
    test_reorder_buffer_deletions_kabat_no_matches()
    test_reorder_buffer_deletions_kabat_single_match_at_buffer_end()
    test_reorder_buffer_deletions_kabat_match_on_each_buffer_end()
    test_reorder_buffer_deletions_kabat_two_matches_on_buffer_right_end()
    test_reorder_buffer_deletions_kabat_three_matches_on_buffer_right_end()
    test_smooth_cdr_junctions_kabat_insertions_and_deletions_light_chain()
