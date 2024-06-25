from riot_na.data.model import AlignmentString, ChainType, Scheme
from riot_na.schemes.smooth_alignment import (
    _reorder_cdr_deletions_imgt,
    smooth_cdr_junctions,
)

SCHEME = Scheme.IMGT
CHAIN_TYPE = ChainType.HEAVY


def test_smooth_cdr_junctions_imgt_no_indels():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == given_alignment_str
    )


def test_smooth_cdr_junctions_deletions_not_in_cdr():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMDMD",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "DMMMMMMMMMMDMMMMD",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDD",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "DMMMMDMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMDMD",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "DMMMMMMMMMMDMMMMD",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDD",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "DMMMMDMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_deletions_in_cdrs():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMDMD",
        "cdr1": "MMDMMMMMMDMM",
        "fwr2": "DMMMMMMMMMMDMMMMD",
        "cdr2": "MMMMMMMMMD",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDD",
        "cdr3": "MMMMMMMMMMDMM",
        "fwr4": "DMMMMDMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMDMD",
        "cdr1": "MMMMMDDMMMMM",
        "fwr2": "DMMMMMMMMMMDMMMMD",
        "cdr2": "MMMMMDMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDD",
        "cdr3": "MMMMMMDMMMMMM",
        "fwr4": "DMMMMDMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_insertions_in_frameworks():
    regions = {
        "fwr1": "MIIMMMMMMMMMMMMMMIMMMMMMMMMIMIM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MIIIMMMMMMIMMMMIMMMMMIIM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MIIMMMMMMMMMMMMMIIMMMMMMMMMMMMMIIMMMMMMMMMMMIIM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == given_alignment_str
    )


def test_smooth_cdr_junctions_imgt_insertions_in_cdrs():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMMII",
        "cdr1": "MIMMMMMMMMIMMM",
        "fwr2": "IIMMMMMMMMMMMMMMMMM",
        "cdr2": "MMIMMMMMMMMI",
        "fwr3": "IMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "IMMIMMMMMMMMMIIMM",
        "fwr4": "IIMMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMIIIIIIMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMIIIMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMIIIIIIMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_insertions_and_deletions():
    regions = {
        "fwr1": "MMIMDMMIMMMMMMMMMMMDMMMMMMMD",
        "cdr1": "MIMMMDMMMMMIMM",
        "fwr2": "IMIMMMMMMMMMMMMMMMM",
        "cdr2": "IMMMDMMMMMIIMI",
        "fwr3": "DMMMMMMMMDMMMMMMMMMMIMMMMMDMMMMMIMMMMMMDD",
        "cdr3": "DMIIMDMMDDMMIMDMI",
        "fwr4": "MIMMMMMIMDMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "MMIMDMMIMMMMMMMMMMMDMMMMMMMD",
        "cdr1": "MMMMMMIIMMMMMM",
        "fwr2": "MIMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMIIIMMMMM",
        "fwr3": "DMMMMMMMMDMMMMMMMMMMIMMMMMDMMMMMIMMMMMMDD",
        "cdr3": "MMMMMMDMMMMMM",
        "fwr4": "MIMMMMMIMDMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_start_on_last_fwr1():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_start_on_cdr1_start():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_start_on_cdr1():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNNMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MMMMDDDDDMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_start_on_cdr1_end():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNNNNNNNNM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "MDDDDDDDDDDD",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_start_on_fwr2():
    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNNNNNNNNN",
        "fwr2": "MDMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "NNNNNNNNNNNNNNNNNNNNNNNNNN",
        "cdr1": "NNNNNNNNNNNN",
        "fwr2": "MDMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMMMMMMMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_short_fwr4():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 121

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMMMMMMMM",
        "fwr4": "MMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 121

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_on_fwr4_junction():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMDMMMMMMMMMM",
        "fwr4": "DM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 119

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMDMMMMMM",
        "fwr4": "DM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 119

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_in_fwr4_junction():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMDMMMMMMMMMM",
        "fwr4": "M",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 118

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMDMMMMMM",
        "fwr4": "M",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 118

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_no_fwr4():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMDMMMMMMMMMM",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 117

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMMDMMMMMM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 117

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


if __name__ == "__main__":
    test_smooth_cdr_junctions_imgt_no_fwr4()


def test_smooth_cdr_junctions_imgt_end_on_cdr3():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMDMMMMMMM",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 114

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMMMMDDDDMMMM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 117

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_on_cdr3_before_indel_pos():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMDM",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 108

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MMDDDDDDDDDDM",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 117

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_single_match_on_cdr3():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "M",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 105

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "MDDDDDDDDDDDD",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 117

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_on_fwr3_junction():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 104

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 104

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_on_fwr3_junction_last_missing():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 103

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 103

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_on_fwr3():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMDMMM",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 98

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMMMM",
        "fwr3": "MMMMMMMMMMMMMMMMMMMMMMMMMMMMMDMMM",
        "cdr3": "",
        "fwr4": "",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 98

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=CHAIN_TYPE, scheme=SCHEME)
        == then_alignment_str
    )


def test_smooth_cdr_junctions_imgt_end_on_cdr2():
    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMMMDM",
        "fwr3": "",
        "cdr3": "",
        "fwr4": "",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 64

    regions = {
        "fwr1": "MMMMMMMMMMMMMMMMMMMMMMMMMM",
        "cdr1": "MMMMMMMMMMMM",
        "fwr2": "MMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMDDMMMM",
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


def test_reorder_buffer_deletions_imgt_no_matches():
    # even buffer len
    relative_indel_position = 6
    buffer = 14 * "D"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 15 * "D"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == buffer


def test_reorder_buffer_deletions_imgt_single_match_on_right_junction():
    # even buffer len
    relative_indel_position = 6
    buffer = 13 * "D" + "M"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 14 * "D" + "M"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == buffer


if __name__ == "__main__":
    test_reorder_buffer_deletions_imgt_single_match_on_right_junction()


def test_reorder_buffer_deletions_imgt_two_matches_on_right_junction():
    # even buffer len
    relative_indel_position = 6
    buffer = 12 * "D" + "MM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == buffer

    # odd buffer len
    relative_indel_position = 7
    buffer = 13 * "D" + "MM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == buffer


def test_reorder_buffer_deletions_imgt_three_matches_on_right_junction():
    # even buffer len
    relative_indel_position = 6
    buffer = 11 * "D" + "MMM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == "M" + 11 * "D" + "MM"

    # odd buffer len
    relative_indel_position = 7
    buffer = 12 * "D" + "MMM"

    matches = buffer.count("M")
    deletions = buffer.count("D")

    assert _reorder_cdr_deletions_imgt(relative_indel_position, matches, deletions) == "M" + 12 * "D" + "MM"


def test_smooth_cdr_junctions_imgt_insertions_and_deletions_light_chain():
    regions = {
        "fwr1": "MMIMDMMIMMMMMMMMMMMDMMMMMMMD",
        "cdr1": "MIMMMDMMMMMIMM",
        "fwr2": "IMIMMMMMMMMMMMMMMMM",
        "cdr2": "IMMMDMMMMMIIMI",
        "fwr3": "DMMMMMMMMDMMMMMMMMMMIMMMMMDMMMMMIMMMMMMDD",
        "cdr3": "DMIIMDMMDDMMIMDMI",
        "fwr4": "MIMMMMMIMDMMM",
    }
    given_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(given_alignment_str.count, ["M", "D", "N"])) == 128

    regions = {
        "fwr1": "MMIMDMMIMMMMMMMMMMMDMMMMMMMD",
        "cdr1": "MMMMMMIIMMMMMM",
        "fwr2": "MIMMMMMMMMMMMMMMMM",
        "cdr2": "MMMMMIIIMMMMM",
        "fwr3": "DMMMMMMMMDMMMMMMMMMMIMMMMMDMMMMMIMMMMMMDD",
        "cdr3": "MMMMMMDMMMMMM",
        "fwr4": "MIMMMMMIMDMMM",
    }
    then_alignment_str = AlignmentString("".join(regions.values()))
    assert sum(map(then_alignment_str.count, ["M", "D", "N"])) == 128

    assert (
        smooth_cdr_junctions(alignment_str=given_alignment_str, chain_type=ChainType.LIGHT, scheme=SCHEME)
        == then_alignment_str
    )


if __name__ == "__main__":
    test_smooth_cdr_junctions_imgt_no_indels()
    test_smooth_cdr_junctions_imgt_insertions_and_deletions()
    test_smooth_cdr_junctions_imgt_start_on_last_fwr1()
    test_smooth_cdr_junctions_imgt_start_on_cdr1_start()
    test_smooth_cdr_junctions_imgt_start_on_cdr1()
    test_smooth_cdr_junctions_imgt_start_on_cdr1_end()
    test_smooth_cdr_junctions_imgt_start_on_fwr2()
    test_smooth_cdr_junctions_imgt_short_fwr4()
    test_smooth_cdr_junctions_imgt_end_on_fwr4_junction()
    test_smooth_cdr_junctions_imgt_end_in_fwr4_junction()
    test_smooth_cdr_junctions_imgt_no_fwr4()
    test_smooth_cdr_junctions_imgt_end_on_cdr3()
    test_smooth_cdr_junctions_imgt_end_on_cdr3_before_indel_pos()
    test_smooth_cdr_junctions_imgt_single_match_on_cdr3()
    test_smooth_cdr_junctions_imgt_end_on_fwr3_junction()
    test_smooth_cdr_junctions_imgt_end_on_fwr3_junction_last_missing()
    test_smooth_cdr_junctions_imgt_end_on_fwr3()
    test_smooth_cdr_junctions_imgt_end_on_cdr2()
    test_smooth_cdr_junctions_imgt_insertions_and_deletions_light_chain()
