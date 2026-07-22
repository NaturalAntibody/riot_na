import pytest

from riot_na import Scheme, get_or_create_riot_aa
from riot_na.api.utils import (
    get_primary_seq,
    get_region_position_indices,
    get_regions_position_indices,
    get_scheme_residue_mapping_by_region,
    get_scheme_residue_mapping_insertions_as_letters,
    int_to_str_insertion,
    map_insertion_number_to_letter,
    scheme_positions_to_index,
    str_to_int_insertion,
)
from riot_na.data.model import AirrRearrangementEntryAA, ShortRegion

HEAVY_CHAIN_WITH_INSERTIONS = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYNMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLSLQMNSLRDEDTAVYYCARAYYQPGKPGKQQYGMDVWGQGTTVTVSS"


@pytest.fixture(scope="module", name="airr_heavy")
def heavy_chain_airr() -> AirrRearrangementEntryAA:
    return get_or_create_riot_aa().run_on_sequence("", HEAVY_CHAIN_WITH_INSERTIONS, scheme=Scheme.IMGT)


@pytest.mark.parametrize(
    ("n", "expected"),
    [
        (1, "A"),
        (2, "B"),
        (26, "Z"),
        (27, "AA"),
        (28, "AB"),
        (52, "AZ"),
        (53, "BA"),
        (100, "CV"),
    ],
)
def test_int_to_str_insertion(n, expected):
    assert int_to_str_insertion(n) == expected


def test_int_to_str_insertion_rejects_non_positive():
    with pytest.raises(ValueError, match="positive integer"):
        int_to_str_insertion(0)


@pytest.mark.parametrize(
    ("position", "expected"),
    [
        ("111", 111),
        ("111.1", 111),
        ("112.2", 112),
    ],
)
def test_str_to_int_insertion(position, expected):
    assert str_to_int_insertion(position) == expected


def test_map_insertion_number_to_letter_no_decimal():
    # given
    position = "111"

    # when
    result = map_insertion_number_to_letter(position)

    # then
    assert result == "111"


def test_map_insertion_number_to_letter_decimal():
    # given
    position = "111.1"

    # when
    result = map_insertion_number_to_letter(position)

    # then
    assert result == "111A"


def test_map_insertion_number_to_letter_z():
    # given
    position = "111.26"

    # when
    result = map_insertion_number_to_letter(position)

    # then
    assert result == "111Z"


def test_map_insertion_number_to_letter_double_letter():
    # given
    position = "111.27"

    # when
    result = map_insertion_number_to_letter(position)

    # then
    assert result == "111AA"


def test_map_insertion_number_to_letter_large_insertion():
    # given
    position = "111.100"

    # when
    result = map_insertion_number_to_letter(position)

    # then
    assert result == "111CV"


def test_get_region_position_indices(airr_heavy):
    assert get_region_position_indices(airr_heavy, ShortRegion.CDR1) == list(range(25, 33))
    assert get_region_position_indices(airr_heavy, ShortRegion.CDR3) == list(range(96, 115))


def test_get_regions_position_indices(airr_heavy):
    result = get_regions_position_indices(airr_heavy, [ShortRegion.CDR1, ShortRegion.CDR3])

    assert result == sorted(set(range(25, 33)) | set(range(96, 115)))


def test_scheme_positions_to_index(airr_heavy):
    assert scheme_positions_to_index(airr_heavy, ["27", "111.1", "111.2"]) == [25, 103, 104]


def test_get_primary_seq(airr_heavy):
    assert get_primary_seq(airr_heavy) == HEAVY_CHAIN_WITH_INSERTIONS


def test_get_scheme_residue_mapping_insertions_as_letters(airr_heavy):
    result = get_scheme_residue_mapping_insertions_as_letters(airr_heavy)

    assert result["111"] == "P"
    assert result["111A"] == "G"
    assert result["111B"] == "K"
    assert result["111C"] == "P"
    assert result["112A"] == "Q"
    assert result["112B"] == "K"
    assert result["112C"] == "G"
    assert result["112"] == "Q"


def test_get_scheme_residue_mapping_by_region(airr_heavy):
    result = get_scheme_residue_mapping_by_region(airr_heavy)

    assert len(result[ShortRegion.FW1]) == 25
    assert result[ShortRegion.CDR3] == {
        "105": "A",
        "106": "R",
        "107": "A",
        "108": "Y",
        "109": "Y",
        "110": "Q",
        "111": "P",
        "111.1": "G",
        "111.2": "K",
        "111.3": "P",
        "112.3": "G",
        "112.2": "K",
        "112.1": "Q",
        "112": "Q",
        "113": "Y",
        "114": "G",
        "115": "M",
        "116": "D",
        "117": "V",
    }
