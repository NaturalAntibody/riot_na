from riot_na.api.utils import map_insertion_number_to_letter


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
