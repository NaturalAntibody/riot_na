import math

import pytest

from riot_na.airr.float_formatting import format_evalue_string_scientific


@pytest.mark.parametrize(
    "log10_e_value, expected",
    [
        (math.log10(2.26e-141), "2.26e-141"),
        (math.log10(9.02e-18), "9.02e-18"),
        (math.log10(24.1875), "2.42e1"),
        (0.0, "1.00e0"),
        (-400.0, "1.00e-400"),
        (-5.0000001, "1.00e-5"),
        (-4.99999999, "1.00e-5"),
        (None, None),
    ],
)
def test_format_evalue_string_scientific(log10_e_value, expected):
    assert format_evalue_string_scientific(log10_e_value) == expected
