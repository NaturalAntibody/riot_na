import math
from typing import Optional


def format_evalue_string_scientific(log10_e_value: Optional[float]) -> Optional[str]:
    """Convert a log10 E-value back to a scientific notation string, e.g. -140.645 -> '2.26e-141'."""
    if log10_e_value is None:
        return None
    exponent = math.floor(log10_e_value)
    mantissa = round(10 ** (log10_e_value - exponent), 2)
    if mantissa >= 10:
        mantissa /= 10
        exponent += 1
    return f"{mantissa:.2f}e{exponent}"
