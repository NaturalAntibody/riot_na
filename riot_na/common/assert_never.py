from typing import NoReturn


def assert_never(value: NoReturn) -> NoReturn:
    raise ValueError(f"Value {value} not handled")
