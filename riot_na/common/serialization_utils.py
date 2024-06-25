import base64


def base64_encode(val: str) -> str:
    return base64.b64encode(val.encode("utf-8")).decode("utf-8")


def base64_decode(val: str) -> str:
    return base64.b64decode(val.encode("utf-8")).decode("utf-8")
