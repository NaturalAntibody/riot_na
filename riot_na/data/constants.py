import json
from pathlib import Path

AMINO_ACIDS = set(list("QWERTYIPASDFGHKLCVNM"))


BLOSUM_62 = json.load((Path(__file__).parent / "blosum_62.json").open("r"))
