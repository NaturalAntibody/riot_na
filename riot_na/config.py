import json
import os
from pathlib import Path

from dotenv import load_dotenv

if bool(json.loads(os.getenv("LOAD_DOTENV", "true"))):
    # this is somewhat a hack - on production code is run directly from a .zip file. Because of that dotenv is not able
    # to find an .env file and will crash. Whole env is being set
    try:
        load_dotenv()
    except OSError:
        pass

GENE_DB_DIR = Path(os.getenv("GENE_DB_DIR", Path(__file__).parent / "databases"))
