import math
from inspect import getframeinfo

import pandas as pd  # type: ignore

from riot_na.airr.airr_validator import REGIONS
from riot_na.data.model import AirrRearrangementEntryNT


def get_sequence_alignment_with_masked_regions(record: pd.Series) -> str:
    def _get_mask_char(region_name: str) -> str:
        return "_" if "fwr" in region_name else "-"

    seq = record["sequence_alignment"]
    regions = record[REGIONS].replace(math.nan, None).to_dict().items()  # type: ignore
    for region_name, region in regions:
        if region is not None and len(region) > 0:
            seq = seq.replace(region, _get_mask_char(region_name) * len(region))
    return seq


def get_asserion_error_msg(currentframe, rearrangement: AirrRearrangementEntryNT) -> str:
    frameinfo = getframeinfo(currentframe)
    return (
        f"{frameinfo.filename}:{frameinfo.lineno} "
        f"sequence_id: {rearrangement.sequence_header.split(' ', maxsplit=1)[0]}"
    )
