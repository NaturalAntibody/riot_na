import csv
from dataclasses import fields
from typing import Any, Dict, Iterable, Optional, Type

from riot_na.data.model import AirrRearrangementEntry_co, serialize_airr_entry


class AirrRearrangementEntryWriter:
    def __init__(
        self,
        file_handle: Any,
        cls: Type[AirrRearrangementEntry_co],
        dialect: str = "excel",
        fmtparams: Optional[Dict[str, Any]] = None,
    ):
        if not file_handle:
            raise ValueError("The f argument is required")

        if fmtparams is None:
            fmtparams = {}

        self._fieldnames = [x.name for x in fields(cls)]
        self._writer = csv.DictWriter(file_handle, dialect=dialect, fieldnames=self._fieldnames, **fmtparams)

    def write(
        self,
        data: Iterable[AirrRearrangementEntry_co] | Iterable[list[AirrRearrangementEntry_co]],
        skip_header: bool = False,
    ):
        if not skip_header:
            self._writer.writeheader()

        for item in data:
            if isinstance(item, list):
                for subitem in item:
                    row = serialize_airr_entry(subitem)
                    self._writer.writerow(row)
            else:
                row = serialize_airr_entry(item)
                self._writer.writerow(row)
