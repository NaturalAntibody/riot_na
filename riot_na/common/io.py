from pathlib import Path
from typing import Iterable, Type

from Bio import SeqIO

from riot_na.common.airr_csv_writer import AirrRearrangementEntryWriter
from riot_na.data.model import AirrRearrangementEntry_co  # type: ignore


def read_fasta(path: Path) -> tuple[list[str], list[str]]:
    headers = []
    sequences = []
    for record in SeqIO.parse(path, "fasta"):
        headers.append(record.description)
        sequences.append(str(record.seq))
    return headers, sequences


def count_fasta_records(input_fasta_path: Path, input_format: str = "fasta") -> int:
    total = 0
    start_token = ">" if input_format == "fasta" else "@"
    with open(input_fasta_path) as input_file:
        for line in input_file:
            if line.startswith(start_token):
                total += 1
    return total


def write_airr_iter_to_csv(
    output_file_path: Path,
    cls: Type[AirrRearrangementEntry_co],
    airr_iter: Iterable[AirrRearrangementEntry_co] | Iterable[list[AirrRearrangementEntry_co]],
):
    output_file_path.parent.mkdir(exist_ok=True, parents=True)
    with open(output_file_path, "w") as output:
        AirrRearrangementEntryWriter(output, cls).write(airr_iter)
