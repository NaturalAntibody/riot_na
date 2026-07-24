from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, TextIO, Union


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """Minimal FASTA record (description + sequence)."""

    description: str
    seq: str


def parse_fasta(path: Union[Path, str, TextIO]) -> Iterator[FastaRecord]:
    """Yield FASTA records from a file path or text stream.

    Header lines start with ``>``; the description is the remainder of that line
    (stripped). Sequence lines are concatenated with whitespace removed.
    """
    if isinstance(path, (str, Path)):
        with open(path) as handle:
            yield from _parse_fasta_handle(handle)
    else:
        yield from _parse_fasta_handle(path)


def _parse_fasta_handle(handle: TextIO) -> Iterator[FastaRecord]:
    description: str | None = None
    seq_parts: list[str] = []

    for line in handle:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if description is not None:
                yield FastaRecord(description=description, seq="".join(seq_parts))
            description = line[1:].strip()
            seq_parts = []
        else:
            if description is None:
                raise ValueError("FASTA file has sequence data before the first header line")
            seq_parts.append(line.strip())

    if description is not None:
        yield FastaRecord(description=description, seq="".join(seq_parts))
