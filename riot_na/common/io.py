import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence, Type

from riot_na.api.molecule_annotation import (
    BINDING_SITE_FIELDS,
    LAYOUT_FIELDS,
    binding_site_rows,
    molecule_layout_rows,
)
from riot_na.common.airr_csv_writer import AirrRearrangementEntryWriter
from riot_na.common.fasta import parse_fasta
from riot_na.data.model import AirrRearrangementEntry_co  # type: ignore
from riot_na.data.model import MoleculeAnnotation


def read_fasta(path: Path) -> tuple[list[str], list[str]]:
    headers = []
    sequences = []
    for record in parse_fasta(path):
        headers.append(record.description)
        sequences.append(record.seq)
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


def write_molecule_layout_csv(output_file_path: Path, molecules: Iterable[MoleculeAnnotation]):
    """Domains of all chains of the molecules, one row per domain, with binding sites of V domains."""
    _write_molecule_rows(output_file_path, molecules, LAYOUT_FIELDS, molecule_layout_rows)


def write_binding_sites_csv(output_file_path: Path, molecules: Iterable[MoleculeAnnotation]):
    """Binding sites of the molecules, one row per binding site."""
    _write_molecule_rows(output_file_path, molecules, BINDING_SITE_FIELDS, binding_site_rows)


def write_molecules_json(output_file_path: Path, molecules: Iterable[MoleculeAnnotation]):
    """Molecules with all details: chain annotations (rearrangements, numbered constant domains) and binding sites."""
    output_file_path.parent.mkdir(exist_ok=True, parents=True)
    with open(output_file_path, "w") as output:
        json.dump([asdict(molecule) for molecule in molecules], output, indent=2)


def _write_molecule_rows(
    output_file_path: Path,
    molecules: Iterable[MoleculeAnnotation],
    fields: Sequence[str],
    rows: Callable[[MoleculeAnnotation], list[dict[str, Any]]],
):
    output_file_path.parent.mkdir(exist_ok=True, parents=True)
    with open(output_file_path, "w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for molecule in molecules:
            writer.writerows(rows(molecule))
