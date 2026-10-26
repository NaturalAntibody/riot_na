import csv
import json
from collections import defaultdict
from dataclasses import asdict
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence, Type

from riot_na.api.molecule_annotation import (
    BINDING_SITE_FIELDS,
    LAYOUT_FIELDS,
    VDomainPair,
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


V_DOMAIN_PAIR_FIELDS = ("molecule_id", "chain_1", "v_domain_1", "chain_2", "v_domain_2")


def read_v_domain_pairs(input_file_path: Path) -> dict[str, list[VDomainPair]]:
    """
    V domain pairs of molecules (pairs of annotate_molecule) by molecule ID, from a CSV file (TSV, if the file
    extension is .tsv) with columns molecule_id, chain_1, v_domain_1, chain_2, v_domain_2: one row per binding site,
    chains and V domains by number (1-based position of the chain in the molecule chains and of the V domain among the
    V domains of the chain).
    """
    delimiter = "\t" if input_file_path.suffix.lower() == ".tsv" else ","
    v_domain_pairs: dict[str, list[VDomainPair]] = defaultdict(list)
    with open(input_file_path, newline="") as input_file:
        reader = csv.DictReader(input_file, delimiter=delimiter)
        missing_fields = [field for field in V_DOMAIN_PAIR_FIELDS if field not in (reader.fieldnames or ())]
        if missing_fields:
            raise ValueError(f"Missing columns in {input_file_path}: {', '.join(missing_fields)}")
        for row in reader:
            try:
                chain_1, v_domain_1, chain_2, v_domain_2 = (int(row[field]) for field in V_DOMAIN_PAIR_FIELDS[1:])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"Invalid V domain pair in line {reader.line_num} of {input_file_path}, "
                    "chains and V domains must be numbers"
                ) from error
            v_domain_pairs[row["molecule_id"]].append((chain_1, v_domain_1, chain_2, v_domain_2))
    return dict(v_domain_pairs)


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
