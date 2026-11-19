"""
Build EU, Kabat and IMGT scheme mappings of human IGHG1-4, IGKC and IGLC C germlines
(riot_na/databases/scheme_mappings/human/<c_scheme>/c_scheme_mapping.csv).

Source: IMGT Scientific chart, correspondence between C numberings, extracted to imgt_c_numbering.tsv:
- https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGHGnber.html - CH1, CH2 and CH3 of IGHG1 and the hinges
  of all four IGHG genes,
- https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGKCnber.html - IGKC,
- https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGLCnber.html - IGLC1, IMGT and Kabat only, EU numbering
  is not defined for lambda chains.

The tables number one gene of each isotype. CH domains of IGHG2-4 have the same lengths as in IGHG1 and IGLC2-7
have the same length as IGLC1, so they take positions of these genes. The exception is IGHG2 CH2, which lacks one
lower hinge residue: EU and Kabat place the gap at EU 236 (IGHG2 P233, V234, A235, no G236), as conventionally used
in Fc engineering; in IMGT the residues before CH2 position 1 are insertions numbered towards it (1.5-1.1).

A scheme mapping is a cigar over the scheme positions (riot_na.data.scheme_definitions):
- M - residue on a scheme position,
- D - scheme position without a residue,
- I - insertion: residue on an IMGT position with a dot (e.g. 1.4), on a lettered Kabat position (e.g. 107A, 241A)
  or without EU position (IGHG3 hinge). Missing insertion positions are skipped.

Germlines in the database are membrane-bound forms, so CHS (EU 446-447) is numbered only when the germline has it.

Run from the repository root:
    poetry run python -m notebooks.data_processing.c_scheme_mappings
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from string import ascii_uppercase
from typing import Optional

from notebooks.data_processing.c_scheme_mapping_validators import (
    validate_numbered_germline,
    validate_scheme_cigar,
)
from riot_na.common.fasta import parse_fasta
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import ChainType, CScheme, Organism
from riot_na.data.scheme_definitions import get_c_scheme_positions

IMGT_TABLE_PATH = Path(__file__).parent / "imgt_c_numbering.tsv"
# germline - gene whose positions are given in the IMGT tables, hinges of IGHG genes are given for each gene
IMGT_TABLE_GENE_BY_GERMLINE = {
    "IGHG1": "IGHG1",
    "IGHG2": "IGHG1",
    "IGHG3": "IGHG1",
    "IGHG4": "IGHG1",
    "IGKC": "IGKC",
    "IGLC1": "IGLC1",
    "IGLC2": "IGLC1",
    "IGLC3": "IGLC1",
    "IGLC6": "IGLC1",
    "IGLC7": "IGLC1",
}
CHS_IMGT_POSITIONS = ("129", "130")
# IMGT insertions 1.x and 85.x precede their position, the others follow it
IMGT_POSITIONS_WITH_PRECEDING_INSERTIONS = ("1", "85")


@dataclass(frozen=True)
class ImgtCorrespondenceRow:
    """Residue of an IMGT table of correspondence between C numberings with its positions in the schemes."""

    gene: str
    # CH1, H (hinge), H1-H4 (IGHG3 hinge exons), CH2, CH3 or CL
    region: str
    imgt_position: str
    residue: str
    # "-" if the residue has no EU position
    eu: str
    kabat: str


@dataclass(frozen=True)
class NumberedResidue:
    residue: str
    # positions as written in the IMGT tables, None if the residue has no position in a scheme
    positions: dict[CScheme, Optional[str]]


def read_imgt_table() -> list[ImgtCorrespondenceRow]:
    with IMGT_TABLE_PATH.open() as table_file:
        return [ImgtCorrespondenceRow(**row) for row in csv.DictReader(table_file, delimiter="\t")]


def read_germlines(path: Path) -> dict[str, str]:
    return {record.description.split("\t")[0]: record.seq for record in parse_fasta(path)}


def scheme_positions(row: ImgtCorrespondenceRow, imgt_position: str) -> dict[CScheme, Optional[str]]:
    return {CScheme.EU: row.eu if row.eu != "-" else None, CScheme.KABAT: row.kabat, CScheme.IMGT: imgt_position}


def number_domain(germline: str, rows: list[ImgtCorrespondenceRow], domain: str) -> list[NumberedResidue]:
    if domain == "CH3":
        chs_rows = [row for row in rows if row.imgt_position in CHS_IMGT_POSITIONS]
        rows = [row for row in rows if row.imgt_position not in CHS_IMGT_POSITIONS]
        chs = germline[len(rows) : len(rows) + len(chs_rows)]
        if chs == "".join(row.residue for row in chs_rows):
            rows = rows + chs_rows

    assert len(germline) >= len(rows), f"{domain} is longer than the germline"
    return [
        NumberedResidue(residue, scheme_positions(row, f"{domain}_{row.imgt_position}"))
        for residue, row in zip(germline, rows)
    ]


def number_ighg2_ch2(germline: str, rows: list[ImgtCorrespondenceRow]) -> list[NumberedResidue]:
    imgt_positions = [f"CH2_{row.imgt_position}" for row in rows if row.imgt_position != "1.6"]
    eu_rows = [row for row in rows if row.eu != "236"]
    return [
        NumberedResidue(residue, scheme_positions(eu_row, imgt_position))
        for residue, imgt_position, eu_row in zip(germline, imgt_positions, eu_rows)
    ]


def number_heavy_germline(gene: str, germline: str, table: list[ImgtCorrespondenceRow]) -> list[NumberedResidue]:
    domain_rows = {
        domain: [row for row in table if row.gene == IMGT_TABLE_GENE_BY_GERMLINE[gene] and row.region == domain]
        for domain in ("CH1", "CH2", "CH3")
    }
    hinge_rows = [row for row in table if row.gene == gene and row.region.startswith("H")]

    numbered: list[NumberedResidue] = []
    numbered += number_domain(germline, domain_rows["CH1"], "CH1")

    hinge = germline[len(numbered) : len(numbered) + len(hinge_rows)]
    assert hinge == "".join(row.residue for row in hinge_rows), f"{gene} hinge differs from IMGT: {hinge}"
    numbered += [
        NumberedResidue(row.residue, scheme_positions(row, f"{row.region}_{row.imgt_position}")) for row in hinge_rows
    ]

    if gene == "IGHG2":
        numbered += number_ighg2_ch2(germline[len(numbered) :], domain_rows["CH2"])
    else:
        numbered += number_domain(germline[len(numbered) :], domain_rows["CH2"], "CH2")

    numbered += number_domain(germline[len(numbered) :], domain_rows["CH3"], "CH3")
    return numbered


def number_light_germline(gene: str, germline: str, table: list[ImgtCorrespondenceRow]) -> list[NumberedResidue]:
    rows = [row for row in table if row.gene == IMGT_TABLE_GENE_BY_GERMLINE[gene]]
    assert len(germline) == len(rows), f"{gene} length differs from {IMGT_TABLE_GENE_BY_GERMLINE[gene]}"
    return number_domain(germline, rows, "CL")


def insertion_anchor(position: Optional[str], position_indices: dict[str, int], last_index: int) -> int:
    """
    Index of the scheme position the insertion is attached to: the position it follows, or for IMGT 1.x and 85.x
    the position before the one they precede.
    """
    if position is None:
        return last_index

    if "_" in position:
        domain, imgt_position = position.split("_")
        number = imgt_position.split(".")[0]
        index = position_indices[f"{domain}_{number}"]
        return index - 1 if number in IMGT_POSITIONS_WITH_PRECEDING_INSERTIONS else index

    return position_indices[position.rstrip(ascii_uppercase)]


def scheme_cigar(numbered: list[NumberedResidue], c_scheme: CScheme, chain_type: ChainType) -> str:
    position_indices = {position: i for i, position in enumerate(get_c_scheme_positions(c_scheme, chain_type))}

    cigar = ""
    last_index = -1
    for residue in numbered:
        position = residue.positions[c_scheme]

        if position in position_indices:
            index = position_indices[position]
            assert index > last_index, f"Position {position} out of order"
            cigar += "D" * (index - last_index - 1) + "M"
        else:
            index = insertion_anchor(position, position_indices, last_index)
            assert index >= last_index, f"Insertion {position} out of order"
            cigar += "D" * (index - last_index) + "I"

        last_index = index

    return cigar


def riot_position(position: Optional[str], c_scheme: CScheme) -> Optional[str]:
    """Position as numbered by riot: lettered Kabat insertions A-Z, AA, BB, ... are written as .1-.26, .27, ..."""
    if c_scheme != CScheme.KABAT or position is None or position.isdigit():
        return position

    number = position.rstrip(ascii_uppercase)
    letters = position[len(number) :]
    return f"{number}.{ascii_uppercase.index(letters[0]) + 1 + 26 * (len(letters) - 1)}"


def main():
    table = read_imgt_table()

    scheme_mappings: dict[CScheme, dict[str, str]] = {c_scheme: {} for c_scheme in CScheme}
    for gene in IMGT_TABLE_GENE_BY_GERMLINE:
        locus = gene[:3].lower()
        germline = read_germlines(GENE_DB_DIR / "gene_db" / "aa_genes" / "c_genes" / "human" / f"{locus}.fasta")[gene]
        deduplicated_germline = read_germlines(
            GENE_DB_DIR / "gene_db" / "aa_genes_deduplicated" / "c_genes" / "human" / f"{locus}.fasta"
        )[gene]
        assert germline == deduplicated_germline, f"{gene} differs between germline databases"

        if locus == "igh":
            chain_type = ChainType.HEAVY
            numbered = number_heavy_germline(gene, germline, table)
        else:
            chain_type = ChainType.LIGHT
            numbered = number_light_germline(gene, germline, table)

        imgt_table = {f"{row.region}_{row.imgt_position}": row.residue for row in table if row.gene == gene}
        validate_numbered_germline(
            gene,
            germline,
            [(residue.residue, residue.positions[CScheme.IMGT]) for residue in numbered],
            imgt_table,
        )

        for c_scheme in CScheme:
            positions = [riot_position(residue.positions[c_scheme], c_scheme) for residue in numbered]
            # EU numbering is not defined for IGLC
            if not any(positions):
                continue

            cigar = scheme_cigar(numbered, c_scheme, chain_type)
            validate_scheme_cigar(
                gene, "".join(residue.residue for residue in numbered), positions, cigar, c_scheme, chain_type
            )
            scheme_mappings[c_scheme][gene] = cigar

        print(f"{gene}: {len(numbered)} numbered residues, unnumbered tail: {germline[len(numbered):][:10]}")

    for c_scheme, gene_cigars in scheme_mappings.items():
        output_path = GENE_DB_DIR / "scheme_mappings" / Organism.HOMO_SAPIENS.value / c_scheme.value
        output_path.mkdir(parents=True, exist_ok=True)

        with (output_path / "c_scheme_mapping.csv").open("w") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=["gene_id", "scheme_cigar"])
            writer.writeheader()
            for gene, cigar in gene_cigars.items():
                writer.writerow({"gene_id": gene, "scheme_cigar": cigar})


if __name__ == "__main__":
    main()
