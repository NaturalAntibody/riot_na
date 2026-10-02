"""
Build EU, Kabat and IMGT scheme mappings of human IGHG1-4, IGKC and IGLC C germlines.

Source: IMGT Scientific chart, correspondence between C numberings, extracted to imgt_c_numbering.tsv:
- https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGHGnber.html - CH1, CH2 and CH3 of IGHG1 and the hinges
  of all four IGHG genes,
- https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGKCnber.html - IGKC,
- https://www.imgt.org/IMGTScientificChart/Numbering/Hu_IGLCnber.html - IGLC1, IMGT and Kabat only, EU numbering
  is not defined for lambda chains.

CH domains of IGHG2-4 have the same lengths as in IGHG1 and take its positions. The exception is IGHG2 CH2,
which lacks one lower hinge residue: IMGT places the gap at CH2 1.4 (IMGT protein display of human IGHC),
EU and Kabat at EU 236 (IGHG2 P233, V234, A235, no G236), as conventionally used in Fc engineering.
IGLC2, IGLC3, IGLC6 and IGLC7 have the same length as IGLC1 and take its positions.

Residues without a position in a scheme (IGHG3 hinge in EU, Kabat 241A-241SS) become insertions after the
preceding position. Kabat lettered positions are written like riot Fv insertions (241A - 241.1), lambda C starts
at Kabat 107A (107.1). Germlines in the database are membrane-bound forms, so CHS (EU 446-447) is numbered only
when the germline has it.

Run from the repository root:
    poetry run python notebooks/data_processing/c_scheme_mappings.py
"""

import csv
from dataclasses import dataclass
from pathlib import Path
from string import ascii_uppercase
from typing import Optional

from riot_na.common.fasta import parse_fasta
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import ChainType, CScheme, Organism
from riot_na.data.scheme_definitions import get_c_scheme_positions

IMGT_TABLE_PATH = Path(__file__).parent / "imgt_c_numbering.tsv"
HEAVY_GENES = ("IGHG1", "IGHG2", "IGHG3", "IGHG4")
GENES_BY_LOCUS = {
    "igh": HEAVY_GENES,
    "igk": ("IGKC",),
    "igl": ("IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGLC7"),
}
# light chain gene - gene numbered in the IMGT table
LIGHT_GENES = {
    "IGKC": "IGKC",
    "IGLC1": "IGLC1",
    "IGLC2": "IGLC1",
    "IGLC3": "IGLC1",
    "IGLC6": "IGLC1",
    "IGLC7": "IGLC1",
}
CHS_IMGT_POSITIONS = ("129", "130")


@dataclass(frozen=True)
class TableRow:
    gene: str
    region: str
    imgt_position: str
    residue: str
    eu: str
    kabat: str


@dataclass(frozen=True)
class NumberedResidue:
    residue: str
    # None - no position in a scheme, the residue is an insertion after the preceding one
    positions: dict[CScheme, Optional[str]]


def read_imgt_table() -> list[TableRow]:
    with IMGT_TABLE_PATH.open() as table_file:
        return [TableRow(**row) for row in csv.DictReader(table_file, delimiter="\t")]


def read_germlines(path: Path) -> dict[str, str]:
    return {record.description.split("\t")[0]: record.seq for record in parse_fasta(path)}


def kabat_position(kabat: str) -> str:
    # lettered Kabat positions are insertions: A-Z, then AA, BB, ... - written as .1-.26, .27, ...
    number = kabat.rstrip(ascii_uppercase)
    letters = kabat[len(number) :]
    if not letters:
        return kabat
    return f"{number}.{ascii_uppercase.index(letters[0]) + 1 + 26 * (len(letters) - 1)}"


def scheme_positions(row: TableRow, imgt_position: str) -> dict[CScheme, Optional[str]]:
    eu = row.eu if row.eu != "-" else None
    return {CScheme.EU: eu, CScheme.KABAT: kabat_position(row.kabat), CScheme.IMGT: imgt_position}


def number_domain(germline: str, rows: list[TableRow], domain: str) -> list[NumberedResidue]:
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


def number_ighg2_ch2(germline: str, rows: list[TableRow]) -> list[NumberedResidue]:
    imgt_positions = [f"CH2_{row.imgt_position}" for row in rows if row.imgt_position != "1.4"]
    eu_rows = [row for row in rows if row.eu != "236"]
    return [
        NumberedResidue(residue, scheme_positions(eu_row, imgt_position))
        for residue, imgt_position, eu_row in zip(germline, imgt_positions, eu_rows)
    ]


def number_heavy_germline(gene: str, germline: str, table: list[TableRow]) -> list[NumberedResidue]:
    ighg1_domains = {
        domain: [row for row in table if row.gene == "IGHG1" and row.region == domain]
        for domain in ("CH1", "CH2", "CH3")
    }
    hinge_rows = [row for row in table if row.gene == gene and row.region.startswith("H")]

    numbered: list[NumberedResidue] = []
    numbered += number_domain(germline, ighg1_domains["CH1"], "CH1")

    hinge = germline[len(numbered) : len(numbered) + len(hinge_rows)]
    assert hinge == "".join(row.residue for row in hinge_rows), f"{gene} hinge differs from IMGT: {hinge}"
    numbered += [
        NumberedResidue(row.residue, scheme_positions(row, f"{row.region}_{row.imgt_position}")) for row in hinge_rows
    ]

    if gene == "IGHG2":
        numbered += number_ighg2_ch2(germline[len(numbered) :], ighg1_domains["CH2"])
    else:
        numbered += number_domain(germline[len(numbered) :], ighg1_domains["CH2"], "CH2")

    numbered += number_domain(germline[len(numbered) :], ighg1_domains["CH3"], "CH3")
    return numbered


def number_light_germline(gene: str, germline: str, table: list[TableRow]) -> list[NumberedResidue]:
    rows = [row for row in table if row.gene == LIGHT_GENES[gene]]
    assert len(germline) == len(rows), f"{gene} length differs from {LIGHT_GENES[gene]}"
    return number_domain(germline, rows, "CL")


def validate(gene: str, germline: str, numbered: list[NumberedResidue], table: list[TableRow]):
    assert germline.startswith("".join(residue.residue for residue in numbered))

    # every residue has an IMGT position
    imgt_numbering = {residue.positions[CScheme.IMGT]: residue.residue for residue in numbered}
    assert None not in imgt_numbering and len(imgt_numbering) == len(numbered)

    if gene in ("IGHG1", "IGKC", "IGLC1"):
        imgt_table = {f"{row.region}_{row.imgt_position}": row.residue for row in table if row.gene == gene}
        assert all(imgt_table[position] == residue for position, residue in imgt_numbering.items() if position)

    # conserved cysteines of C domains
    domains = ("CH1", "CH2", "CH3") if gene in HEAVY_GENES else ("CL",)
    for domain in domains:
        for position in ("23", "104"):
            assert imgt_numbering[f"{domain}_{position}"] == "C", f"{gene} {domain} {position}"


def scheme_cigar(numbered: list[NumberedResidue], c_scheme: CScheme, chain_type: ChainType) -> str:
    scheme_positions_list = get_c_scheme_positions(c_scheme, chain_type)
    position_indices = {position: i for i, position in enumerate(scheme_positions_list)}

    cigar = ""
    last_index = -1
    insertions = 0
    for residue in numbered:
        position = residue.positions[c_scheme]
        if position not in position_indices:
            assert last_index >= 0, "Insertion before the first scheme position"
            insertions += 1
            # lettered Kabat positions have to match the insertion labels produced by numbering
            assert position is None or position == f"{scheme_positions_list[last_index]}.{insertions}", position
            cigar += "I"
            continue

        index = position_indices[position]
        assert index > last_index, f"Position {position} out of order"
        cigar += "D" * (index - last_index - 1) + "M"
        last_index = index
        insertions = 0

    return cigar


def main():
    table = read_imgt_table()

    numbered_germlines: dict[str, tuple[ChainType, list[NumberedResidue]]] = {}
    for locus, genes in GENES_BY_LOCUS.items():
        germlines = read_germlines(GENE_DB_DIR / "gene_db" / "aa_genes" / "c_genes" / "human" / f"{locus}.fasta")
        deduplicated_germlines = read_germlines(
            GENE_DB_DIR / "gene_db" / "aa_genes_deduplicated" / "c_genes" / "human" / f"{locus}.fasta"
        )

        for gene in genes:
            germline = germlines[gene]
            assert germline == deduplicated_germlines[gene], f"{gene} differs between germline databases"

            if gene in HEAVY_GENES:
                chain_type = ChainType.HEAVY
                numbered = number_heavy_germline(gene, germline, table)
            else:
                chain_type = ChainType.LIGHT
                numbered = number_light_germline(gene, germline, table)

            validate(gene, germline, numbered, table)
            numbered_germlines[gene] = (chain_type, numbered)
            print(f"{gene}: {len(numbered)} numbered residues, unnumbered tail: {germline[len(numbered):][:10]}")

    for c_scheme in CScheme:
        output_path = GENE_DB_DIR / "scheme_mappings" / Organism.HOMO_SAPIENS.value / c_scheme.value
        output_path.mkdir(parents=True, exist_ok=True)

        with (output_path / "c_scheme_mapping.csv").open("w") as output_file:
            writer = csv.DictWriter(output_file, fieldnames=["gene_id", "scheme_cigar"])
            writer.writeheader()
            for gene, (chain_type, numbered) in numbered_germlines.items():
                # EU numbering is not defined for IGLC
                if all(residue.positions[c_scheme] is None for residue in numbered):
                    continue
                writer.writerow({"gene_id": gene, "scheme_cigar": scheme_cigar(numbered, c_scheme, chain_type)})


if __name__ == "__main__":
    main()
