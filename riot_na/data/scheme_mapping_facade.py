import csv
from pathlib import Path
from typing import Optional

from riot_na.config import GENE_DB_DIR
from riot_na.data.model import AlignmentString, CScheme, GeneId, Organism, Scheme


class SchemeMappingFacade:
    mappings: dict[GeneId, AlignmentString]

    def __init__(
        self, scheme: Scheme, allowed_species: Optional[tuple[Organism, ...]] = None, db_dir: Path = GENE_DB_DIR
    ):
        if not allowed_species:
            allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

        self.mappings = {}

        for species in allowed_species:
            file_path = db_dir / "scheme_mappings" / species.value / scheme.value / "scheme_mapping.csv"

            with file_path.open() as scheme_mapping_file:
                scheme_mapping_reader = csv.DictReader(scheme_mapping_file)
                for row in scheme_mapping_reader:
                    self.mappings[species.value + "|" + row["gene_id"]] = AlignmentString(row["scheme_cigar"])

    def get_mapping(self, organism: Organism, gene_id: str) -> AlignmentString:
        gene_id = gene_id.split("|")[-1]
        return AlignmentString(self.mappings[organism.value + "|" + gene_id])


class CSchemeMappingFacade:
    """
    Scheme mappings of C germlines. Only some C germlines have them (human IGHG1-4, IGKC and IGLC, no EU for IGLC),
    for other germlines get_mapping returns None.
    """

    mappings: dict[GeneId, AlignmentString]

    def __init__(self, c_scheme: CScheme, db_dir: Path = GENE_DB_DIR):
        self.mappings = {}

        for species in Organism:
            file_path = db_dir / "scheme_mappings" / species.value / c_scheme.value / "c_scheme_mapping.csv"
            if not file_path.exists():
                continue

            with file_path.open() as scheme_mapping_file:
                scheme_mapping_reader = csv.DictReader(scheme_mapping_file)
                for row in scheme_mapping_reader:
                    self.mappings[species.value + "|" + row["gene_id"]] = AlignmentString(row["scheme_cigar"])

    def get_mapping(self, organism: Organism, gene_id: str) -> Optional[AlignmentString]:
        # species without own C germlines are aligned to human ones (see create_vjc_aligner_aa)
        if organism not in (Organism.HOMO_SAPIENS, Organism.VICUGNA_PACOS, Organism.CUSTOM):
            organism = Organism.HOMO_SAPIENS

        gene_id = gene_id.split("|")[-1]
        return self.mappings.get(organism.value + "|" + gene_id)


if __name__ == "__main__":
    scheme_mapping_data_path = GENE_DB_DIR
    scheme_mapping_db = SchemeMappingFacade(
        scheme=Scheme.IMGT, allowed_species=(Organism.HOMO_SAPIENS,), db_dir=scheme_mapping_data_path
    )

    print(scheme_mapping_db.mappings["IGHV1-2*07"])
    print(scheme_mapping_db.mappings["IGKJ4*01"])
