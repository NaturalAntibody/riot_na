import csv
from pathlib import Path
from typing import Optional

from riot_na.config import GENE_DB_DIR
from riot_na.data.model import AlignmentString, GeneId, Organism, Scheme


class SchemeMappingFacade:
    mappings: dict[GeneId, AlignmentString]

    def __init__(self, scheme: Scheme, allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR):
        if not allowed_species:
            allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS]

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


if __name__ == "__main__":
    scheme_mapping_data_path = GENE_DB_DIR
    scheme_mapping_db = SchemeMappingFacade(
        scheme=Scheme.IMGT, allowed_species=[Organism.HOMO_SAPIENS], db_dir=scheme_mapping_data_path
    )

    print(scheme_mapping_db.mappings["IGHV1-2*07"])
    print(scheme_mapping_db.mappings["IGKJ4*01"])
