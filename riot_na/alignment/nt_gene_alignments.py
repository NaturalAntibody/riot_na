from pathlib import Path
from typing import Optional

from riot_na.alignment.alignment_utils import offset_alignments
from riot_na.alignment.gene_aligner import (
    GeneAligner,
    create_aligner,
    create_v_gene_aligner,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import AlignmentsNT, GermlineGene, Locus, Organism


class VDJCAlignerNT:
    def __init__(
        self,
        v_aligner: GeneAligner,
        d_aligners: dict[Organism, GeneAligner],
        j_aligners: dict[Organism, dict[Locus, GeneAligner]],
        c_aligners: dict[Locus, GeneAligner],
    ) -> None:
        self.v_aligner = v_aligner
        self.d_aligners = d_aligners
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners  # at the time of writing, there is only one c germline database

    def produce_nt_alignments(self, query: str) -> AlignmentsNT:
        result = AlignmentsNT()

        best_v_alignment = self.v_aligner.align(query, both_strains=True)

        if best_v_alignment is None:
            return result

        query_v = best_v_alignment.q_seq
        result.v = best_v_alignment

        # mask alignments
        query_j = best_v_alignment.q_seq[best_v_alignment.q_end :]

        # align j genes
        species = best_v_alignment.species
        locus = best_v_alignment.locus

        j_aligner = self.j_aligners[species][locus]
        best_j_alignment = j_aligner.align(query_j, both_strains=False)

        if best_j_alignment is None:
            return result

        best_j_alignment = offset_alignments(best_v_alignment.q_end, best_j_alignment)
        assert best_j_alignment is not None
        result.j = best_j_alignment

        # mask c input (query_db_clean -> valn_end + jaln_end)
        query_c = query_v[best_j_alignment.q_end :]

        # align c
        c_aligner = self.c_aligners[locus]
        best_c_alignment = c_aligner.align(query_c, both_strains=False)

        if best_c_alignment is not None:
            best_c_alignment = offset_alignments(best_j_alignment.q_end, best_c_alignment)
            result.c = best_c_alignment

        # mask d input (query_db_clean -> valn_end + jaln_end)
        query_d = query_v[best_v_alignment.q_end : best_j_alignment.q_start]

        # align c genes
        d_aligner = self.d_aligners[species]
        best_d_alignment = d_aligner.align(query_d, both_strains=False)

        if best_d_alignment is not None:
            best_d_alignment = offset_alignments(best_v_alignment.q_end, best_d_alignment)
            result.d = best_d_alignment

        return result


def create_vdjc_aligner_nt(allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR):
    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS]

    v_aligner = create_v_gene_aligner(allowed_species=allowed_species, db_dir=db_dir)

    d_aligners = {}
    for organism in allowed_species:
        d_aligner = create_aligner(
            allowed_species=organism, germline_gene=GermlineGene.D, locus=Locus.IGH, db_dir=db_dir
        )
        d_aligners[organism] = d_aligner

    j_aligners: dict[Organism, dict[Locus, GeneAligner]] = {}

    for organism in allowed_species:
        for locus in Locus:
            j_aligner = create_aligner(
                allowed_species=organism, germline_gene=GermlineGene.J, locus=locus, db_dir=db_dir
            )
            organism_aligners = j_aligners.get(organism, {})
            organism_aligners[locus] = j_aligner
            j_aligners[organism] = organism_aligners

    c_aligners = {}

    for locus in Locus:
        c_organism = Organism.HOMO_SAPIENS if organism != Organism.CUSTOM else Organism.CUSTOM

        c_aligner = create_aligner(allowed_species=c_organism, germline_gene=GermlineGene.C, locus=locus, db_dir=db_dir)
        c_aligners[locus] = c_aligner

    vdjc_aligner_nt = VDJCAlignerNT(v_aligner, d_aligners, j_aligners, c_aligners)

    return vdjc_aligner_nt
