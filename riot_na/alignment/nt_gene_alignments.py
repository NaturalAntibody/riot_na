import logging
from pathlib import Path
from typing import Optional

from riot_na.alignment.alignment_utils import offset_alignments
from riot_na.alignment.gene_aligner import (
    GeneAligner,
    create_aligner,
    create_v_gene_aligner,
)
from riot_na.alignment.segment_gene_aligner import SegmentGeneAligner
from riot_na.alignment.segment_gene_aligner import (
    create_v_gene_aligner as create_segment_v_gene_aligner,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentSegment,
    AlignmentsNT,
    GermlineGene,
    Locus,
    Organism,
)

logger = logging.getLogger(__name__)


class VDJCAlignerNT:
    def __init__(
        self,
        v_aligner: SegmentGeneAligner | GeneAligner,
        d_aligners: dict[Organism, GeneAligner],
        j_aligners: dict[Organism, dict[Locus, GeneAligner]],
        c_aligners: dict[Organism, dict[Locus, GeneAligner]],
    ) -> None:
        self.v_aligner = v_aligner
        self.d_aligners = d_aligners
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners  # at the time of writing, there is only one c germline database

    def produce_nt_alignments(self, query: str) -> list[AlignmentsNT]:
        results = []

        v_alignments = self.v_aligner.align(query, both_strains=True)
        for v_alignment in v_alignments:
            try:
                aligments_nt = self._produce_vdjc_nt_alignment(v_alignment)
                results.append(aligments_nt)
            except ValueError as e:
                logger.error(
                    "Error aligning %s: %s",
                    v_alignment.best_alignment.q_seq if v_alignment.best_alignment else "unknown sequence",
                    e,
                )
        return results

    def _produce_vdjc_nt_alignment(self, v_alignment_segment: AlignmentSegment) -> AlignmentsNT:
        aligments_nt = AlignmentsNT()
        aligments_nt.segment_start = v_alignment_segment.segment_start
        aligments_nt.segment_end = v_alignment_segment.segment_end

        v_alignment = v_alignment_segment.best_alignment
        if v_alignment is None:
            return aligments_nt

        query_v = v_alignment.q_seq
        aligments_nt.v = v_alignment

        # mask alignments
        query_j = v_alignment.q_seq[v_alignment.q_end :]

        # align j genes
        species = v_alignment.species
        locus = v_alignment.locus

        j_aligner = self.j_aligners[species][locus]
        j_alignment = j_aligner.align(query_j, both_strains=False)[0].best_alignment

        if j_alignment is None:
            return aligments_nt

        j_alignment = offset_alignments(v_alignment.q_end, j_alignment)
        assert j_alignment is not None
        aligments_nt.j = j_alignment

        # mask c input (query_db_clean -> valn_end + jaln_end)
        query_c = query_v[j_alignment.q_end :]

        if query_c:
            # align c
            c_aligner = self.c_aligners[species][locus]
            c_alignment = c_aligner.align(query_c, both_strains=False)[0].best_alignment

            if c_alignment is not None:
                c_alignment = offset_alignments(j_alignment.q_end, c_alignment)
                aligments_nt.c = c_alignment

        if locus == Locus.IGH:
            # mask d input (query_db_clean -> valn_end + jaln_end)
            query_d = query_v[v_alignment.q_end : j_alignment.q_start]
            if query_d:
                # align d genes
                d_aligner = self.d_aligners[species]
                d_alignment = d_aligner.align(query_d, both_strains=False)[0].best_alignment

                if d_alignment is not None:
                    d_alignment = offset_alignments(v_alignment.q_end, d_alignment)
                    aligments_nt.d = d_alignment

        return aligments_nt


def create_vdjc_aligner_nt(
    allowed_species: Optional[list[Organism]] = None, db_dir: Path = GENE_DB_DIR, use_segment_aligner: bool = False
):
    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS]

    v_aligner = (
        create_segment_v_gene_aligner(allowed_species=allowed_species, db_dir=db_dir)
        if use_segment_aligner
        else create_v_gene_aligner(allowed_species=allowed_species, db_dir=db_dir)
    )

    d_aligners = {}
    for organism in allowed_species:
        d_aligner = create_aligner(organism=organism, germline_gene=GermlineGene.D, locus=Locus.IGH, db_dir=db_dir)
        d_aligners[organism] = d_aligner

    j_aligners: dict[Organism, dict[Locus, GeneAligner]] = {}

    for organism in allowed_species:
        organism_aligners = j_aligners.get(organism, {})

        for locus in Locus:
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            j_aligner = create_aligner(organism=organism, germline_gene=GermlineGene.J, locus=locus, db_dir=db_dir)
            organism_aligners[locus] = j_aligner

        j_aligners[organism] = organism_aligners

    c_aligners: dict[Organism, dict[Locus, GeneAligner]] = {}

    for organism in allowed_species:

        organism_aligners = c_aligners.get(organism, {})

        if organism not in [Organism.HOMO_SAPIENS, Organism.VICUGNA_PACOS, Organism.CUSTOM]:
            c_organism = Organism.HOMO_SAPIENS
        else:
            c_organism = organism

        for locus in Locus:
            if organism == Organism.VICUGNA_PACOS and locus != Locus.IGH:
                continue

            c_aligner = create_aligner(organism=c_organism, germline_gene=GermlineGene.C, locus=locus, db_dir=db_dir)
            organism_aligners[locus] = c_aligner

        c_aligners[organism] = organism_aligners

    vdjc_aligner_nt = VDJCAlignerNT(v_aligner, d_aligners, j_aligners, c_aligners)

    return vdjc_aligner_nt
