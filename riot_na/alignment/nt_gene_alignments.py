import logging
from pathlib import Path
from typing import Optional

from riot_na.alignment.alignment_utils import offset_alignments
from riot_na.alignment.gene_aligner import GeneAligner, create_aligner
from riot_na.alignment.segment_gene_aligner import (
    SegmentAlignmentResult,
    SegmentGeneAligner,
    create_v_gene_aligner,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import AlignmentsNT, GermlineGene, Locus, Organism

logger = logging.getLogger(__name__)


class VDJCAlignerNT:
    def __init__(
        self,
        v_aligner: SegmentGeneAligner,
        d_aligners: dict[Organism, GeneAligner],
        j_aligners: dict[Organism, dict[Locus, GeneAligner]],
        c_aligners: dict[Locus, GeneAligner],
    ) -> None:
        self.v_aligner = v_aligner
        self.d_aligners = d_aligners
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners  # at the time of writing, there is only one c germline database

    def produce_nt_alignments(self, query: str) -> list[AlignmentsNT]:
        results = []

        v_segment_alignments: SegmentAlignmentResult = self.v_aligner.align_segments(query, both_strains=True)
        v_segments_alignments = v_segment_alignments.segment_alignments
        for v_segment_alignment in v_segments_alignments:
            result = AlignmentsNT()
            try:
                result.segment_start = v_segment_alignment.segment_start
                result.segment_end = v_segment_alignment.segment_end

                v_alignment = v_segment_alignment.best_alignment
                if v_alignment is None:
                    continue

                query_v = v_alignment.q_seq
                result.v = v_alignment

                # mask alignments
                query_j = v_alignment.q_seq[v_alignment.q_end :]

                # align j genes
                species = v_alignment.species
                locus = v_alignment.locus

                j_aligner = self.j_aligners[species][locus]
                j_alignment = j_aligner.align(query_j, both_strains=False)

                if j_alignment is None:
                    continue

                j_alignment = offset_alignments(v_alignment.q_end, j_alignment)
                assert j_alignment is not None
                result.j = j_alignment

                # mask c input (query_db_clean -> valn_end + jaln_end)
                query_c = query_v[j_alignment.q_end :]

                if query_c:
                    # align c
                    c_aligner = self.c_aligners[locus]
                    c_alignment = c_aligner.align(query_c, both_strains=False)

                    if c_alignment is not None:
                        c_alignment = offset_alignments(j_alignment.q_end, c_alignment)
                        result.c = c_alignment

                # mask d input (query_db_clean -> valn_end + jaln_end)
                query_d = query_v[v_alignment.q_end : j_alignment.q_start]
                if query_d:
                    # align d genes
                    d_aligner = self.d_aligners[species]
                    d_alignment = d_aligner.align(query_d, both_strains=False)

                    if d_alignment is not None:
                        d_alignment = offset_alignments(v_alignment.q_end, d_alignment)
                        result.d = d_alignment

            except ValueError as e:
                logger.error(
                    "Error aligning segment %s:%s: %s",
                    v_segment_alignment.segment_start,
                    v_segment_alignment.segment_end,
                    e,
                )
                continue
            finally:
                if result.v or result.j or result.c:
                    results.append(result)

        return results


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
