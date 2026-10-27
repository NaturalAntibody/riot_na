import logging
from pathlib import Path
from typing import Optional

from riot_na.alignment.gene_aligner import GeneAligner
from riot_na.alignment.nt_gene_alignments import (
    DEFAULT_SPECIES,
    align_djc_nt,
    create_djc_aligners_nt,
)
from riot_na.alignment.segment_gene_aligner import (
    SegmentGeneAligner,
    create_v_gene_aligner,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import Locus, Organism, SegmentAlignmentsNT

logger = logging.getLogger(__name__)


class SegmentVDJCAlignerNT:
    """Aligns V, D, J and C genes to every domain (segment) of a multiple domain query."""

    def __init__(
        self,
        v_aligner: SegmentGeneAligner,
        d_aligners: dict[Organism, GeneAligner],
        j_aligners: dict[Organism, dict[Locus, GeneAligner]],
        c_aligners: dict[Organism, dict[Locus, GeneAligner]],
    ) -> None:
        self.v_aligner = v_aligner
        self.d_aligners = d_aligners
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners

    def produce_nt_alignments(self, query: str) -> list[SegmentAlignmentsNT]:
        results = []

        for v_segment_alignment in self.v_aligner.align(query, both_strains=True):
            try:
                alignments = align_djc_nt(
                    v_segment_alignment.best_alignment, self.d_aligners, self.j_aligners, self.c_aligners
                )
            except ValueError as e:
                logger.error("Error aligning %s: %s", v_segment_alignment.best_alignment.q_seq, e)
                continue

            results.append(
                SegmentAlignmentsNT(
                    alignments=alignments,
                    segment_start=v_segment_alignment.segment_start,
                    segment_end=v_segment_alignment.segment_end,
                )
            )
        return results


def create_segment_vdjc_aligner_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> SegmentVDJCAlignerNT:
    allowed_species = allowed_species or DEFAULT_SPECIES

    v_aligner = create_v_gene_aligner(allowed_species=allowed_species, db_dir=db_dir)
    d_aligners, j_aligners, c_aligners = create_djc_aligners_nt(allowed_species, db_dir)

    return SegmentVDJCAlignerNT(v_aligner, d_aligners, j_aligners, c_aligners)
