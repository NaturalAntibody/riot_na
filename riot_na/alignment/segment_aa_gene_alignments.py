import logging
from pathlib import Path
from typing import Optional

from riot_na.alignment.aa_gene_alignments import (
    DEFAULT_SPECIES,
    aa_genes_dir,
    align_jc_aa,
    create_jc_aligners_aa,
)
from riot_na.alignment.gene_aligner import GeneAlignerAA
from riot_na.alignment.segment_gene_aligner import (
    SegmentGeneAlignerAA,
    create_aa_v_gene_aligner,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import Locus, Organism, SegmentAlignmentsAA

logger = logging.getLogger(__name__)


class SegmentVJCAlignerAA:
    """Aligns V, J and C genes to every domain (segment) of a multiple domain query."""

    def __init__(
        self,
        v_aligner: SegmentGeneAlignerAA,
        j_aligners: dict[Organism, dict[Locus, GeneAlignerAA]],
        c_aligners: Optional[dict[Organism, dict[Locus, GeneAlignerAA]]] = None,
    ) -> None:
        self.v_aligner = v_aligner
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners

    def produce_aa_alignments(self, query: str) -> list[SegmentAlignmentsAA]:
        results = []

        for v_segment_alignment in self.v_aligner.align(query):
            try:
                alignments = align_jc_aa(v_segment_alignment.best_alignment, self.j_aligners, self.c_aligners)
            except ValueError as e:
                logger.error("Error aligning segment %s due to %s", v_segment_alignment.best_alignment.q_seq, e)
                continue

            results.append(
                SegmentAlignmentsAA(
                    alignments=alignments,
                    segment_start=v_segment_alignment.segment_start,
                    segment_end=v_segment_alignment.segment_end,
                )
            )
        return results


def create_segment_vjc_aligner_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
    c_genes_organism: Optional[Organism] = None,
) -> SegmentVJCAlignerAA:
    allowed_species = allowed_species or DEFAULT_SPECIES
    genes_dir = aa_genes_dir(db_dir)

    v_aligner = create_aa_v_gene_aligner(allowed_species=allowed_species, aa_genes_dir=genes_dir)
    j_aligners, c_aligners = create_jc_aligners_aa(allowed_species, genes_dir, c_genes_organism)

    return SegmentVJCAlignerAA(v_aligner, j_aligners, c_aligners)
