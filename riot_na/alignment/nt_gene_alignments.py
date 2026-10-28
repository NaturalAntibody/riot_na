import logging
from pathlib import Path
from typing import Optional

from riot_na.alignment.alignment_utils import align_c, offset_alignments
from riot_na.alignment.gene_aligner import (
    GeneAligner,
    create_aligner,
    create_v_gene_aligner,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntryNT,
    AlignmentsNT,
    GermlineGene,
    Locus,
    Organism,
)

logger = logging.getLogger(__name__)


def align_djc_nt(
    v_alignment: AlignmentEntryNT,
    d_aligners: dict[Organism, GeneAligner],
    j_aligners: dict[Organism, dict[Locus, GeneAligner]],
    c_aligners: dict[Organism, dict[Locus, GeneAligner]],
) -> AlignmentsNT:
    """Align J, C and D genes to the part of the query following the V alignment. C and D are aligned only after J."""
    aligments_nt = AlignmentsNT(v=v_alignment)

    query_v = v_alignment.q_seq

    # mask alignments
    query_j = v_alignment.q_seq[v_alignment.q_end :]

    # align j genes
    species = v_alignment.species
    locus = v_alignment.locus

    j_aligner = j_aligners[species][locus]
    j_alignment = j_aligner.align(query_j, both_strains=False)

    if j_alignment is None:
        return aligments_nt

    j_alignment = offset_alignments(v_alignment.q_end, j_alignment)
    aligments_nt.j = j_alignment

    # mask c input (query_db_clean -> valn_end + jaln_end)
    query_c = query_v[j_alignment.q_end :]

    if query_c:
        # align c
        c_alignment = align_c(
            c_aligners[species],
            locus,
            query_c,
            lambda c_aligner: c_aligner.align(query_c, both_strains=False),
        )

        if c_alignment is not None:
            aligments_nt.c = offset_alignments(j_alignment.q_end, c_alignment)

    if locus == Locus.IGH:
        # mask d input (query_db_clean -> valn_end + jaln_end)
        query_d = query_v[v_alignment.q_end : j_alignment.q_start]
        if query_d:
            # align d genes
            d_aligner = d_aligners[species]
            d_alignment = d_aligner.align(query_d, both_strains=False)

            if d_alignment is not None:
                d_alignment = offset_alignments(v_alignment.q_end, d_alignment)
                aligments_nt.d = d_alignment

    return aligments_nt


class VDJCAlignerNT:
    def __init__(
        self,
        v_aligner: GeneAligner,
        d_aligners: dict[Organism, GeneAligner],
        j_aligners: dict[Organism, dict[Locus, GeneAligner]],
        c_aligners: dict[Organism, dict[Locus, GeneAligner]],
    ) -> None:
        self.v_aligner = v_aligner
        self.d_aligners = d_aligners
        self.j_aligners = j_aligners
        self.c_aligners = c_aligners  # at the time of writing, there is only one c germline database

    def produce_nt_alignments(self, query: str) -> AlignmentsNT:
        v_alignment = self.v_aligner.align(query, both_strains=True)
        if v_alignment is None:
            return AlignmentsNT()

        try:
            return align_djc_nt(v_alignment, self.d_aligners, self.j_aligners, self.c_aligners)
        except ValueError as e:
            logger.error("Error aligning %s: %s", v_alignment.q_seq, e)
            return AlignmentsNT()


DEFAULT_SPECIES = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)


def _create_organism_c_aligners_nt(organism: Organism, db_dir: Path) -> dict[Locus, GeneAligner]:
    """C aligners by locus of the organism, empty if there are no C germlines of the organism."""
    if not (Path(db_dir) / "gene_db" / "c_genes" / organism.value).exists():
        return {}

    return {
        locus: create_aligner(organism=organism, germline_gene=GermlineGene.C, locus=locus, db_dir=db_dir)
        for locus in Locus
        # at the moment we only have heavy genes
        if not (organism == Organism.VICUGNA_PACOS and locus != Locus.IGH)
    }


def _create_c_aligners_nt(
    allowed_species: tuple[Organism, ...], db_dir: Path, c_genes_organism: Optional[Organism]
) -> dict[Organism, dict[Locus, GeneAligner]]:
    """C aligners of the species of the V alignment: of its own C germlines or of c_genes_organism, if given."""
    if c_genes_organism is None:
        return {organism: _create_organism_c_aligners_nt(organism, db_dir) for organism in allowed_species}

    c_organism_aligners = _create_organism_c_aligners_nt(c_genes_organism, db_dir)
    if not c_organism_aligners:
        raise ValueError(f"There are no C germlines of {c_genes_organism} in {db_dir}")
    return {organism: c_organism_aligners for organism in allowed_species}


def create_djc_aligners_nt(
    allowed_species: tuple[Organism, ...], db_dir: Path, c_genes_organism: Optional[Organism] = None
) -> tuple[
    dict[Organism, GeneAligner], dict[Organism, dict[Locus, GeneAligner]], dict[Organism, dict[Locus, GeneAligner]]
]:
    """D, J and C aligners of the allowed species (by locus)."""
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

    c_aligners = _create_c_aligners_nt(allowed_species, db_dir, c_genes_organism)

    return d_aligners, j_aligners, c_aligners


def create_vdjc_aligner_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
    c_genes_organism: Optional[Organism] = None,
) -> VDJCAlignerNT:
    allowed_species = allowed_species or DEFAULT_SPECIES

    v_aligner = create_v_gene_aligner(allowed_species=allowed_species, db_dir=db_dir)
    d_aligners, j_aligners, c_aligners = create_djc_aligners_nt(allowed_species, db_dir, c_genes_organism)

    return VDJCAlignerNT(v_aligner, d_aligners, j_aligners, c_aligners)
