from pathlib import Path
from typing import Optional

from cachetools import cached

from riot_na.airr.airr_builder import AirrBuilder
from riot_na.airr.airr_builder_aa import AirrBuilderAA
from riot_na.alignment.aa_gene_alignments import (
    VJCAlignerAA,
    VJCAlignmentTranslatorAA,
    create_vjc_aligner_aa,
    create_vjc_alignment_translator_aa,
)
from riot_na.alignment.nt_gene_alignments import VDJCAlignerNT, create_vdjc_aligner_nt
from riot_na.api.domain_numbering import number_domain_aa, number_domain_nt
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    CScheme,
    Organism,
    Scheme,
)
from riot_na.schemes.scheme_alignment import SchemeAligner


class RiotNumberingNT:
    def __init__(
        self,
        vdjc_aligner_nt: VDJCAlignerNT,
        vjc_alignment_translator_aa: VJCAlignmentTranslatorAA,
        scheme_aligner: SchemeAligner,
    ):
        self.vdjc_aligner = vdjc_aligner_nt
        self.vjc_alignment_translator_aa = vjc_alignment_translator_aa
        self.scheme_aligner = scheme_aligner

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        *,
        c_scheme: Optional[CScheme] = None,
    ) -> AirrRearrangementEntryNT:
        """
        Align and number nucleotide sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Nucleotide sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).
        - extend_alignment: Extend alignments to cover full sequence.
        - c_scheme : Numbering scheme of the constant region (eu, kabat or imgt). Only human IGHG1-4, IGKC and IGLC
        constant regions are numbered, EU numbering is not defined for IGLC. By default constant region is not numbered.

        Returns:
        -----------
        AirrRearrangementEntryNT object with alignment results (Extended AIRR format)
        """
        query_sequence = query_sequence.upper().replace("U", "T")

        # produce alignments for v,d,j,c
        alignments = self.vdjc_aligner.produce_nt_alignments(query_sequence)

        # Query might have been reverse complemented during V alignment.
        sequence = alignments.v.q_seq if alignments.v else query_sequence
        rearrangement = AirrRearrangementEntryNT(
            sequence_header=header, sequence=sequence, numbering_scheme=scheme.value
        )
        return number_domain_nt(
            AirrBuilder(rearrangement),
            alignments,
            self.vjc_alignment_translator_aa,
            self.scheme_aligner,
            scheme,
            extend_alignment,
            c_scheme,
        )


class RiotNumberingAA:
    def __init__(self, vjc_aligner_aa: VJCAlignerAA, scheme_aligner: SchemeAligner):
        self.vjc_aligner_aa = vjc_aligner_aa
        self.scheme_aligner = scheme_aligner

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        *,
        v_gene: Optional[str] = None,
        j_gene: Optional[str] = None,
        c_scheme: Optional[CScheme] = None,
    ) -> AirrRearrangementEntryAA:
        """
        Align and number amino acid sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Amino acid sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).
        - extend_alignment: RIOT uses Striped Smith-Waterman local alignemtn algorithm.
        Due to that, if query sequence has too many differences (due to mutations or other reasons)
        to the germline sequence at 5' or / and 3' ends the result may not show that part of the gene.
        This option appends missing amino acids, so the result covers full gene sequence. Beware that
        those query fragments are not truly a part of alignemnt and are shown only for convenience!
        - v_gene, j_gene : Align to these exact germlines (e.g. IGHV3-7*01, IGHJ4*02) instead of
        searching the germline database. Both have to be given together. Gene databases are not
        prefiltered then, but alignment thresholds still apply. Germlines are looked up in the species
        the numbering was created for, human germlines first.
        - c_scheme : Numbering scheme of the constant region (eu, kabat or imgt). Only human IGHG1-4, IGKC and IGLC
        constant regions are numbered, EU numbering is not defined for IGLC. By default constant region is not numbered.

        Returns:
        -----------
        AirrRearrangementEntryAA object with alignment results (Extended AIRR format)
        """
        assert (v_gene is None) == (j_gene is None), "Both v_gene and j_gene have to be specified"
        query_sequence = query_sequence.upper()

        # produce alignments for v,j,c
        if v_gene and j_gene:
            alignments = self.vjc_aligner_aa.produce_aa_alignments_for_genes(
                query_sequence, v_gene=v_gene, j_gene=j_gene
            )
        else:
            alignments = self.vjc_aligner_aa.produce_aa_alignments(query_sequence)

        rearrangement = AirrRearrangementEntryAA(
            sequence_header=header, sequence_aa=query_sequence, numbering_scheme=scheme.value
        )
        return number_domain_aa(
            AirrBuilderAA(rearrangement), alignments, self.scheme_aligner, scheme, extend_alignment, c_scheme
        )


def create_riot_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
    c_genes_organism: Optional[Organism] = None,
) -> RiotNumberingNT:
    """
    Factory function for creating RiotNumberingNT object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism tuple (default use all species - mouse, human and alpaca).
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).
    - c_genes_organism : Align constant regions to C germlines of this organism only, for V alignments of any
    species. Default is to use C germlines of the species of the V alignment (mouse has none, so its constant regions
    are not aligned). Useful e.g. for humanization analysis.

    Returns:
    -----------
    RiotNumberingNT object used for performing alignment and numbering.
    """
    vdjc_aligner_nt = create_vdjc_aligner_nt(
        allowed_species=allowed_species, db_dir=db_dir, c_genes_organism=c_genes_organism
    )
    vjc_alignment_translator_aa = create_vjc_alignment_translator_aa(
        allowed_species=allowed_species, db_dir=db_dir, c_genes_organism=c_genes_organism
    )
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingNT(vdjc_aligner_nt, vjc_alignment_translator_aa, scheme_aligner)


def create_riot_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
    c_genes_organism: Optional[Organism] = None,
) -> RiotNumberingAA:
    """
    Factory function for creating RiotNumberingAA object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism tuple (default use all species - mouse, human and alpaca).
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).
    - c_genes_organism : Align constant regions to C germlines of this organism only, for V alignments of any
    species. Default is to use C germlines of the species of the V alignment (mouse has none, so its constant regions
    are not aligned). Useful e.g. for humanization analysis.

    Returns:
    -----------
    RiotNumberingAA object used for performing alignment and numbering.
    """
    vjc_aligner_aa = create_vjc_aligner_aa(
        allowed_species=allowed_species, db_dir=db_dir, c_genes_organism=c_genes_organism
    )
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingAA(vjc_aligner_aa, scheme_aligner)


@cached({})
def get_or_create_riot_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
    c_genes_organism: Optional[Organism] = None,
) -> RiotNumberingNT:
    return create_riot_nt(allowed_species=allowed_species, db_dir=db_dir, c_genes_organism=c_genes_organism)


@cached({})
def get_or_create_riot_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
    c_genes_organism: Optional[Organism] = None,
) -> RiotNumberingAA:
    return create_riot_aa(allowed_species=allowed_species, db_dir=db_dir, c_genes_organism=c_genes_organism)
