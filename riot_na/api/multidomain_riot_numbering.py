from pathlib import Path
from typing import Optional

from cachetools import cached

from riot_na.airr.airr_builder import AirrBuilder
from riot_na.airr.airr_builder_aa import AirrBuilderAA
from riot_na.alignment.aa_gene_alignments import (
    VJCAlignmentTranslatorAA,
    create_vjc_alignment_translator_aa,
)
from riot_na.alignment.segment_aa_gene_alignments import (
    SegmentVJCAlignerAA,
    create_segment_vjc_aligner_aa,
)
from riot_na.alignment.segment_nt_gene_alignments import (
    SegmentVDJCAlignerNT,
    create_segment_vdjc_aligner_nt,
)
from riot_na.api.domain_numbering import number_domain_aa, number_domain_nt
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    CScheme,
    Organism,
    Scheme,
    SegmentedAirrRearrangementEntryAA,
    SegmentedAirrRearrangementEntryNT,
)
from riot_na.schemes.scheme_alignment import SchemeAligner


class MultidomainRiotNumberingNT:
    def __init__(
        self,
        segment_vdjc_aligner_nt: SegmentVDJCAlignerNT,
        vjc_alignment_translator_aa: VJCAlignmentTranslatorAA,
        scheme_aligner: SchemeAligner,
    ):
        self.segment_vdjc_aligner = segment_vdjc_aligner_nt
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
    ) -> list[SegmentedAirrRearrangementEntryNT]:
        """
        Align and number all domains of a multiple domain nucleotide sequence.

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
        SegmentedAirrRearrangementEntryNT object of every domain with alignment results (Extended AIRR format) and
        the domain position in the query sequence. A single entry without alignments if no domain is found.
        """
        query_sequence = query_sequence.upper().replace("U", "T")

        domains = self.segment_vdjc_aligner.produce_nt_alignments(query_sequence)
        if not domains:
            rearrangement = SegmentedAirrRearrangementEntryNT(
                sequence_header=header, sequence="", numbering_scheme=scheme.value, query_sequence=query_sequence
            )
            return [AirrBuilder(rearrangement).build()]

        results = []
        for domain in domains:
            # Domain sequence might have been reverse complemented during V alignment.
            sequence = domain.alignments.v.q_seq if domain.alignments.v else ""
            rearrangement = SegmentedAirrRearrangementEntryNT(
                sequence_header=header,
                sequence=sequence,
                numbering_scheme=scheme.value,
                query_sequence=query_sequence,
                segment_start=domain.segment_start + 1,
                segment_end=domain.segment_end,
            )
            results.append(
                number_domain_nt(
                    AirrBuilder(rearrangement),
                    domain.alignments,
                    self.vjc_alignment_translator_aa,
                    self.scheme_aligner,
                    scheme,
                    extend_alignment,
                    c_scheme,
                )
            )
        return results


class MultidomainRiotNumberingAA:
    def __init__(self, segment_vjc_aligner_aa: SegmentVJCAlignerAA, scheme_aligner: SchemeAligner):
        self.segment_vjc_aligner_aa = segment_vjc_aligner_aa
        self.scheme_aligner = scheme_aligner

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        *,
        c_scheme: Optional[CScheme] = None,
    ) -> list[SegmentedAirrRearrangementEntryAA]:
        """
        Align and number all domains of a multiple domain amino acid sequence.

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
        - c_scheme : Numbering scheme of the constant region (eu, kabat or imgt). Only human IGHG1-4, IGKC and IGLC
        constant regions are numbered, EU numbering is not defined for IGLC. By default constant region is not numbered.

        Returns:
        -----------
        SegmentedAirrRearrangementEntryAA object of every domain with alignment results (Extended AIRR format) and
        the domain position in the query sequence. A single entry without alignments if no domain is found.
        """
        query_sequence = query_sequence.upper()

        domains = self.segment_vjc_aligner_aa.produce_aa_alignments(query_sequence)
        if not domains:
            rearrangement = SegmentedAirrRearrangementEntryAA(
                sequence_header=header,
                sequence_aa=query_sequence,
                numbering_scheme=scheme.value,
                query_sequence=query_sequence,
            )
            return [AirrBuilderAA(rearrangement).build()]

        results = []
        for domain in domains:
            sequence_aa = query_sequence[domain.segment_start : domain.segment_end]
            rearrangement = SegmentedAirrRearrangementEntryAA(
                sequence_header=header,
                sequence_aa=sequence_aa,
                numbering_scheme=scheme.value,
                query_sequence=query_sequence,
                segment_start=domain.segment_start + 1,
                segment_end=domain.segment_end,
            )
            results.append(
                number_domain_aa(
                    AirrBuilderAA(rearrangement),
                    domain.alignments,
                    self.scheme_aligner,
                    scheme,
                    extend_alignment,
                    c_scheme,
                )
            )
        return results


def create_multidomain_riot_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> MultidomainRiotNumberingNT:
    """
    Factory function for creating MultidomainRiotNumberingNT object, which numbers all domains of multiple domain
    sequences.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism tuple (default use all species - mouse, human and alpaca).
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    MultidomainRiotNumberingNT object used for performing alignment and numbering.
    """
    segment_vdjc_aligner_nt = create_segment_vdjc_aligner_nt(allowed_species=allowed_species, db_dir=db_dir)
    vjc_alignment_translator_aa = create_vjc_alignment_translator_aa(allowed_species=allowed_species, db_dir=db_dir)
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return MultidomainRiotNumberingNT(segment_vdjc_aligner_nt, vjc_alignment_translator_aa, scheme_aligner)


def create_multidomain_riot_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> MultidomainRiotNumberingAA:
    """
    Factory function for creating MultidomainRiotNumberingAA object, which numbers all domains of multiple domain
    sequences.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism tuple (default use all species - mouse, human and alpaca).
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    MultidomainRiotNumberingAA object used for performing alignment and numbering.
    """
    segment_vjc_aligner_aa = create_segment_vjc_aligner_aa(allowed_species=allowed_species, db_dir=db_dir)
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return MultidomainRiotNumberingAA(segment_vjc_aligner_aa, scheme_aligner)


@cached({})
def get_or_create_multidomain_riot_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> MultidomainRiotNumberingNT:
    return create_multidomain_riot_nt(allowed_species=allowed_species, db_dir=db_dir)


@cached({})
def get_or_create_multidomain_riot_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    db_dir: Path = GENE_DB_DIR,
) -> MultidomainRiotNumberingAA:
    return create_multidomain_riot_aa(allowed_species=allowed_species, db_dir=db_dir)
