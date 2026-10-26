# pylint: disable=too-many-branches, too-many-statements

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, Literal, Mapping, Optional, Sequence, cast, overload

from cachetools import cached

from riot_na.airr.airr_builder import AirrBuilder, SegmentedAirrBuilder
from riot_na.airr.airr_builder_aa import AirrBuilderAA, SegmentedAirrBuilderAA
from riot_na.alignment.aa_gene_alignments import (
    VJCAlignerAA,
    VJCAlignmentTranslatorAA,
    create_vjc_aligner_aa,
    create_vjc_alignment_translator_aa,
)
from riot_na.alignment.alignment_utils import v_domain_end
from riot_na.alignment.constant_region_aligner import ConstantRegionAligner
from riot_na.alignment.gene_aligner import GeneAlignerAA
from riot_na.alignment.nt_gene_alignments import VDJCAlignerNT, create_vdjc_aligner_nt
from riot_na.alignment.segment_gene_aligner import (
    SegmentGeneAligner,
    SegmentGeneAlignerAA,
)
from riot_na.api.chain_annotation import (
    TranslatedDomainAlignment,
    VDomainSpan,
    align_constant_domains_aa,
    align_constant_domains_nt,
    annotate_constant_domains,
    chain_domains,
    chain_layout,
    reverse_complement,
)
from riot_na.api.molecule_annotation import PairingTarget, pair_binding_sites
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    AlignmentsAA,
    AlignmentsNT,
    ChainAnnotation,
    ChainType,
    CScheme,
    Locus,
    MoleculeAnnotation,
    Organism,
    Scheme,
    SegmentedAirrRearrangementEntry,
    SegmentedAirrRearrangementEntryAA,
    SegmentedAirrRearrangementEntryNT,
)
from riot_na.schemes.region_offsets import (
    infer_aa_region_offsets,
    infer_nt_region_offsets,
    infer_region_offsets,
)
from riot_na.schemes.scheme_alignment import (
    SchemeAligner,
    fix_imgt_cdrs_numbering,
    get_positional_scheme_mapping,
    get_scheme_residue_mapping,
    produce_numbering,
    sort_imgt_numbering,
)

MULTIPLE_DOMAINS_NOT_ENABLED = (
    "return_all_domains requires numbering created with return_all_domains=True, "
    "e.g. create_riot_aa(return_all_domains=True)"
)
CHAIN_ANNOTATION_NOT_ENABLED = (
    "annotate_chain requires numbering created with return_all_domains=True, "
    "e.g. create_riot_aa(return_all_domains=True)"
)


def create_constant_region_aligner(
    c_aligners: Optional[dict[Organism, dict[Locus, GeneAlignerAA]]], scheme_aligner: SchemeAligner
) -> Optional[ConstantRegionAligner]:
    if not c_aligners:
        return None
    return ConstantRegionAligner(c_aligners, scheme_aligner.c_scheme_mapping_facades[CScheme.IMGT])


def v_domain_span_aa(alignments: AlignmentsAA) -> VDomainSpan:
    assert alignments.v
    segment_start = alignments.segment_start or 0
    return VDomainSpan(
        segment_start + alignments.v.q_start,
        # J is aligned to the query left after V alignment
        segment_start + v_domain_end(alignments.v, alignments.j, j_query_start=alignments.v.q_end),
        alignments.v.species,
        alignments.v.locus,
    )


def v_domain_span_nt(alignments: AlignmentsNT) -> VDomainSpan:
    assert alignments.v
    segment_start = alignments.segment_start or 0
    return VDomainSpan(
        segment_start + alignments.v.q_start,
        segment_start + v_domain_end(alignments.v, alignments.j),
        alignments.v.species,
        alignments.v.locus,
    )


def number_c_region(
    scheme_aligner: SchemeAligner, query: str, aa_alignments: AlignmentsAA, c_scheme: CScheme
) -> Optional[tuple[dict[str, str], dict[int, str]]]:
    """
    Number C region of the query. Returns scheme residue mapping and positional scheme mapping
    or None if there is no C alignment or the C gene has no scheme mapping.
    """
    if not (aa_alignments.v and aa_alignments.c):
        return None

    # C is aligned to the query left after V and J alignments (or V alignment, if J was not found)
    c_query_offset = v_domain_end(aa_alignments.v, aa_alignments.j, j_query_start=aa_alignments.v.q_end)
    return scheme_aligner.number_c_alignment(query, aa_alignments.c, c_scheme, c_query_offset)


class _RiotNumbering(ABC, Generic[SegmentedAirrRearrangementEntry]):
    """Annotation of chains and molecules, shared by nucleotide and amino acid numbering."""

    # length of residues in chain sequence positions
    residue_length: int
    scheme_aligner: SchemeAligner

    @abstractmethod
    def annotate_chain(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        c_scheme: Optional[CScheme] = None,
    ) -> ChainAnnotation[SegmentedAirrRearrangementEntry]: ...

    def annotate_molecule(
        self,
        molecule_id: str,
        chains: Mapping[str, str],
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        c_scheme: Optional[CScheme] = None,
        pairs: Optional[Sequence[tuple[PairingTarget, PairingTarget]]] = None,
    ) -> MoleculeAnnotation[SegmentedAirrRearrangementEntry]:
        """
        Annotate all chains of a molecule (annotate_chain) and pair their V domains into binding sites:
        scFv, diabodies, Fab (including domain crossovers and V domains stacked on them), single domains.
        Requires numbering created with return_all_domains=True.

        Parameters:
        -----------
        - molecule_id : Molecule ID, chain headers are molecule_id|chain_id.
        - chains : Sequences of the molecule chains (nucleotide or amino acid sequences, as of the numbering) by chain
        ID, every distinct chain once.
        - scheme : Numbering scheme of V domains (default IMGT).
        - extend_alignment: Extend alignments of V domains to cover full sequence.
        - c_scheme : Numbering scheme of constant domains (eu, kabat or imgt). By default constant domains are not
        numbered.
        - pairs : Pairing hints, pairs of chains or V domains (chain ID and domain index, 1-based position among all
        domains of the chain) which pair, e.g. [("HC1", "LC1"), ("HC2", "LC2")] or [(("HC", 1), "LC1"), (("HC", 3),
        "LC2")]. V domains given in pairs pair only with V domains they are given with, except light chain Fab arms
        (followed by CL) and light V domains paired by pairs only, which can be shared by heavy chains. By default
        V domains with more than one possible partner are ambiguous and V domains not paired by their format (e.g.
        V domain only chains) are single domains or unpaired.

        Returns:
        -----------
        MoleculeAnnotation with chain annotations and binding sites.
        """
        chain_annotations = {
            chain_id: self.annotate_chain(f"{molecule_id}|{chain_id}", sequence, scheme, extend_alignment, c_scheme)
            for chain_id, sequence in chains.items()
        }
        return MoleculeAnnotation(
            molecule_id=molecule_id,
            chains=chain_annotations,
            binding_sites=pair_binding_sites(chain_annotations, residue_length=self.residue_length, pairs=pairs),
        )

    def _chain_annotation(  # pylint: disable=too-many-arguments
        self,
        header: str,
        sequence: str,
        rev_comp: bool,
        rearrangements: list[SegmentedAirrRearrangementEntry],
        v_domains: Sequence[VDomainSpan],
        domain_alignments: Sequence[TranslatedDomainAlignment],
        c_scheme: Optional[CScheme],
    ) -> ChainAnnotation[SegmentedAirrRearrangementEntry]:
        """Chain annotation with V domains and constant domains of the domain alignments, numbered in c_scheme."""
        constant_domains = annotate_constant_domains(self.scheme_aligner, domain_alignments, v_domains, c_scheme)
        domains = chain_domains(v_domains, constant_domains)
        return ChainAnnotation(
            sequence_header=header,
            sequence=sequence,
            rev_comp=rev_comp,
            layout=chain_layout(len(sequence), domains),
            rearrangements=rearrangements,
            constant_domains=constant_domains,
            domains=domains,
        )


class RiotNumberingNT(_RiotNumbering[SegmentedAirrRearrangementEntryNT]):
    residue_length = 3

    def __init__(
        self,
        vdjc_aligner_nt: VDJCAlignerNT,
        vjc_alignment_translator_aa: VJCAlignmentTranslatorAA,
        scheme_aligner: SchemeAligner,
    ):
        self.vdjc_aligner = vdjc_aligner_nt
        self.vjc_alignment_translator_aa = vjc_alignment_translator_aa
        self.scheme_aligner = scheme_aligner
        self.constant_region_aligner = create_constant_region_aligner(
            vjc_alignment_translator_aa.c_aligners, scheme_aligner
        )

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        *,
        c_scheme: Optional[CScheme] = ...,
    ) -> AirrRearrangementEntryNT: ...

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        return_all_domains: Literal[True] = ...,
        *,
        c_scheme: Optional[CScheme] = ...,
    ) -> list[AirrRearrangementEntryNT]: ...

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        return_all_domains: Literal[False] = ...,
        *,
        c_scheme: Optional[CScheme] = ...,
    ) -> AirrRearrangementEntryNT: ...

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        return_all_domains: bool = ...,
        *,
        c_scheme: Optional[CScheme] = ...,
    ) -> AirrRearrangementEntryNT | list[AirrRearrangementEntryNT]: ...

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        return_all_domains: bool = False,
        c_scheme: Optional[CScheme] = None,
    ) -> AirrRearrangementEntryNT | list[AirrRearrangementEntryNT]:  # pylint: disable=too-many-statements
        """
        Align and number nucleotide sequence.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Nucleotide sequence to number.
        - scheme : Which numbering scheme to use (default IMGT).
        - return_all_domains: If True, return all domains of multiple domain proteins. Requires numbering created
        with return_all_domains=True, which returns the best scoring domain only if return_all_domains is False.
        - extend_alignment: Extend alignments to cover full sequence.
        - c_scheme : Numbering scheme of the constant region (eu, kabat or imgt). Only human IGHG1-4, IGKC and IGLC
        constant regions are numbered, EU numbering is not defined for IGLC. By default constant region is not numbered.

        Returns:
        -----------
        AirrRearrangementEntryNT object with alignment results (Extended AIRR format)
        """
        if return_all_domains and not isinstance(self.vdjc_aligner.v_aligner, SegmentGeneAligner):
            raise ValueError(MULTIPLE_DOMAINS_NOT_ENABLED)
        query_sequence = query_sequence.upper().replace("U", "T")
        # produce alignments for v,d,j,c
        alignments = self.vdjc_aligner.produce_nt_alignments(query_sequence)
        results = self._build_airr_entries(
            header, query_sequence, alignments, scheme, extend_alignment, return_all_domains, c_scheme
        )

        if not results:
            results = [
                (
                    AirrBuilder(header, "", scheme).build()
                    if not return_all_domains
                    else SegmentedAirrBuilder(header, "", scheme, query_sequence).build()
                )
            ]

        if return_all_domains:
            return results

        return sorted(results)[-1]

    def annotate_chain(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        c_scheme: Optional[CScheme] = None,
    ) -> ChainAnnotation[SegmentedAirrRearrangementEntryNT]:
        """
        Annotate all domains of a nucleotide chain in sequence order: rearrangements (V domains, as returned with
        return_all_domains=True) and constant domains. Constant domains are searched for in all parts of the chain not
        covered by V domains, in all reading frames (and both strands, if the chain has no V domains).
        Requires numbering created with return_all_domains=True.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Nucleotide sequence to annotate.
        - scheme : Numbering scheme of V domains (default IMGT).
        - extend_alignment: Extend alignments to cover full sequence.
        - c_scheme : Numbering scheme of constant domains (eu, kabat or imgt). By default constant domains are not
        numbered.

        Returns:
        -----------
        ChainAnnotation with rearrangements and constant domains. Constant domain positions are nucleotide positions
        in the chain sequence, which is reverse complemented if the chain is on the reverse strand.
        """
        if not isinstance(self.vdjc_aligner.v_aligner, SegmentGeneAligner) or not self.constant_region_aligner:
            raise ValueError(CHAIN_ANNOTATION_NOT_ENABLED)
        query_sequence = query_sequence.upper().replace("U", "T")

        alignments = [alignment for alignment in self.vdjc_aligner.produce_nt_alignments(query_sequence) if alignment.v]
        rearrangements = self._build_airr_entries(
            header, query_sequence, alignments, scheme, extend_alignment, True, c_scheme
        )
        v_domains = [v_domain_span_nt(alignment) for alignment in alignments]

        if alignments:
            assert alignments[0].v
            rev_comp = alignments[0].v.rev_comp
            sequence = reverse_complement(query_sequence) if rev_comp else query_sequence
            domain_alignments = align_constant_domains_nt(self.constant_region_aligner, sequence, v_domains)
        else:
            # without V domains the chain is on the strand with better constant domain alignments
            rev_comp, sequence, domain_alignments = max(
                (
                    (rev_comp, sequence, align_constant_domains_nt(self.constant_region_aligner, sequence, []))
                    for rev_comp, sequence in ((False, query_sequence), (True, reverse_complement(query_sequence)))
                ),
                key=lambda strand: sum(
                    domain_alignment.domain_alignment.alignment.alignment_score or 0 for domain_alignment in strand[2]
                ),
            )

        return self._chain_annotation(
            header,
            sequence,
            rev_comp,
            cast(list[SegmentedAirrRearrangementEntryNT], rearrangements),
            v_domains,
            domain_alignments,
            c_scheme,
        )

    def _build_airr_entries(
        self,
        header: str,
        query_sequence: str,
        alignments: list[AlignmentsNT],
        scheme: Scheme,
        extend_alignment: bool,
        return_all_domains: bool,
        c_scheme: Optional[CScheme],
    ) -> list[AirrRearrangementEntryNT]:
        return [
            self._build_airr_entry(
                header, query_sequence, alignment, scheme, extend_alignment, return_all_domains, c_scheme
            )
            for alignment in alignments
        ]

    def _build_airr_entry(
        self,
        header: str,
        query_sequence: str,
        alignment: AlignmentsNT,
        scheme: Scheme,
        extend_alignment: bool,
        return_all_domains: bool,
        c_scheme: Optional[CScheme],
    ) -> AirrRearrangementEntryNT:
        """AIRR record of one domain of the query, every domain is numbered separately."""
        sequence = ""
        segment_start = -1
        segment_end = -1
        translated_alignments = None
        nt_offsets = None
        aa_offsets = None
        sch_alignment = None
        scheme_residue_mapping = None
        positional_scheme_mapping = None
        c_numbering = None
        scheme_alignment_exc = None

        # Query might have been reverse complemented during V alignment.
        if alignment.v:
            sequence = alignment.v.q_seq
            segment_start = alignment.segment_start or 0
            segment_end = alignment.segment_end or len(sequence)
            try:
                translated_alignments = self.vjc_alignment_translator_aa.translate_nt_alignments(
                    alignment, extend_alignment=extend_alignment
                )

                if translated_alignments:
                    sch_alignment = self.scheme_aligner.align_to_scheme(
                        translated_alignments.aa_alignments, scheme, extend_alignment=False
                    )

                    # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
                    assert sch_alignment
                    assert alignment.v
                    chain_type = ChainType.from_locus(alignment.v.locus)

                    region_offsets = infer_region_offsets(sch_alignment, scheme, chain_type)
                    aa_offsets = infer_aa_region_offsets(region_offsets)

                    alignment_start = alignment.v.q_start

                    nt_offsets = infer_nt_region_offsets(
                        region_offsets, alignment_start, translated_alignments.reading_frame
                    )
                    numbering = produce_numbering(
                        translated_alignments.translated_query, sch_alignment.alignment_str, sch_alignment.q_start
                    )

                    if scheme == Scheme.IMGT:
                        numbering = fix_imgt_cdrs_numbering(numbering)
                        numbering = sort_imgt_numbering(numbering)

                    if numbering:
                        scheme_residue_mapping = get_scheme_residue_mapping(numbering)
                        positional_scheme_mapping = get_positional_scheme_mapping(
                            scheme_residue_mapping, sch_alignment.q_start
                        )

                    if c_scheme:
                        c_numbering = number_c_region(
                            self.scheme_aligner,
                            translated_alignments.translated_query,
                            translated_alignments.aa_alignments,
                            c_scheme,
                        )
            except AssertionError as err:
                scheme_alignment_exc = str(err)

        # produce airr result
        # this is separated on purpose to enable other output formats
        airr_builder = AirrBuilder(header, sequence or query_sequence, scheme)
        if return_all_domains:
            airr_builder = SegmentedAirrBuilder(header, sequence, scheme, query_sequence)
            airr_builder.with_segment_start_end(segment_start, segment_end)

        if alignment.v:
            airr_builder.with_v_gene_alignment(alignment.v)

        if translated_alignments:
            airr_builder.with_aa_alignments(translated_alignments)

        if alignment.j:
            airr_builder.with_j_gene_alignment(alignment.j)

        if alignment.d:
            airr_builder.with_d_gene_alignment(alignment.d)

        if alignment.c:
            airr_builder.with_c_gene_alignment(alignment.c)

        if sch_alignment:
            assert nt_offsets
            assert aa_offsets
            airr_builder.with_nt_region_offsets(nt_offsets)
            airr_builder.with_aa_region_offsets(aa_offsets)
            airr_builder.with_aa_scheme_alignment(sch_alignment)

        if scheme_residue_mapping is not None and positional_scheme_mapping is not None:
            airr_builder.with_scheme_residue_mapping(scheme_residue_mapping).with_positional_scheme_mapping(
                positional_scheme_mapping
            )

        if c_numbering:
            assert c_scheme
            airr_builder.with_c_scheme_numbering(c_scheme, *c_numbering)

        if scheme_alignment_exc:
            airr_builder.with_scheme_alignment_exc(scheme_alignment_exc)

        return airr_builder.build()


class RiotNumberingAA(_RiotNumbering[SegmentedAirrRearrangementEntryAA]):
    residue_length = 1

    def __init__(self, vjc_aligner_aa: VJCAlignerAA, scheme_aligner: SchemeAligner):
        self.vjc_aligner_aa = vjc_aligner_aa
        self.scheme_aligner = scheme_aligner
        self.constant_region_aligner = create_constant_region_aligner(vjc_aligner_aa.c_aligners, scheme_aligner)

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        *,
        v_gene: Optional[str] = ...,
        j_gene: Optional[str] = ...,
        c_scheme: Optional[CScheme] = ...,
    ) -> AirrRearrangementEntryAA: ...

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        return_all_domains: Literal[True] = ...,
        *,
        v_gene: Optional[str] = ...,
        j_gene: Optional[str] = ...,
        c_scheme: Optional[CScheme] = ...,
    ) -> list[AirrRearrangementEntryAA]: ...

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        return_all_domains: Literal[False] = ...,
        *,
        v_gene: Optional[str] = ...,
        j_gene: Optional[str] = ...,
        c_scheme: Optional[CScheme] = ...,
    ) -> AirrRearrangementEntryAA: ...

    @overload
    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = ...,
        extend_alignment: bool = ...,
        return_all_domains: bool = ...,
        *,
        v_gene: Optional[str] = ...,
        j_gene: Optional[str] = ...,
        c_scheme: Optional[CScheme] = ...,
    ) -> AirrRearrangementEntryAA | list[AirrRearrangementEntryAA]: ...

    def run_on_sequence(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        return_all_domains: bool = False,
        v_gene: Optional[str] = None,
        j_gene: Optional[str] = None,
        c_scheme: Optional[CScheme] = None,
    ) -> AirrRearrangementEntryAA | list[AirrRearrangementEntryAA]:  # pylint: disable=too-many-statements
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
        - return_all_domains: If True, return all domains of multiple domain proteins. Requires numbering created
        with return_all_domains=True, which returns the best scoring domain only if return_all_domains is False.
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
        if return_all_domains and not isinstance(self.vjc_aligner_aa.v_aligner, SegmentGeneAlignerAA):
            raise ValueError(MULTIPLE_DOMAINS_NOT_ENABLED)
        query_sequence = query_sequence.upper()
        # produce alignments for v,j,c
        if v_gene and j_gene:
            alignments = self.vjc_aligner_aa.produce_aa_alignments_for_genes(
                query_sequence, v_gene=v_gene, j_gene=j_gene
            )
        else:
            alignments = self.vjc_aligner_aa.produce_aa_alignments(query_sequence)

        results = self._build_airr_entries(
            header, query_sequence, alignments, scheme, extend_alignment, return_all_domains, c_scheme
        )

        if not results:
            results = [
                (
                    AirrBuilderAA(header, query_sequence, scheme).build()
                    if not return_all_domains
                    else SegmentedAirrBuilderAA(header, query_sequence, scheme, query_sequence).build()
                )
            ]

        if return_all_domains:
            return results

        return sorted(results)[-1]

    def annotate_chain(
        self,
        header: str,
        query_sequence: str,
        scheme: Scheme = Scheme.IMGT,
        extend_alignment: bool = True,
        c_scheme: Optional[CScheme] = None,
    ) -> ChainAnnotation[SegmentedAirrRearrangementEntryAA]:
        """
        Annotate all domains of an amino acid chain in sequence order: rearrangements (V domains, as returned with
        return_all_domains=True) and constant domains. Constant domains are searched for in all parts of the chain
        not covered by V domains. Requires numbering created with return_all_domains=True.

        Parameters:
        -----------
        - header : Fasta sequence header, usually this is simply sequence ID.
        - query_sequence : Amino acid sequence to annotate.
        - scheme : Numbering scheme of V domains (default IMGT).
        - extend_alignment: Extend alignments of V domains to cover full sequence.
        - c_scheme : Numbering scheme of constant domains (eu, kabat or imgt). By default constant domains are not
        numbered.

        Returns:
        -----------
        ChainAnnotation with rearrangements and constant domains.
        """
        if not isinstance(self.vjc_aligner_aa.v_aligner, SegmentGeneAlignerAA) or not self.constant_region_aligner:
            raise ValueError(CHAIN_ANNOTATION_NOT_ENABLED)
        query_sequence = query_sequence.upper()

        alignments = [
            alignment for alignment in self.vjc_aligner_aa.produce_aa_alignments(query_sequence) if alignment.v
        ]
        rearrangements = self._build_airr_entries(
            header, query_sequence, alignments, scheme, extend_alignment, True, c_scheme
        )
        v_domains = [v_domain_span_aa(alignment) for alignment in alignments]
        domain_alignments = align_constant_domains_aa(self.constant_region_aligner, query_sequence, v_domains)
        return self._chain_annotation(
            header,
            query_sequence,
            False,
            cast(list[SegmentedAirrRearrangementEntryAA], rearrangements),
            v_domains,
            domain_alignments,
            c_scheme,
        )

    def _build_airr_entries(
        self,
        header: str,
        query_sequence: str,
        alignments: list[AlignmentsAA],
        scheme: Scheme,
        extend_alignment: bool,
        return_all_domains: bool,
        c_scheme: Optional[CScheme],
    ) -> list[AirrRearrangementEntryAA]:
        return [
            self._build_airr_entry(
                header, query_sequence, alignment, scheme, extend_alignment, return_all_domains, c_scheme
            )
            for alignment in alignments
        ]

    def _build_airr_entry(
        self,
        header: str,
        query_sequence: str,
        alignment: AlignmentsAA,
        scheme: Scheme,
        extend_alignment: bool,
        return_all_domains: bool,
        c_scheme: Optional[CScheme],
    ) -> AirrRearrangementEntryAA:
        """AIRR record of one domain of the query, every domain is numbered separately."""
        sequence_aa = query_sequence
        sch_alignment = None
        aa_offsets = None
        scheme_residue_mapping = None
        positional_scheme_mapping = None
        c_numbering = None
        scheme_alignment_exc = None
        segment_start = -1
        segment_end = -1

        if alignment.v:
            segment_start = alignment.segment_start or 0
            segment_end = alignment.segment_end or len(alignment.v.q_seq)

            sequence_aa = alignment.v.q_seq
            try:
                sch_alignment = self.scheme_aligner.align_to_scheme(
                    alignment, scheme, extend_alignment=extend_alignment
                )

                # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
                assert sch_alignment
                chain_type = ChainType.from_locus(alignment.v.locus)

                region_offsets = infer_region_offsets(sch_alignment, scheme, chain_type)
                aa_offsets = infer_aa_region_offsets(region_offsets)

                numbering = produce_numbering(alignment.v.q_seq, sch_alignment.alignment_str, sch_alignment.q_start)

                if scheme == Scheme.IMGT:
                    numbering = fix_imgt_cdrs_numbering(numbering)
                    numbering = sort_imgt_numbering(numbering)

                if numbering:
                    scheme_residue_mapping = get_scheme_residue_mapping(numbering)
                    positional_scheme_mapping = get_positional_scheme_mapping(
                        scheme_residue_mapping, sch_alignment.q_start
                    )
            except AssertionError as err:
                scheme_alignment_exc = str(err)

            if c_scheme:
                c_numbering = number_c_region(self.scheme_aligner, alignment.v.q_seq, alignment, c_scheme)

        # produce airr result
        # this is separated on purpose to enable other output formats
        airr_builder = AirrBuilderAA(header, query_sequence, scheme)
        if return_all_domains:
            airr_builder = SegmentedAirrBuilderAA(header, sequence_aa, scheme, query_sequence)
            airr_builder.with_segment_start_end(segment_start, segment_end)

        if alignment.v:
            airr_builder.with_v_gene_alignment_aa(alignment.v)

        if alignment.j:
            airr_builder.with_j_gene_alignment_aa(alignment.j)

        if alignment.c:
            airr_builder.with_c_gene_alignment_aa(alignment.c)

        if sch_alignment:
            assert aa_offsets
            airr_builder.with_aa_region_offsets(aa_offsets)
            airr_builder.with_aa_scheme_alignment(sch_alignment, extend_alignment=extend_alignment)

        if scheme_residue_mapping is not None and positional_scheme_mapping is not None:
            airr_builder.with_scheme_residue_mapping(scheme_residue_mapping).with_positional_scheme_mapping(
                positional_scheme_mapping
            )

        if c_numbering:
            assert c_scheme
            airr_builder.with_c_scheme_numbering(c_scheme, *c_numbering)

        if scheme_alignment_exc:
            airr_builder.with_scheme_alignment_exc(scheme_alignment_exc)

        return airr_builder.build()


def create_riot_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    return_all_domains: bool = False,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingNT:
    """
    Factory function for creating RiotNumberingNT object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism tuple (default use all species - mouse, human and alpaca).
    - return_all_domains: If True, return all domains of multiple domain proteins.
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    RiotNumberingNT object used for performing alignment and numbering.
    """
    vdjc_aligner_nt = create_vdjc_aligner_nt(
        allowed_species=allowed_species, db_dir=db_dir, use_segment_aligner=return_all_domains
    )
    vjc_alignment_translator_aa = create_vjc_alignment_translator_aa(allowed_species=allowed_species, db_dir=db_dir)
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingNT(vdjc_aligner_nt, vjc_alignment_translator_aa, scheme_aligner)


def create_riot_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    return_all_domains: bool = False,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingAA:
    """
    Factory function for creating RiotNumberingAA object.

    Parameters:
    -----------
    - allowed_species : Limit gene database to specified Organism tuple (default use all species - mouse, human and alpaca).
    - return_all_domains: If True, return all domains of multiple domain proteins.
    - db_dir : Path to gene and scheme mappings database directory (default use embeded database).

    Returns:
    -----------
    RiotNumberingAA object used for performing alignment and numbering.
    """
    vjc_aligner_aa = create_vjc_aligner_aa(
        allowed_species=allowed_species,
        db_dir=db_dir,
        use_segment_aligner=return_all_domains,
    )
    scheme_aligner = SchemeAligner(allowed_species=allowed_species, db_dir=db_dir)
    return RiotNumberingAA(vjc_aligner_aa, scheme_aligner)


@cached({})
def get_or_create_riot_nt(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    return_all_domains: bool = False,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingNT:
    return create_riot_nt(allowed_species=allowed_species, return_all_domains=return_all_domains, db_dir=db_dir)


@cached({})
def get_or_create_riot_aa(
    allowed_species: Optional[tuple[Organism, ...]] = None,
    return_all_domains: bool = False,
    db_dir: Path = GENE_DB_DIR,
) -> RiotNumberingAA:
    return create_riot_aa(allowed_species=allowed_species, return_all_domains=return_all_domains, db_dir=db_dir)


if __name__ == "__main__":
    import json
    from dataclasses import asdict

    # nt_numbering = create_riot_nt(return_all_domains=True)
    # NT_QUERY = "caagugcaauuggucgagagcgggggcggggucgugcaaccgggcaggagccugcgcuugagcugugccgcgagccaguucacguucgggagcuauggcaugcacugggucaggcaaauacccgggaaagggcucgaguggguggccacgaucucguacgauggcacgaaaaaguaccaugccgauagcgugugggaucgcuuuaucauaagcagggacaauagcaagaacacgcucuuucuccagaugaauagcuugcgcccggaagauaccgcgcucuauuucugugucaaggaccagcgccaggacgagugcgaggaaugguggucggacuauuacgauuucgggcgcaggcucccgugcaggaaaucgcgcgggcuagcgggcauauucgaugucuggggccauggcacgauggugaccgugagcucg"
    # NT_QUERY = "CAGGCCCACCTGGAGCAAAGCGGCTCCGGGGTGAAGAAACCCGGAGCTTCTGTCAGAGTTAGCTGCTGGTCCTCTGAAGACATCTTCGAGCGGACCGAACTCATTCATTGGGTGCGCCAGGCCCCTGGCCAGGGGCTGGAGTGGATCGGATGGGTGAAGACAGTCACGGGCGCGGTGAACTTTGGGAGCCTTGATTTCAGACACCGGATTTCCCTGACCCGCGACAGAGATCTCTCTACAGCTTACATGGACATCCGGGGACTGACCCAGGATGACACAGCCACCTATTTTTGTGCCCGCCAAAAATTCGCTAGCAGATACTCCGGCGATCAGGGCAGCTATTTTGACCTTTGGGGGCGGGGAACACTGATTGTTGTGTCTTCC"
    # NT_QUERY = (
    #         # VH domain (nucleotide)
    #         "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTCCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTTAGCAGCTATGCCATGAGCTGGGTCCGCCAGGCTCCAGGCAAGGGGCTGGAGTGGGTGGCAGTTATATCATATGATGGAAGTAATAAATACTATGCAGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGATCGGGGGAGGAACGGTTATGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG"
    #         # (GGGGS)3 linker (nucleotide)
    #         "GGTGGCGGTTCAGGCGGAGGTGGCTCTGGCGGTGGCGGATCG"
    #         # VL domain (nucleotide)
    #         "GACATCCAGATGACCCAGTCTCCATCCTCCCTGTCTGCATCTGTAGGAGACAGAGTCACCATCACTTGCCGGGCAAGTCAGGGCATTAGCAGTTGGCTGGCCTGGTATCAGCAGAAACCAGGGAAAGCCCCTAAGCTCCTGATCTATGCTGCATCCAGTTTGCAAAGTGGGGTCCCATCAAGGTTCAGTGGCAGTGGATCTGGGACAGATTTCACTCTCACCATCAGCAGTCTGCAACCTGAAGATTTTGCAACTTACTACTGTCAGCAAGCTAACAGCTTCCCTTATACGTTCGGCCAAGGGACCAAGGTGGAAATCAAACGG"
    #     )
    # nt_sample_result = nt_numbering.run_on_sequence("header", NT_QUERY, Scheme.IMGT)
    # if isinstance(nt_sample_result, list):
    #     print(json.dumps([asdict(result) for result in nt_sample_result], indent=4))
    # else:
    #     print(json.dumps(asdict(nt_sample_result), indent=4))

    aa_numbering = create_riot_aa(return_all_domains=True)
    AA_QUERY = "EVQLVESGGGLVQPGGSLRLSCAASGSILSILSMAWYRQAPGKGRELVAGILSDGSTSYADSVKGRFTISRDNAKNTLYLQMNSLRPEDTAVYYCAALVDTENGAKLYWGQGTLVTVSS"

    aa_sample_result = aa_numbering.run_on_sequence("header", AA_QUERY, Scheme.IMGT, extend_alignment=True)
    # # Handle both single result and list result
    if isinstance(aa_sample_result, list):
        print(json.dumps([asdict(result) for result in aa_sample_result], indent=4))
    else:
        print(json.dumps(asdict(aa_sample_result), indent=4))
