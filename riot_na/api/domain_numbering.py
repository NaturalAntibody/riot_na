from typing import Optional

from riot_na.airr.airr_builder import AirrBuilder
from riot_na.airr.airr_builder import RearrangementT as RearrangementNTT
from riot_na.airr.airr_builder_aa import AirrBuilderAA
from riot_na.airr.airr_builder_aa import RearrangementT as RearrangementAAT
from riot_na.alignment.aa_gene_alignments import VJCAlignmentTranslatorAA
from riot_na.data.model import AlignmentsAA, AlignmentsNT, ChainType, CScheme, Scheme
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
    produce_c_numbering,
    produce_numbering,
    sort_imgt_numbering,
)


def number_c_region(
    scheme_aligner: SchemeAligner, query: str, aa_alignments: AlignmentsAA, c_scheme: CScheme
) -> Optional[tuple[dict[str, str], dict[int, str]]]:
    """
    Number C region of the query. Returns scheme residue mapping and positional scheme mapping
    or None if the C gene has no scheme mapping.
    """
    c_sch_alignment = scheme_aligner.align_c_to_scheme(aa_alignments, c_scheme)
    if not c_sch_alignment:
        return None

    assert aa_alignments.c
    chain_type = ChainType.from_locus(aa_alignments.c.locus)

    c_scheme_residue_mapping = produce_c_numbering(query, c_sch_alignment, c_scheme, chain_type)
    if not c_scheme_residue_mapping:
        return None

    c_positional_scheme_mapping = get_positional_scheme_mapping(c_scheme_residue_mapping, c_sch_alignment.q_start)
    return c_scheme_residue_mapping, c_positional_scheme_mapping


def number_domain_aa(  # pylint: disable=too-many-arguments
    builder: AirrBuilderAA[RearrangementAAT],
    alignments: AlignmentsAA,
    scheme_aligner: SchemeAligner,
    scheme: Scheme,
    extend_alignment: bool,
    c_scheme: Optional[CScheme],
) -> RearrangementAAT:
    """
    Number a single domain aligned to V, J and C germlines and build its AIRR rearrangement.
    Scheme alignment errors are reported in the rearrangement exc field.
    """
    if alignments.v is None:
        return builder.build()

    builder.with_v_gene_alignment_aa(alignments.v)
    if alignments.j:
        builder.with_j_gene_alignment_aa(alignments.j)
    if alignments.c:
        builder.with_c_gene_alignment_aa(alignments.c)

    try:
        sch_alignment = scheme_aligner.align_to_scheme(alignments, scheme, extend_alignment=extend_alignment)
        assert sch_alignment

        # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
        chain_type = ChainType.from_locus(alignments.v.locus)
        region_offsets = infer_region_offsets(sch_alignment, scheme, chain_type)
        builder.with_aa_region_offsets(infer_aa_region_offsets(region_offsets))
        builder.with_aa_scheme_alignment(sch_alignment, extend_alignment=extend_alignment)

        numbering = produce_numbering(alignments.v.q_seq, sch_alignment.alignment_str, sch_alignment.q_start)
        if scheme == Scheme.IMGT:
            numbering = fix_imgt_cdrs_numbering(numbering)
            numbering = sort_imgt_numbering(numbering)

        if numbering:
            scheme_residue_mapping = get_scheme_residue_mapping(numbering)
            builder.with_scheme_residue_mapping(scheme_residue_mapping).with_positional_scheme_mapping(
                get_positional_scheme_mapping(scheme_residue_mapping, sch_alignment.q_start)
            )
    except AssertionError as err:
        builder.with_scheme_alignment_exc(str(err))

    if c_scheme:
        c_numbering = number_c_region(scheme_aligner, alignments.v.q_seq, alignments, c_scheme)
        if c_numbering:
            builder.with_c_scheme_numbering(c_scheme, *c_numbering)

    return builder.build()


def number_domain_nt(  # pylint: disable=too-many-arguments
    builder: AirrBuilder[RearrangementNTT],
    alignments: AlignmentsNT,
    alignment_translator: VJCAlignmentTranslatorAA,
    scheme_aligner: SchemeAligner,
    scheme: Scheme,
    extend_alignment: bool,
    c_scheme: Optional[CScheme],
) -> RearrangementNTT:
    """
    Translate and number a single domain aligned to V, D, J and C germlines and build its AIRR rearrangement.
    Translation and scheme alignment errors are reported in the rearrangement exc field.
    """
    if alignments.v is None:
        return builder.build()

    builder.with_v_gene_alignment(alignments.v)
    if alignments.j:
        builder.with_j_gene_alignment(alignments.j)
    if alignments.d:
        builder.with_d_gene_alignment(alignments.d)
    if alignments.c:
        builder.with_c_gene_alignment(alignments.c)

    try:
        translated_alignments = alignment_translator.translate_nt_alignments(
            alignments, extend_alignment=extend_alignment
        )
        if translated_alignments is not None:
            builder.with_aa_alignments(translated_alignments)

            sch_alignment = scheme_aligner.align_to_scheme(
                translated_alignments.aa_alignments, scheme, extend_alignment=False
            )
            assert sch_alignment

            # calculate aa offsets based on the numbering (iterate over cigar to produce positional-> imgt mapping)
            chain_type = ChainType.from_locus(alignments.v.locus)
            region_offsets = infer_region_offsets(sch_alignment, scheme, chain_type)
            builder.with_nt_region_offsets(
                infer_nt_region_offsets(region_offsets, alignments.v.q_start, translated_alignments.reading_frame)
            )
            builder.with_aa_region_offsets(infer_aa_region_offsets(region_offsets))
            builder.with_aa_scheme_alignment(sch_alignment)

            numbering = produce_numbering(
                translated_alignments.translated_query, sch_alignment.alignment_str, sch_alignment.q_start
            )
            if scheme == Scheme.IMGT:
                numbering = fix_imgt_cdrs_numbering(numbering)
                numbering = sort_imgt_numbering(numbering)

            if numbering:
                scheme_residue_mapping = get_scheme_residue_mapping(numbering)
                builder.with_scheme_residue_mapping(scheme_residue_mapping).with_positional_scheme_mapping(
                    get_positional_scheme_mapping(scheme_residue_mapping, sch_alignment.q_start)
                )

            if c_scheme:
                c_numbering = number_c_region(
                    scheme_aligner,
                    translated_alignments.translated_query,
                    translated_alignments.aa_alignments,
                    c_scheme,
                )
                if c_numbering:
                    builder.with_c_scheme_numbering(c_scheme, *c_numbering)
    except AssertionError as err:
        builder.with_scheme_alignment_exc(str(err))

    return builder.build()
