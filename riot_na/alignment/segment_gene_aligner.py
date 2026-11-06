from functools import partial
from pathlib import Path
from typing import Callable, Optional, Sequence

from ssw_aligner import StripedSmithWaterman

from riot_na.alignment.alignment_metrics import log10_evalue_threshold
from riot_na.alignment.gene_aligner import (
    get_gene_parsing_function,
    parse_gene_aa,
    read_genes,
)
from riot_na.alignment.params import (
    DEFAULT_ALIGNER_PARAMS,
    AASSWParams,
    AlignerParams,
    NTSSWParams,
)
from riot_na.alignment.skbio_alignment import align, align_aa
from riot_na.common.gene_match_utils import create_gene_lookup
from riot_na.common.multi_species_segment_prefiltering import (
    MultiSpeciesSegmentPrefiltering,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntry,
    AlignmentEntryAA,
    AlignmentEntryNT,
    AlignmentSegment,
    Gene,
    GeneAA,
    GermlineGene,
    Organism,
    SegmentMatch,
    SpeciesPrefilteringSegmentResult,
)


class SegmentGeneAligner:  # pylint: disable=too-many-instance-attributes
    """
    Segment-based gene aligner that provides alignments for each detected domain/segment.

    This aligner uses the segment-centric prefiltering approach to first identify
    non-overlapping segments in the query sequence, then performs detailed alignments
    for each segment against its best matching genes.

    Compatible with GeneAligner interface for easy replacement.
    """

    def __init__(  # pylint: disable=too-many-arguments
        self,
        genes: Sequence[Gene],
        kmer_size: int = 6,
        distance_threshold: int = 10,
        top_n: int = 10,
        modulo_n: int = 1,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
        min_prefiltering_coverage: int = 30,
        min_segment_length: int = 30,
        min_seq_identity: Optional[float] = None,
        ssw_params: NTSSWParams = DEFAULT_ALIGNER_PARAMS.nt_ssw,
    ) -> None:
        self.ssw_params = ssw_params
        self.min_seq_identity = min_seq_identity
        # Initialize the Rust-based prefiltering
        self.prefiltering = MultiSpeciesSegmentPrefiltering(
            all_genes=genes,
            kmer_size=kmer_size,
            distance_threshold=distance_threshold,
            top_n=top_n,
            modulo_n=modulo_n,
            min_segment_length=min_segment_length,
            min_coverage=min_prefiltering_coverage,
        )

        # Store gene sequences for alignment
        self.gene_lookup = create_gene_lookup(genes)
        self.db_length = sum(len(gene.sequence) for gene in genes)

        # Alignment parameters
        self.e_value_threshold = e_value_threshold
        self.log10_e_value_threshold = log10_evalue_threshold(e_value_threshold)
        self.alignment_length_threshold = alignment_length_threshold

    # Compatibility method to match GeneAligner interface
    def _prefilter(self, query: str, both_strains: bool) -> SpeciesPrefilteringSegmentResult:
        """Prefilter method compatible with GeneAligner interface"""
        if both_strains:
            return self.prefiltering.calculate_segment_matches_with_rev_comp(query)
        return self.prefiltering.calculate_segment_matches(query)

    def _align_segment_genes(
        self, segment: SegmentMatch, query: str, rev_comp_query: str
    ) -> Sequence[AlignmentEntryNT]:
        """Align the top genes for a specific segment"""
        alignments = []

        for gene_info in segment.matching_genes:
            # Get the appropriate query sequence based on strand
            current_query = rev_comp_query if gene_info.rev_comp else query

            # Create aligner for this query
            aligner = StripedSmithWaterman(current_query, **self.ssw_params.kwargs())

            # Get target sequence
            target_gene = self.gene_lookup[gene_info.species_gene_id]

            # Perform alignment
            alignment = align(
                aligner=aligner,
                target_id=gene_info.gene_id,
                target=target_gene.sequence,
                query=current_query,
                db_length=self.db_length,
                rev_comp=gene_info.rev_comp,
            )

            # Apply filtering criteria
            alignment_length = alignment.t_end - alignment.t_start

            if self.log10_e_value_threshold is not None:
                if alignment.log10_e_value > self.log10_e_value_threshold:
                    continue

            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
                    continue

            if self.min_seq_identity is not None and alignment.seq_identity < self.min_seq_identity:
                continue

            # Convert to AlignmentEntryNT
            alignment_entry = AlignmentEntryNT(
                target_id=alignment.target_id,
                alignment_score=alignment.alignment_score,
                seq_identity=alignment.seq_identity,
                log10_e_value=alignment.log10_e_value,
                q_start=alignment.q_start,
                q_end=alignment.q_end,
                q_len=alignment.q_len,
                t_start=alignment.t_start,
                t_end=alignment.t_end,
                t_len=alignment.t_len,
                cigar=alignment.cigar,
                rev_comp=alignment.rev_comp,
                species=target_gene.species,
                locus=target_gene.locus,
                q_seq=current_query,
                t_seq=target_gene.sequence,
                reading_frame=target_gene.reading_frame if isinstance(target_gene, Gene) else None,
            )

            alignments.append(alignment_entry)

        return alignments

    # Compatibility method to match GeneAligner interface
    def _align_sequences(
        self, query: str, prefiltering_result: SpeciesPrefilteringSegmentResult
    ) -> Sequence[AlignmentEntryNT]:
        """Align sequences method compatible with GeneAligner interface"""
        all_alignments: list[AlignmentEntryNT] = []

        # Align each segment separately and collect all alignments
        for segment in prefiltering_result.segments:
            segment_alignments = self._align_segment_genes(segment, query, prefiltering_result.rev_comp_query)
            all_alignments.extend(segment_alignments)

        return all_alignments

    def align(self, query: str, both_strains: bool = True) -> Sequence[AlignmentSegment[AlignmentEntryNT]]:
        """
        Perform segment-based alignment of the query sequence.

        Args:
            query: Query sequence to align
            both_strains: Whether to consider both forward and reverse complement

        Returns:
            AlignmentResult containing alignments for each detected segment
        """
        # Get segment-centric prefiltering results
        prefiltering_result = self._prefilter(query, both_strains)

        return align_segments(
            prefiltering_result.segments,
            len(query),
            lambda segment, start, end: self._align_segment_genes(
                segment, query[start:end], prefiltering_result.rev_comp_query[start:end]
            ),
        )


class SegmentGeneAlignerAA:  # pylint: disable=too-many-instance-attributes
    """
    Amino acid version of the segment-based gene aligner.

    Compatible with GeneAlignerAA interface for easy replacement.
    """

    def __init__(  # pylint: disable=too-many-arguments
        self,
        genes: Sequence[GeneAA],  # gene_id -> AA sequence mapping
        kmer_size: int = 3,  # Smaller k-mer for amino acids
        distance_threshold: int = 5,
        top_n: int = 10,
        modulo_n: int = 1,
        e_value_threshold: Optional[float] = None,
        alignment_length_threshold: Optional[int] = None,
        max_cdr3_length: Optional[int] = None,  # Added for compatibility with GeneAlignerAA
        min_prefiltering_coverage: int = 20,
        min_segment_length: int = 60,
        min_seq_identity: Optional[float] = None,
        ssw_params: AASSWParams = DEFAULT_ALIGNER_PARAMS.aa_ssw,
    ) -> None:
        self.ssw_params = ssw_params
        self.min_seq_identity = min_seq_identity
        self.prefiltering = MultiSpeciesSegmentPrefiltering(
            genes,
            kmer_size,
            distance_threshold,
            top_n,
            modulo_n,
            min_segment_length=min_segment_length,
            min_coverage=min_prefiltering_coverage,
        )

        # Store gene sequences for alignment
        self.db_length = sum(len(gene.sequence) for gene in genes)
        self.gene_lookup = create_gene_lookup(genes)

        # Alignment parameters
        self.e_value_threshold = e_value_threshold
        self.log10_e_value_threshold = log10_evalue_threshold(e_value_threshold)
        self.alignment_length_threshold = alignment_length_threshold
        self.max_cdr3_length = max_cdr3_length  # Added for compatibility

    # Compatibility method to match GeneAlignerAA interface
    def _prefilter(self, query: str) -> SpeciesPrefilteringSegmentResult:
        """Prefilter method compatible with GeneAlignerAA interface"""
        return self.prefiltering.calculate_segment_matches(query)

    def _align_segment_genes_aa(self, segment: SegmentMatch, query: str) -> Sequence[AlignmentEntryAA]:
        """Align the top genes for a specific segment (amino acid version)"""
        alignments = []

        # Create aligner for this query
        aligner = StripedSmithWaterman(query, **self.ssw_params.kwargs())

        for gene_info in segment.matching_genes:
            # Note: AA version doesn't support reverse complement
            if gene_info.rev_comp:
                continue

            # Get target sequence
            target_sequence = self.gene_lookup[gene_info.species_gene_id].sequence

            # Perform alignment
            alignment = align_aa(
                aligner=aligner,
                target_id=gene_info.gene_id,
                target=target_sequence,
                query=query,
                db_length=self.db_length,
            )

            # Apply filtering criteria
            alignment_length = alignment.t_end - alignment.t_start

            if self.log10_e_value_threshold is not None:
                if alignment.log10_e_value is not None and alignment.log10_e_value > self.log10_e_value_threshold:
                    continue

            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
                    continue

            if (
                self.min_seq_identity is not None
                and alignment.seq_identity is not None
                and alignment.seq_identity < self.min_seq_identity
            ):
                continue

            # Apply max_cdr3_length filter for compatibility with GeneAlignerAA
            if self.max_cdr3_length is not None:
                if alignment.q_start > self.max_cdr3_length:
                    continue

            # Convert to AlignmentEntryAA
            alignment_entry = AlignmentEntryAA(
                target_id=alignment.target_id,
                alignment_score=alignment.alignment_score,
                seq_identity=alignment.seq_identity,
                log10_e_value=alignment.log10_e_value,
                q_start=alignment.q_start,
                q_end=alignment.q_end,
                t_start=alignment.t_start,
                t_end=alignment.t_end,
                cigar=alignment.cigar,
                species=gene_info.species,
                locus=gene_info.locus,
                q_seq=query,
                t_seq=target_sequence,
            )

            alignments.append(alignment_entry)

        return alignments

    # Compatibility method to match GeneAlignerAA interface
    def _align_sequences(
        self, query: str, prefiltering_result: SpeciesPrefilteringSegmentResult
    ) -> Sequence[AlignmentEntryAA]:
        """Align sequences method compatible with GeneAlignerAA interface"""
        all_alignments: list[AlignmentEntryAA] = []

        # Align each segment separately and collect all alignments
        for segment in prefiltering_result.segments:
            segment_alignments = self._align_segment_genes_aa(segment, query)
            all_alignments.extend(segment_alignments)

        return all_alignments

    def align(self, query: str) -> Sequence[AlignmentSegment[AlignmentEntryAA]]:
        """
        Perform segment-based amino acid alignment of the query sequence.

        Args:
            query: Query AA sequence to align

        Returns:
            List[AlignmentSegment] containing alignments for each detected segment
        """
        # Get segment-centric prefiltering results (AA doesn't use rev comp)
        prefiltering_result = self._prefilter(query)

        return align_segments(
            prefiltering_result.segments,
            len(query),
            lambda segment, start, end: self._align_segment_genes_aa(segment, query[start:end]),
        )


def align_segments(
    segments: Sequence[SegmentMatch],
    query_length: int,
    align_segment: Callable[[SegmentMatch, int, int], Sequence[AlignmentEntry]],
) -> list[AlignmentSegment[AlignmentEntry]]:
    """
    Align each segment from its start to the start of the next segment.

    Segments are aligned starting from the last one, so a segment without alignments (e.g. a spurious
    prefiltering match in a constant region) is dropped and the preceding segment extends over it.
    """
    segment_alignments: list[AlignmentSegment[AlignmentEntry]] = []
    segment_end = query_length

    for segment in reversed(segments):
        if segment.segment_start >= segment_end:
            continue

        alignments = align_segment(segment, segment.segment_start, segment_end)
        if alignments:
            segment_alignments.append(
                AlignmentSegment(segment_start=segment.segment_start, segment_end=segment_end, alignments=alignments)
            )
            segment_end = segment.segment_start

    return segment_alignments[::-1]


def create_v_gene_aligner(
    allowed_species: Sequence[Organism],
    db_dir: Path = GENE_DB_DIR,
    params: AlignerParams = DEFAULT_ALIGNER_PARAMS,
) -> SegmentGeneAligner:
    """Create V gene segment aligner"""
    genes = []

    if not allowed_species:
        allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

    gene_parsing_function = get_gene_parsing_function(GermlineGene.V)

    for species in allowed_species:
        input_path = Path(db_dir) / "gene_db" / "v_genes" / f"{species.value}.fasta"
        species_genes = read_genes(input_path, gene_parsing_function)
        genes.extend(species_genes)

    genes_aligner = SegmentGeneAligner(genes, **params.segment_v_nt.kwargs(), ssw_params=params.nt_ssw)
    return genes_aligner


def create_aa_v_gene_aligner(
    allowed_species: Sequence[Organism], aa_genes_dir: Path, params: AlignerParams = DEFAULT_ALIGNER_PARAMS
) -> SegmentGeneAlignerAA:
    """Create amino acid V gene segment aligner"""
    genes = []

    if not allowed_species:
        allowed_species = (Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS, Organism.VICUGNA_PACOS)

    for species in allowed_species:
        input_path = aa_genes_dir / "v_genes" / f"{species.value}.fasta"

        assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

        species_genes = read_genes(input_path, partial(parse_gene_aa, species=species))
        genes.extend(species_genes)

    genes_aligner = SegmentGeneAlignerAA(genes, **params.segment_v_aa.kwargs(), ssw_params=params.aa_ssw)
    return genes_aligner
