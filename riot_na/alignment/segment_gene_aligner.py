#!/usr/bin/env python3
"""
Segment-based Gene Aligner

This module provides segment-centric alignment functionality that leverages the new
segment-centric prefiltering to return alignments for each domain/segment separately.
This is particularly useful for:
- scFv sequences (separate VL and VH alignments)
- Fusion proteins (individual domain alignments)
- Multi-domain constructs
- Domain-specific alignment quality assessment
"""

import os
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import List, Optional, Sequence

from skbio.alignment import StripedSmithWaterman  # type: ignore

from riot_na.alignment.gene_aligner import AA_ALIGNER_PARAMS, ALIGNER_PARAMS
from riot_na.alignment.gene_aligner import (
    get_aa_aligner_params as gene_aa_aligner_params,
)
from riot_na.alignment.gene_aligner import get_aligner_params as gene_aligner_params
from riot_na.alignment.gene_aligner import (
    get_gene_parsing_function,
    parse_gene_aa,
    read_genes,
)
from riot_na.alignment.skbio_alignment import align, align_aa
from riot_na.common.gene_match_utils import create_gene_lookup
from riot_na.common.multi_species_segment_prefiltering import (
    MultiSpeciesSegmentPrefiltering,
)
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AlignmentEntryAA,
    AlignmentEntryNT,
    Gene,
    GeneAA,
    GermlineGene,
    Locus,
    Organism,
    SegmentMatch,
    SpeciesPrefilteringSegmentResult,
)


@dataclass
class SegmentAlignment:
    """Alignment result for a single segment/domain"""

    segment_start: int
    segment_end: int
    alignments: List[AlignmentEntryNT]  # Alignments for genes in this segment
    best_alignment: Optional[AlignmentEntryNT] = None  # Best alignment for this segment

    def __post_init__(self):
        """Set the best alignment after initialization"""
        if self.alignments:
            self.best_alignment = min(self.alignments)  # Lowest e-value first


@dataclass
class SegmentAlignmentAA:
    """Amino acid alignment result for a single segment/domain"""

    segment_start: int
    segment_end: int
    alignments: List[AlignmentEntryAA]  # Alignments for genes in this segment
    best_alignment: Optional[AlignmentEntryAA] = None  # Best alignment for this segment

    def __post_init__(self):
        """Set the best alignment after initialization"""
        if self.alignments:
            self.best_alignment = min(self.alignments)  # Lowest e-value first


@dataclass
class SegmentAlignmentResult:
    """Complete segment-based alignment result"""

    query: str
    rev_comp_query: str
    segment_alignments: List[SegmentAlignment]
    prefiltering_result: SpeciesPrefilteringSegmentResult

    @property
    def total_segments(self) -> int:
        """Number of segments detected"""
        return len(self.segment_alignments)

    @property
    def successful_alignments(self) -> int:
        """Number of segments with successful alignments"""
        return sum(1 for seg_align in self.segment_alignments if seg_align.best_alignment)


@dataclass
class SegmentAlignmentResultAA:
    """Complete segment-based amino acid alignment result"""

    query: str
    segment_alignments: List[SegmentAlignmentAA]
    prefiltering_result: SpeciesPrefilteringSegmentResult

    @property
    def total_segments(self) -> int:
        """Number of segments detected"""
        return len(self.segment_alignments)

    @property
    def successful_alignments(self) -> int:
        """Number of segments with successful alignments"""
        return sum(1 for seg_align in self.segment_alignments if seg_align.best_alignment)


class SegmentGeneAligner:
    """
    Segment-based gene aligner that provides alignments for each detected domain/segment.

    This aligner uses the segment-centric prefiltering approach to first identify
    non-overlapping segments in the query sequence, then performs detailed alignments
    for each segment against its best matching genes.

    Compatible with GeneAligner interface for easy replacement.
    """

    def __init__(
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
    ) -> None:
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
        self.alignment_length_threshold = alignment_length_threshold

    # Compatibility method to match GeneAligner interface
    def _prefilter(self, query: str, both_strains: bool) -> SpeciesPrefilteringSegmentResult:
        """Prefilter method compatible with GeneAligner interface"""
        if both_strains:
            return self.prefiltering.calculate_segment_matches_with_rev_comp(query)
        return self.prefiltering.calculate_segment_matches(query)

    def _align_segment_genes(self, segment: SegmentMatch, query: str, rev_comp_query: str) -> List[AlignmentEntryNT]:
        """Align the top genes for a specific segment"""
        alignments = []

        for gene_info in segment.matching_genes:
            # Get the appropriate query sequence based on strand
            current_query = rev_comp_query if gene_info.rev_comp else query

            # Create aligner for this query
            aligner = StripedSmithWaterman(current_query, **ALIGNER_PARAMS)

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

            if self.e_value_threshold is not None:
                if alignment.e_value > self.e_value_threshold:
                    continue

            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
                    continue

            # Convert to AlignmentEntryNT
            alignment_entry = AlignmentEntryNT(
                target_id=alignment.target_id,
                alignment_score=alignment.alignment_score,
                seq_identity=alignment.seq_identity,
                e_value=alignment.e_value,
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
    ) -> List[AlignmentEntryNT]:
        """Align sequences method compatible with GeneAligner interface"""
        all_alignments = []

        # Align each segment separately and collect all alignments
        for segment in prefiltering_result.segments:
            segment_alignments = self._align_segment_genes(segment, query, prefiltering_result.rev_comp_query)
            all_alignments.extend(segment_alignments)

        return all_alignments

    # Compatibility method to match GeneAligner interface
    def align(self, query: str, both_strains: bool = True) -> Optional[AlignmentEntryNT]:
        """
        Compatibility method that returns the best single alignment across all segments.
        This matches the GeneAligner interface for drop-in replacement.

        Args:
            query: Query sequence to align
            both_strains: Whether to consider both forward and reverse complement

        Returns:
            Best alignment across all segments, or None if no alignments found
        """
        prefiltering_result = self._prefilter(query, both_strains)
        alignments = self._align_sequences(query, prefiltering_result)

        if not alignments:
            return None

        # Sort alignments and return the best one (lowest e-value)
        alignments.sort()
        return alignments[0]

    def align_segments(self, query: str, both_strains: bool = True) -> SegmentAlignmentResult:
        """
        Perform segment-based alignment of the query sequence.

        Args:
            query: Query sequence to align
            both_strains: Whether to consider both forward and reverse complement

        Returns:
            SegmentAlignmentResult containing alignments for each detected segment
        """
        # Get segment-centric prefiltering results
        prefiltering_result = self._prefilter(query, both_strains)

        segment_alignments = []

        # Align each segment separately
        for segment_no, segment in enumerate(prefiltering_result.segments):
            current_segment_algn_start = segment.segment_start
            # has next segment
            if segment_no + 1 < len(prefiltering_result.segments):
                next_segment_algn_start = prefiltering_result.segments[segment_no + 1].segment_start
                assert (
                    next_segment_algn_start > current_segment_algn_start
                ), f"Next segment start is before current segment start: {next_segment_algn_start} < {current_segment_algn_start}"
            else:
                next_segment_algn_start = len(query)

            alignments = self._align_segment_genes(
                segment,
                query[current_segment_algn_start:next_segment_algn_start],
                prefiltering_result.rev_comp_query[current_segment_algn_start:next_segment_algn_start],
            )

            if alignments:
                segment_alignment = SegmentAlignment(
                    segment_start=current_segment_algn_start, segment_end=next_segment_algn_start, alignments=alignments
                )

                segment_alignments.append(segment_alignment)

        return SegmentAlignmentResult(
            query=query,
            rev_comp_query=prefiltering_result.rev_comp_query,
            segment_alignments=segment_alignments,
            prefiltering_result=prefiltering_result,
        )


class SegmentGeneAlignerAA:
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
    ) -> None:
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
        self.alignment_length_threshold = alignment_length_threshold
        self.max_cdr3_length = max_cdr3_length  # Added for compatibility

    # Compatibility method to match GeneAlignerAA interface
    def _prefilter(self, query: str) -> SpeciesPrefilteringSegmentResult:
        """Prefilter method compatible with GeneAlignerAA interface"""
        return self.prefiltering.calculate_segment_matches(query)

    def _align_segment_genes_aa(self, segment: SegmentMatch, query: str) -> List[AlignmentEntryAA]:
        """Align the top genes for a specific segment (amino acid version)"""
        alignments = []

        # Create aligner for this query
        aligner = StripedSmithWaterman(query, **AA_ALIGNER_PARAMS)

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

            if self.e_value_threshold is not None:
                if alignment.e_value and alignment.e_value > self.e_value_threshold:
                    continue

            if self.alignment_length_threshold is not None:
                if alignment_length < self.alignment_length_threshold:
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
                e_value=alignment.e_value,
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
    ) -> List[AlignmentEntryAA]:
        """Align sequences method compatible with GeneAlignerAA interface"""
        all_alignments = []

        # Align each segment separately and collect all alignments
        for segment in prefiltering_result.segments:
            segment_alignments = self._align_segment_genes_aa(segment, query)
            all_alignments.extend(segment_alignments)

        return all_alignments

    # Compatibility method to match GeneAlignerAA interface
    def align(self, query: str) -> Optional[AlignmentEntryAA]:
        """
        Compatibility method that returns the best single alignment across all segments.
        This matches the GeneAlignerAA interface for drop-in replacement.

        Args:
            query: Query AA sequence to align

        Returns:
            Best alignment across all segments, or None if no alignments found
        """
        prefiltering_result = self._prefilter(query)
        alignments = self._align_sequences(query, prefiltering_result)

        if not alignments:
            return None

        # Sort alignments and return the best one (lowest e-value)
        alignments.sort()
        return alignments[0]

    def align_segments(self, query: str) -> SegmentAlignmentResultAA:
        """
        Perform segment-based amino acid alignment of the query sequence.

        Args:
            query: Query AA sequence to align

        Returns:
            SegmentAlignmentResultAA containing alignments for each detected segment
        """
        # Get segment-centric prefiltering results (AA doesn't use rev comp)
        prefiltering_result = self._prefilter(query)

        segment_alignments = []
        # Align each segment separately
        for segment_no, segment in enumerate(prefiltering_result.segments):
            current_segment_algn_start = segment.segment_start
            if segment_no + 1 < len(prefiltering_result.segments):
                next_segment_algn_start = prefiltering_result.segments[segment_no + 1].segment_start
            else:
                next_segment_algn_start = len(query)

            alignments = self._align_segment_genes_aa(
                segment, query[current_segment_algn_start:next_segment_algn_start]
            )

            segment_alignment = SegmentAlignmentAA(
                segment_start=current_segment_algn_start, segment_end=next_segment_algn_start, alignments=alignments
            )
            segment_alignments.append(segment_alignment)

        return SegmentAlignmentResultAA(
            query=query, segment_alignments=segment_alignments, prefiltering_result=prefiltering_result
        )


def get_aligner_params(germline_gene: GermlineGene, locus: Optional[Locus]) -> dict:
    """Get aligner parameters for different germline gene types"""
    # ideally this should be read from a file or env

    match germline_gene:
        case GermlineGene.V:
            return {
                "kmer_size": 9,
                "distance_threshold": 13,
                "top_n": int(os.environ.get("TOP_N", 12)),
                "modulo_n": 2,
                "e_value_threshold": 0.001,
                "alignment_length_threshold": 80,
                "min_prefiltering_coverage": 60,
                "min_segment_length": 120,
            }
        case _:
            return gene_aligner_params(germline_gene, locus)


def get_aa_aligner_params(germline_gene: GermlineGene) -> dict:
    """Get amino acid aligner parameters for different germline gene types"""
    # ideally this should be read from a file or env

    match germline_gene:
        case GermlineGene.V:
            return {
                "kmer_size": 3,
                "distance_threshold": 4,
                "top_n": int(os.environ.get("TOP_N", 100)),
                "modulo_n": 1,
                "e_value_threshold": 1e-55,
                "alignment_length_threshold": 80,
                "min_prefiltering_coverage": 20,
                "min_segment_length": 60,
            }
        case _:
            return gene_aa_aligner_params(germline_gene)


def create_v_gene_aligner(allowed_species: List[Organism], db_dir: Path = GENE_DB_DIR) -> SegmentGeneAligner:
    """Create V gene segment aligner"""
    genes = []

    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS]

    gene_parsing_function = get_gene_parsing_function(GermlineGene.V)

    for species in allowed_species:
        input_path = Path(db_dir) / "gene_db" / "v_genes" / f"{species.value}.fasta"
        species_genes = read_genes(input_path, gene_parsing_function)
        genes.extend(species_genes)

    aligner_params = get_aligner_params(GermlineGene.V, None)
    genes_aligner = SegmentGeneAligner(genes, **aligner_params)
    return genes_aligner


def create_aa_v_gene_aligner(allowed_species: List[Organism], aa_genes_dir: Path) -> SegmentGeneAlignerAA:
    """Create amino acid V gene segment aligner"""
    genes = []

    if not allowed_species:
        allowed_species = [Organism.HOMO_SAPIENS, Organism.MUS_MUSCULUS]

    for species in allowed_species:
        input_path = aa_genes_dir / "v_genes" / f"{species.value}.fasta"

        assert input_path.exists(), f"Input germline file: {str(input_path)} does not exists"

        species_genes = read_genes(input_path, partial(parse_gene_aa, species=species))
        genes.extend(species_genes)

    aligner_params = get_aa_aligner_params(GermlineGene.V)
    genes_aligner = SegmentGeneAlignerAA(genes, **aligner_params)
    return genes_aligner
