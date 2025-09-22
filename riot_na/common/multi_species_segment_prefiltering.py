from typing import Sequence

from riot_na import Prefiltering
from riot_na.data.model import (
    Gene,
    GeneAA,
    Locus,
    Organism,
    SegmentMatch,
    SpeciesGeneMatch,
    SpeciesPrefilteringSegmentResult,
)
from riot_na.riot_na import (  # pylint: disable=no-name-in-module
    PrefilteringSegmentResult,
)


class MultiSpeciesSegmentPrefiltering:
    def __init__(
        self,
        all_genes: Sequence[Gene | GeneAA],
        kmer_size: int,
        distance_threshold: int,
        top_n: int,
        modulo_n: int,
        min_segment_length: int = 30,
        min_coverage: int = 20,
    ):
        self.all_genes = dict(
            map(lambda gene: (gene.species + "|" + gene.locus + "|" + gene.name, gene.sequence), all_genes)
        )

        self.prefiltering = Prefiltering(
            self.all_genes,
            top_n=top_n,
            kmer_size=kmer_size,
            distance_threshold=distance_threshold,
            modulo_n=modulo_n,
            min_segment_length=min_segment_length,
            min_coverage=min_coverage,
        )

    def calculate_segment_matches_with_rev_comp(self, query: str) -> SpeciesPrefilteringSegmentResult:
        try:
            raw_result = self.prefiltering.calculate_segment_matches_with_rev_comp(query)
        except ValueError:
            return SpeciesPrefilteringSegmentResult(query=query, rev_comp_query="", segments=[])

        return self._map_raw_result_to_segments(raw_result)

    def calculate_segment_matches(self, query: str) -> SpeciesPrefilteringSegmentResult:
        try:
            raw_result = self.prefiltering.calculate_segment_matches(query)
        except ValueError:
            return SpeciesPrefilteringSegmentResult(query=query, rev_comp_query="", segments=[])

        return self._map_raw_result_to_segments(raw_result)

    def _map_raw_result_to_segments(self, raw_result: PrefilteringSegmentResult) -> SpeciesPrefilteringSegmentResult:
        segments = []
        for segment in raw_result.segments:
            mapped_matches = []
            for match in segment.matching_genes:
                (species, locus, name) = match.gene_id.split("|")
                mapped_matches.append(
                    SpeciesGeneMatch(
                        species_gene_id=match.gene_id,
                        gene_id=name,
                        rev_comp=match.rev_comp,
                        coverage=match.coverage,
                        species=Organism(species),
                        locus=Locus(locus),
                    )
                )
            segments.append(
                SegmentMatch(
                    segment_start=segment.segment_start,
                    segment_length=segment.query_length(),
                    query_start=segment.query_start,
                    query_end=segment.query_end,
                    coverage=segment.coverage,
                    match_count=segment.match_count,
                    matching_genes=mapped_matches,
                )
            )
        return SpeciesPrefilteringSegmentResult(
            query=raw_result.query, rev_comp_query=raw_result.rev_comp_query, segments=segments
        )
