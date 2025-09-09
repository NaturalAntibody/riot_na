from typing import Sequence

from riot_na import Prefiltering
from riot_na.data.model import (
    Gene,
    GeneAA,
    Locus,
    Organism,
    SpeciesGeneMatch,
    SpeciesPrefilteringResult,
)
from riot_na.riot_na import PrefilteringResult  # pylint: disable=no-name-in-module


class MultiSpeciesPrefiltering:
    def __init__(
        self, all_genes: Sequence[Gene | GeneAA], kmer_size: int, distance_threshold: int, top_n: int, modulo_n: int
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
        )

    def calculate_top_matches_with_rev_comp(self, query: str) -> SpeciesPrefilteringResult:

        try:
            raw_result = self.prefiltering.calculate_top_matches_with_rev_comp(query)
        except ValueError:
            return SpeciesPrefilteringResult(query=query, rev_comp_query="", top_matches=[])

        return self._map_raw_result_to_matches(raw_result)

    def calculate_top_matches(self, query: str) -> SpeciesPrefilteringResult:

        try:
            raw_result = self.prefiltering.calculate_top_matches(query)
        except ValueError:
            return SpeciesPrefilteringResult(query=query, rev_comp_query="", top_matches=[])

        return self._map_raw_result_to_matches(raw_result)

    def _map_raw_result_to_matches(self, raw_result: PrefilteringResult) -> SpeciesPrefilteringResult:
        mapped_matches = []
        for match in raw_result.top_matches:
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
        return SpeciesPrefilteringResult(
            query=raw_result.query, rev_comp_query=raw_result.rev_comp_query, top_matches=mapped_matches
        )
