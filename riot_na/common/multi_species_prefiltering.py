from riot_na import Prefiltering
from riot_na.data.model import (
    Gene,
    Locus,
    Organism,
    SpeciesGeneMatch,
    SpeciesPrefilteringResult,
)


class MultiSpeciesPrefiltering:
    def __init__(self, all_genes: list[Gene], kmer_size: int, distance_threshold: int, top_n: int, modulo_n: int):
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
        mapped_matches = []
        try:
            raw_result = self.prefiltering.calculate_top_matches_with_rev_comp(query)
        except ValueError:
            return SpeciesPrefilteringResult(query=query, rev_comp_query="", top_matches=[])

        for match in raw_result.top_matches:
            (species, locus, name) = match.gene_id.split("|")
            mapped_matches.append(
                SpeciesGeneMatch(
                    gene_id=name,
                    rev_comp=match.rev_comp,
                    coverage=match.coverage,
                    species=Organism(species),
                    locus=Locus(locus),
                )
            )
        return SpeciesPrefilteringResult(
            query=query, rev_comp_query=raw_result.rev_comp_query, top_matches=mapped_matches
        )

    def calculate_top_matches(self, query: str) -> SpeciesPrefilteringResult:
        mapped_matches = []
        try:
            raw_result = self.prefiltering.calculate_top_matches(query)

        except ValueError:
            return SpeciesPrefilteringResult(query=query, rev_comp_query="", top_matches=[])

        for match in raw_result.top_matches:
            (species, locus, name) = match.gene_id.split("|")
            mapped_matches.append(
                SpeciesGeneMatch(
                    gene_id=name,
                    rev_comp=match.rev_comp,
                    coverage=match.coverage,
                    species=Organism(species),
                    locus=Locus(locus),
                )
            )

        return SpeciesPrefilteringResult(
            query=query, rev_comp_query=raw_result.rev_comp_query, top_matches=mapped_matches
        )
