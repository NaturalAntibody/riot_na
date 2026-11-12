"""
Utility functions for handling gene matches and species information.

This module provides shared functionality for converting between different gene match formats
and extracting species information, avoiding the need for proxy classes and code duplication.
"""

from typing import Sequence

from riot_na.data.model import Gene, GeneAA, Locus, Organism


def gene_lookup_key(species: Organism, locus: Locus, gene_name: str) -> str:
    return f"{species.value}|{locus.value}|{gene_name}"


def create_gene_lookup(genes: Sequence[Gene | GeneAA]) -> dict[str, Gene | GeneAA]:
    """Create a lookup dictionary for genes using species|locus|name as key."""
    return {gene_lookup_key(gene.species, gene.locus, gene.name): gene for gene in genes}
