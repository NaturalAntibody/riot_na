"""
Utility functions for handling gene matches and species information.

This module provides shared functionality for converting between different gene match formats
and extracting species information, avoiding the need for proxy classes and code duplication.
"""

from typing import Sequence

from riot_na.data.model import Gene, GeneAA, Locus, Organism, SpeciesGeneMatch
from riot_na.riot_na import GeneMatch  # pylint: disable=no-name-in-module


def create_gene_lookup(genes: Sequence[Gene | GeneAA]) -> dict[str, Gene | GeneAA]:
    """Create a lookup dictionary for genes using species|name as key."""
    return {f"{gene.species.value}|{gene.locus.value}|{gene.name}": gene for gene in genes}


def parse_species_gene_match_simple(gene_match: GeneMatch) -> SpeciesGeneMatch:
    """
    Parse a GeneMatch to extract species information from a simple gene_id format.

    This function handles the format used by SegmentGeneAligner where gene_ids
    can be in formats like "species|name" or "species|locus|name".

    Args:
        gene_match: The raw gene match from prefiltering

    Returns:
        SpeciesGeneMatch with parsed species and locus information
    """
    parts = gene_match.gene_id.split("|")

    if len(parts) >= 3:
        # Format: "species|locus|name"
        species_str, locus_str, name = parts[0], parts[1], parts[2]
        locus = Locus(locus_str)
    elif len(parts) == 2:
        # Format: "species|name"
        species_str, name = parts[0], parts[1]
        locus = Locus.IGH  # Default locus
    else:
        # Single part - assume it's just the name
        species_str = "HOMO_SAPIENS"  # Default species
        name = gene_match.gene_id
        locus = Locus.IGH  # Default locus

    return SpeciesGeneMatch(
        species_gene_id=gene_match.gene_id,
        gene_id=name,
        rev_comp=gene_match.rev_comp,
        coverage=gene_match.coverage,
        species=Organism(species_str),
        locus=locus,
    )
