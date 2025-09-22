"""
Utility functions for handling gene matches and species information.

This module provides shared functionality for converting between different gene match formats
and extracting species information, avoiding the need for proxy classes and code duplication.
"""

from typing import Sequence

from riot_na.data.model import Gene, GeneAA


def create_gene_lookup(genes: Sequence[Gene | GeneAA]) -> dict[str, Gene | GeneAA]:
    """Create a lookup dictionary for genes using species|name as key."""
    return {f"{gene.species.value}|{gene.locus.value}|{gene.name}": gene for gene in genes}
