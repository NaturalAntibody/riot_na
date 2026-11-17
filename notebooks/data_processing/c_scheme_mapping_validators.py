"""Validation of C germline numbering and scheme mappings built by c_scheme_mappings.py."""

from typing import Optional, Sequence

from riot_na.data.model import AlignmentString, ChainType, CScheme, SchemeAlignment
from riot_na.schemes.scheme_alignment import produce_c_numbering

CONSERVED_CYSTEINE_POSITIONS = ("23", "104")


def validate_numbered_germline(
    gene: str, germline: str, imgt_numbering: Sequence[tuple[str, str]], imgt_table: dict[str, str]
) -> None:
    """
    imgt_numbering - numbered germline residues with their IMGT positions,
    imgt_table - IMGT positions and residues of the gene, if IMGT gives them.

    Numbered residues start the germline (only membrane exons are left), every residue has its own IMGT position,
    residues equal the ones in the IMGT table and C domains have the conserved cysteines.
    """
    residues = "".join(residue for residue, _ in imgt_numbering)
    assert germline.startswith(residues), f"{gene} numbered residues differ from the germline"

    imgt_residues = dict((position, residue) for residue, position in imgt_numbering)
    assert len(imgt_residues) == len(imgt_numbering), f"{gene} has residues on the same IMGT position"

    for position, residue in imgt_residues.items():
        assert imgt_table.get(position, residue) == residue, f"{gene} {position} differs from IMGT: {residue}"

    domains = {position.split("_")[0] for position in imgt_residues if position.startswith(("CH", "CL"))}
    for domain in domains:
        for position in CONSERVED_CYSTEINE_POSITIONS:
            assert imgt_residues.get(f"{domain}_{position}") == "C", f"{gene} has no cysteine on {domain}_{position}"


def validate_scheme_cigar(
    gene: str,
    residues: str,
    positions: Sequence[Optional[str]],
    cigar: str,
    c_scheme: CScheme,
    chain_type: ChainType,
) -> None:
    """
    Numbering germline residues with the scheme cigar, the same way riot numbers queries, gives their positions.
    Residues without position (insertions without EU position) are only checked to be numbered.
    """
    assert cigar.count("M") + cigar.count("I") == len(residues), f"{gene} {c_scheme.value} cigar length"

    scheme_alignment = SchemeAlignment(
        target_id=c_scheme.value, q_start=0, q_end=len(residues), alignment_str=AlignmentString(cigar)
    )
    numbering = produce_c_numbering(residues, scheme_alignment, c_scheme, chain_type)
    numbered_positions = list(numbering)

    assert "".join(numbering.values()) == residues, f"{gene} {c_scheme.value} numbering skips residues"
    for i, (expected, numbered) in enumerate(zip(positions, numbered_positions)):
        assert expected in (None, numbered), f"{gene} {c_scheme.value} residue {i}: {numbered} instead of {expected}"
