"""Alpaca ighJ-5 lacks the last S of the other J genes (...TVS instead of ...TVSS), the extension adds it."""

import pytest

from riot_na import Organism, RiotNumberingAA, Scheme, create_riot_aa
from riot_na.data.model import AlignmentEntryAA, Cigar, Locus
from riot_na.schemes.scheme_alignment import force_c_terminus_matches

VHH = "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVAAISWSGGSTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA"
CDR3 = "GDRYPQFEY"  # ends with the last residues of ighJ-5 (PQFEY)


@pytest.fixture(name="riot_aa", scope="module")
def fixture_riot_aa() -> RiotNumberingAA:
    return create_riot_aa(allowed_species=(Organism.VICUGNA_PACOS,))


def number(riot_aa: RiotNumberingAA, j_region: str, extend_alignment: bool):
    return riot_aa.run_on_sequence("", VHH + CDR3 + j_region, scheme=Scheme.IMGT, extend_alignment=extend_alignment)


def test_missing_last_residue_is_numbered_with_extension(riot_aa: RiotNumberingAA):
    sequence = VHH + CDR3 + "WGQGTLVTVSS"
    airr = number(riot_aa, "WGQGTLVTVSS", extend_alignment=True)

    assert airr.j_call == "ighJ-5"
    assert airr.scheme_residue_mapping
    assert list(airr.scheme_residue_mapping.items())[-2:] == [("127", "S"), ("128", "S")]
    assert airr.fwr4_aa == "WGQGTLVTVSS"
    assert airr.sequence_end_aa_extended == len(sequence)
    # the alignment itself is not changed, the residue is only added to the extended alignment
    assert airr.sequence_alignment_aa and airr.sequence_alignment_aa.endswith("WGQGTLVTVS")


def test_missing_last_residue_is_not_added_without_extension(riot_aa: RiotNumberingAA):
    airr = number(riot_aa, "WGQGTLVTVSS", extend_alignment=False)

    assert airr.j_call == "ighJ-5"
    assert airr.scheme_residue_mapping
    assert list(airr.scheme_residue_mapping)[-1] == "127"
    assert airr.fwr4_aa == "WGQGTLVTVS"
    assert airr.sequence_alignment_aa_extended is None


@pytest.mark.parametrize("tail", ["", "A", "SS", "HHHHHH"])
def test_only_a_single_missing_residue_is_added(riot_aa: RiotNumberingAA, tail: str):
    airr = number(riot_aa, "WGQGTLVTVS" + tail, extend_alignment=True)

    assert airr.scheme_residue_mapping
    expected_last = "128" if tail.startswith("S") else "127"
    assert list(airr.scheme_residue_mapping)[-1] == expected_last
    assert airr.fwr4_aa == ("WGQGTLVTVSS" if tail.startswith("S") else "WGQGTLVTVS")


def test_other_alpaca_j_gene_is_unchanged(riot_aa: RiotNumberingAA):
    # ighJ-4 has the S on its own
    airr = riot_aa.run_on_sequence("", VHH + "GDRYEYDYWGQGTQVTVSS", scheme=Scheme.IMGT, extend_alignment=True)

    assert airr.j_call == "ighJ-4"
    assert airr.scheme_residue_mapping
    assert list(airr.scheme_residue_mapping.items())[-1] == ("128", "S")
    assert airr.fwr4_aa == "WGQGTQVTVSS"


def alignment(species: Organism, target_id: str, target: str, q_end: int, cigar: str, query: str):
    return AlignmentEntryAA(
        target_id=target_id,
        alignment_score=None,
        seq_identity=None,
        log10_e_value=None,
        q_start=0,
        q_end=q_end,
        t_start=0,
        t_end=q_end,
        cigar=Cigar(cigar),
        species=species,
        locus=Locus.IGH,
        q_seq=query,
        t_seq=target,
    )


@pytest.mark.parametrize(
    "species, target_id, expected_cigar, expected_q_end",
    [
        (Organism.VICUGNA_PACOS, "ighJ-5", "15M0M1I", 16),
        (Organism.VICUGNA_PACOS, "ighJ-4", "15M0M", 15),  # other gene
        (Organism.HOMO_SAPIENS, "ighJ-5", "15M0M", 15),  # other species
    ],
)
def test_force_c_terminus_matches_adds_the_residue_only_for_the_gene(
    species: Organism, target_id: str, expected_cigar: str, expected_q_end: int
):
    germline = "PQFEYWGQGTLVTVS"
    query = germline + "SHHH"
    result = force_c_terminus_matches(query, germline, alignment(species, target_id, germline, 15, "15M", query))

    assert str(result.cigar) == expected_cigar
    assert result.q_end == expected_q_end
    assert result.t_end == 15
