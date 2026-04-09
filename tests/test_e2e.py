from pathlib import Path
from tempfile import TemporaryDirectory

from riot_na import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    ChainType,
    Locus,
    Scheme,
)
from riot_na.api.riot_numbering import get_or_create_riot_aa, get_or_create_riot_nt
from riot_na.common.io import write_airr_iter_to_csv
from riot_na.data.model import SegmentedAirrRearrangementEntryNT

SEQUENCES_NT = {
    "H": "AACAACACATGTCCAATGTCCTCTCCACAGACACTGAACACACTGACTCTAACCATGGGAAGGAGCTGGATCTTTCTCTTCCTCCTGTCAGGAACTGCAGGTGTCCACTCTGAGGTCCAGCTGCAACAGTCTGGACCTGTGCTGGTGAAGCCTGGGGCTTCAGTGAAGATGTCCTGTAAGGCTTCTGGATACACATTCACTGACTACTATATGAACTGGGTGAAGCAGAGCCATGGAAAGAGACTTGAGTGGATTGGAGTTATTAATCCTTACAACGGTGGTACTAACTATAACCAGAAGTTCAAGGGCAAGGCCACATTGACTGTTGACAAGTCCTCCAGCACAGCCTACATGGAGCTCAACAGCCTGACATCTGAGGACTCTGCAGTCTATTACTGTGCAGATGGGATTATTACGAATTGGTATTTCGATGTCTGGGGCACAGGGACCACGGTCACCGTCTCCTCAGCCAAAACGACACCCCCATCTGTCTATCCACTGGCCCCTGGATCTGCTGCCCAAACTAACTCCATGGTGACCCTGGGATGCCTGGTCAAGGGCTATTTCCCTGAGCCAGTGACAGTGACCTGGAACTCTGGATCCCTGTCCAGCGGTGTGCACACCTTCCCAGCTGTCCTGCAGTCTGACCTCTACACTCTGAGCAGCTCAGTGACTGTCCCCTCCAGCACCTGGCCCAGCCAGACCGTCACCTGCAACGTTGCCCACCCGGCCAGCAGCACCAAGGTGGACAAGAAAATTGTGCCCAGGGATTGTGGTTGTAAGCCTTGCATATGTACAGTCCCAGAAGT",
    "L": "GCTGACCAATATTGAAAAGAATAGACCTGGTTTGTGAATTATGGCCTGGATTTCACTTATACTCTCTCTCCTGGCTCTCAGCTCAGGGGCCATTTCCCAGGCTGTTGTGACTCAGGAATCTGCACTCACCACATCACCTGGTGAAACAGTCACACTCACTTGTCGCTCAAGTACTGGGGCTGTTACAACTAGTAACTATGCCAACTGGGTCCAAGAAAAACCAGATCATTTATTCACTGGTCTAATAGGTGGTACCAACAACCGAGCTCCAGGTGTTCCTGCCAGATTCTCAGGCTCCCTGATTGGAGACAAGGCTGCCCTCACCATCACAGGGGCACAGACTGAGGATGAGGCAATATATTTCTGTGCTCTATGGTACAGCAACCATTTCCACAATGACATGTGTAGATGGGGAAGTAGAACAAGAACACTCTGGTACAGTCTCATAACT",
}


def test_e2e_nucleotides():
    riot_numbering_nt = get_or_create_riot_nt()

    for scheme in Scheme:
        airr_heavy: AirrRearrangementEntryNT = riot_numbering_nt.run_on_sequence("", SEQUENCES_NT["H"], scheme=scheme)
        assert airr_heavy.v_call
        assert airr_heavy.j_call
        assert ChainType.from_locus(Locus(airr_heavy.locus)) == ChainType.HEAVY

        airr_light: AirrRearrangementEntryNT = riot_numbering_nt.run_on_sequence("", SEQUENCES_NT["L"], scheme=scheme)
        assert airr_light.v_call
        assert airr_light.j_call
        assert ChainType.from_locus(Locus(airr_light.locus)) == ChainType.LIGHT

    riot_numbering_nt = get_or_create_riot_nt(return_all_domains=True)
    airrs_list: list[SegmentedAirrRearrangementEntryNT] = riot_numbering_nt.run_on_sequence(
        "", SEQUENCES_NT["H"] + "GGGGGGSGGGG" + SEQUENCES_NT["L"], scheme=Scheme.IMGT, return_all_domains=True
    )
    assert airrs_list[0].v_call
    assert airrs_list[0].j_call
    assert ChainType.from_locus(Locus(airrs_list[0].locus)) == ChainType.HEAVY
    assert airrs_list[1].v_call
    assert airrs_list[1].j_call
    assert ChainType.from_locus(Locus(airrs_list[1].locus)) == ChainType.LIGHT

    with TemporaryDirectory() as temp_dir:
        write_airr_iter_to_csv(Path(temp_dir) / "output.csv", SegmentedAirrRearrangementEntryNT, airrs_list)


SEQUENCES_AA = {
    "H": "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHHHHHH",
    "L": "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC",
}


def test_e2e_amino_acids():
    riot_numbering_aa = get_or_create_riot_aa()

    for scheme in Scheme:
        airr_heavy: AirrRearrangementEntryAA = riot_numbering_aa.run_on_sequence("", SEQUENCES_AA["H"], scheme=scheme)
        assert airr_heavy.v_call
        assert airr_heavy.j_call
        assert ChainType.from_locus(Locus(airr_heavy.locus)) == ChainType.HEAVY

        airr_light: AirrRearrangementEntryAA = riot_numbering_aa.run_on_sequence("", SEQUENCES_AA["L"], scheme=scheme)
        assert airr_light.v_call
        assert airr_light.j_call
        assert ChainType.from_locus(Locus(airr_light.locus)) == ChainType.LIGHT


def test_e2e_amino_acids_extend_alignment_fields():
    riot_numbering_aa = get_or_create_riot_aa()

    # extend_alignment=False: all three fields must be None
    airr_no_ext: AirrRearrangementEntryAA = riot_numbering_aa.run_on_sequence(
        "", SEQUENCES_AA["H"], scheme=Scheme.IMGT, extend_alignment=False
    )
    assert airr_no_ext.sequence_start_aa_extended is None
    assert airr_no_ext.sequence_end_aa_extended is None
    assert airr_no_ext.sequence_alignment_aa_extended is None

    # extend_alignment=True (default): fields must be populated
    for seq_key in ("H", "L"):
        airr: AirrRearrangementEntryAA = riot_numbering_aa.run_on_sequence(
            "", SEQUENCES_AA[seq_key], scheme=Scheme.IMGT, extend_alignment=True
        )
        assert airr.sequence_start_aa_extended is not None
        assert airr.sequence_end_aa_extended is not None
        assert airr.sequence_alignment_aa_extended is not None

        # slice must be consistent with start/end
        assert (
            airr.sequence_alignment_aa_extended
            == airr.sequence_aa[airr.sequence_start_aa_extended - 1 : airr.sequence_end_aa_extended]
        )

        # extended start must be <= SW V start
        assert airr.sequence_start_aa_extended <= airr.v_sequence_start_aa

        # extended end must be >= SW J end (or V end when no J)
        j_or_v_end = airr.j_sequence_end_aa if airr.j_sequence_end_aa is not None else airr.v_sequence_end_aa
        assert airr.sequence_end_aa_extended >= j_or_v_end

        # positional_scheme_mapping first key must equal sequence_start_aa_extended - 1 (0-based)
        if airr.positional_scheme_mapping:
            first_key = next(iter(airr.positional_scheme_mapping))
            assert first_key == airr.sequence_start_aa_extended - 1


if __name__ == "__main__":
    test_e2e_nucleotides()
    test_e2e_amino_acids()
    print("All tests passed!")
