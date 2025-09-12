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
        "", SEQUENCES_NT["H"] + "GGGGGGSGGGG" + SEQUENCES_NT["L"], scheme=Scheme.IMGT
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


if __name__ == "__main__":
    test_e2e_nucleotides()
    test_e2e_amino_acids()
    print("All tests passed!")
