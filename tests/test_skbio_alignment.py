from ssw_aligner import StripedSmithWaterman

from riot_na.alignment.gumbel import get_gumbel_params
from riot_na.alignment.params import NTSSWParams
from riot_na.alignment.skbio_alignment import align


def test_align():
    # given
    query = "CTATACTACTATGGTTCGGGGAGTTATTATAGCCTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCTCGT"
    target = "ACTACTTTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"

    # when
    aligner = StripedSmithWaterman(query)
    alignment = align(aligner, "IGHJ4*02", target, 100, query, get_gumbel_params(NTSSWParams()))

    # then

    assert query[alignment.q_start : alignment.q_end] == "TTGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"


if __name__ == "__main__":
    test_align()
