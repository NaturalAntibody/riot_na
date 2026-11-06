from dataclasses import FrozenInstanceError, replace

import pytest

from riot_na import (
    DEFAULT_ALIGNER_PARAMS,
    Organism,
    get_or_create_riot_aa,
    get_or_create_riot_nt,
)
from riot_na.alignment.params import warn_c_genes_organism
from riot_na.data.model import GermlineGene, Locus

SEQUENCE_AA = (
    "QVQLQQWGAGLLKPSETLSLTCAVFGGSFSGYYWSWIRQPPGKGLEWIGEINHRGNTNDNPSLKSRVTISVDTSKNQFALKLSSVTAADTAVYYCARERGYTYGNFDH"
    "WGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKK"
    "VEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCK"
)


def test_defaults_are_the_established_settings():
    params = DEFAULT_ALIGNER_PARAMS

    assert (params.ch_nt.e_value_threshold, params.ch_nt.min_seq_identity) == (1e-05, 0.85)
    assert (params.c_aa.e_value_threshold, params.c_aa.min_seq_identity) == (1e-10, 0.85)
    assert params.nt_ssw.kwargs() == {
        "match_score": 1,
        "mismatch_score": -1,
        "gap_open_penalty": 4,
        "gap_extend_penalty": 1,
    }
    assert (params.aa_ssw.gap_open_penalty, params.aa_ssw.gap_extend_penalty) == (11, 1)
    assert (params.c_aa_ssw.gap_open_penalty, params.c_aa_ssw.gap_extend_penalty) == (11, 5)
    assert params.nt_gene(GermlineGene.J, Locus.IGL).distance_threshold == 7
    assert params.aa_gene(GermlineGene.J).max_cdr3_length == 60


def test_params_are_immutable_and_hashable():
    with pytest.raises(FrozenInstanceError):
        DEFAULT_ALIGNER_PARAMS.v_nt = DEFAULT_ALIGNER_PARAMS.ch_nt  # type: ignore[misc]

    assert hash(DEFAULT_ALIGNER_PARAMS) == hash(replace(DEFAULT_ALIGNER_PARAMS))
    assert hash(DEFAULT_ALIGNER_PARAMS) != hash(DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4))


def test_with_c_min_seq_identity_changes_all_c_aligners_only():
    params = DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4)

    assert [p.min_seq_identity for p in (params.ch_nt, params.ck_nt, params.cl_nt, params.c_aa)] == [0.4] * 4
    assert params.v_nt == DEFAULT_ALIGNER_PARAMS.v_nt
    assert params.c_aa_ssw == DEFAULT_ALIGNER_PARAMS.c_aa_ssw


def test_params_reach_the_aligners():
    params = replace(
        DEFAULT_ALIGNER_PARAMS,
        c_aa=replace(DEFAULT_ALIGNER_PARAMS.c_aa, min_seq_identity=0.5),
        c_aa_ssw=replace(DEFAULT_ALIGNER_PARAMS.c_aa_ssw, gap_extend_penalty=3),
    )

    riot = get_or_create_riot_aa(params=params)

    c_aligner = riot.vjc_aligner_aa.c_aligners[Organism.HOMO_SAPIENS][Locus.IGH]  # type: ignore[index]
    assert c_aligner.min_seq_identity == 0.5
    assert c_aligner.ssw_params.gap_extend_penalty == 3
    assert get_or_create_riot_aa(params=params) is riot
    assert get_or_create_riot_aa() is not riot


def test_c_identity_threshold_decides_if_diverged_c_region_is_aligned():
    # every 4th residue of the constant region replaced, about 75% identity to the germline
    start = SEQUENCE_AA.index("ASTKGPS")
    constant = "".join(
        "A" if i % 4 == 0 and residue != "A" else residue for i, residue in enumerate(SEQUENCE_AA[start:])
    )
    sequence = SEQUENCE_AA[:start] + constant

    default = get_or_create_riot_aa().run_on_sequence("a", sequence)
    loose = get_or_create_riot_aa(params=DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4)).run_on_sequence(
        "a", sequence
    )

    assert default.c_call is None
    assert loose.c_call == "IGHG1"


def test_forcing_c_organism_warns_about_identity_threshold():
    with pytest.warns(
        UserWarning, match=r"minimum sequence identity of C alignments is 85%.*pass it as the last argument"
    ):
        get_or_create_riot_nt(c_genes_organism=Organism.HOMO_SAPIENS, params=replace(DEFAULT_ALIGNER_PARAMS))


def test_no_warning_without_forced_c_organism(recwarn):
    warn_c_genes_organism(None, DEFAULT_ALIGNER_PARAMS)
    warn_c_genes_organism(Organism.HOMO_SAPIENS, DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(None))

    assert not recwarn.list
