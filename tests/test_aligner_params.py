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

    assert (params.ch_nt.thresholds.e_value_threshold, params.ch_nt.thresholds.min_seq_identity) == (1.35e-05, 0.8)
    assert (params.c_aa.thresholds.e_value_threshold, params.c_aa.thresholds.min_seq_identity) == (4e-12, 0.8)
    assert params.v_nt.aligner.kwargs() == {
        "match_score": 1,
        "mismatch_score": -1,
        "gap_open_penalty": 4,
        "gap_extend_penalty": 1,
    }
    assert (params.v_aa.aligner.gap_open_penalty, params.v_aa.aligner.gap_extend_penalty) == (11, 1)
    assert (params.c_aa.aligner.gap_open_penalty, params.c_aa.aligner.gap_extend_penalty) == (11, 5)
    assert params.nt_gene(GermlineGene.J, Locus.IGL).prefiltering.distance_threshold == 7
    assert params.aa_gene(GermlineGene.J).thresholds.max_cdr3_length == 60


def test_params_are_immutable_and_hashable():
    with pytest.raises(FrozenInstanceError):
        DEFAULT_ALIGNER_PARAMS.v_nt = DEFAULT_ALIGNER_PARAMS.ch_nt  # type: ignore[misc]

    assert hash(DEFAULT_ALIGNER_PARAMS) == hash(replace(DEFAULT_ALIGNER_PARAMS))
    assert hash(DEFAULT_ALIGNER_PARAMS) != hash(DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4))


def test_with_c_min_seq_identity_changes_all_c_aligners_only():
    params = DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4)

    assert [p.thresholds.min_seq_identity for p in (params.ch_nt, params.ck_nt, params.cl_nt, params.c_aa)] == [0.4] * 4
    assert params.v_nt == DEFAULT_ALIGNER_PARAMS.v_nt
    assert params.c_aa.aligner == DEFAULT_ALIGNER_PARAMS.c_aa.aligner


def test_params_reach_the_aligners():
    params = DEFAULT_ALIGNER_PARAMS.with_param("c_aa.thresholds.min_seq_identity", 0.5).with_param(
        "c_aa.aligner.gap_extend_penalty", 3
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
        UserWarning, match=r"minimum sequence identity of C alignments is 80%.*pass it as the last argument"
    ):
        get_or_create_riot_nt(c_genes_organism=Organism.HOMO_SAPIENS, params=replace(DEFAULT_ALIGNER_PARAMS))


def test_no_warning_without_forced_c_organism(recwarn):
    warn_c_genes_organism(None, DEFAULT_ALIGNER_PARAMS)
    warn_c_genes_organism(Organism.HOMO_SAPIENS, DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(None))

    assert not recwarn.list


def test_no_warning_when_identity_is_lowered(recwarn):
    warn_c_genes_organism(Organism.HOMO_SAPIENS, DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4))
    get_or_create_riot_aa(
        c_genes_organism=Organism.HOMO_SAPIENS, params=DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.4)
    )

    assert not [w for w in recwarn.list if issubclass(w.category, UserWarning)]


def test_warning_if_identity_is_lowered_only_for_some_c_aligners():
    params = DEFAULT_ALIGNER_PARAMS.with_param("c_aa.thresholds.min_seq_identity", 0.4)

    with pytest.warns(UserWarning, match="minimum sequence identity of C alignments"):
        warn_c_genes_organism(Organism.HOMO_SAPIENS, params)

    # at or above the default is not lowered
    with pytest.warns(UserWarning):
        warn_c_genes_organism(Organism.HOMO_SAPIENS, DEFAULT_ALIGNER_PARAMS.with_c_min_seq_identity(0.9))


def test_with_param_overrides_a_nested_setting():
    params = DEFAULT_ALIGNER_PARAMS.with_param("v_aa.prefiltering.kmer_size", 4)

    assert params.v_aa.prefiltering.kmer_size == 4
    assert DEFAULT_ALIGNER_PARAMS.v_aa.prefiltering.kmer_size == 3
    assert params.v_aa.thresholds == DEFAULT_ALIGNER_PARAMS.v_aa.thresholds
    assert params.c_aa == DEFAULT_ALIGNER_PARAMS.c_aa
    assert params.with_param("v_aa.aligner.gap_extend_penalty", 2).v_aa.aligner.gap_extend_penalty == 2
    assert params.with_param("v_aa.aligner.gap_extend_penalty", 2).v_aa.prefiltering.kmer_size == 4


@pytest.mark.parametrize("path", ["v_aa.prefiltering.kmer", "nope", "v_aa.prefiltering.kmer_size.deeper"])
def test_with_param_rejects_unknown_paths(path):
    with pytest.raises(AttributeError, match="Invalid setting"):
        DEFAULT_ALIGNER_PARAMS.with_param(path, 1)
