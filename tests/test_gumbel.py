# pylint: disable=unused-argument,protected-access
import json

import pytest

from riot_na.alignment import gumbel
from riot_na.alignment.aa_gene_alignments import aa_genes_dir
from riot_na.alignment.alignment_metrics import GumbelParams
from riot_na.alignment.gene_aligner import create_aa_c_gene_aligner
from riot_na.alignment.params import DEFAULT_ALIGNER_PARAMS, AASSWParams, NTSSWParams
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import Locus, Organism


@pytest.fixture(name="cache_path")
def fixture_cache_path(tmp_path, monkeypatch):
    path = tmp_path / "cache.json"
    monkeypatch.setattr(gumbel, "CACHE_PATH", path)
    gumbel.get_gumbel_params.cache_clear()
    yield path
    gumbel.get_gumbel_params.cache_clear()


def test_defaults_are_precomputed(monkeypatch, cache_path):
    def fail(_):
        raise AssertionError("ALP must not run for the default settings")

    monkeypatch.setattr(gumbel, "compute_gumbel_params", fail)
    for ssw_params in DEFAULT_ALIGNER_PARAMS.ssw_params():
        params = gumbel.get_gumbel_params(ssw_params)
        assert params.lambda_ > 0 and params.K > 0
    assert not cache_path.exists()


def test_overridden_settings_are_computed_once_and_cached(monkeypatch, cache_path):
    calls = []

    def fake_compute(ssw_params):
        calls.append(ssw_params)
        return GumbelParams(lambda_=0.3, K=0.1)

    monkeypatch.setattr(gumbel, "compute_gumbel_params", fake_compute)
    custom = AASSWParams(gap_open_penalty=10, gap_extend_penalty=2)

    assert gumbel.get_gumbel_params(custom) == GumbelParams(lambda_=0.3, K=0.1)
    assert json.loads(cache_path.read_text())["params"]["aa|matrix=blosum62|gap_open=10|gap_extend=2"]["K"] == 0.1

    gumbel.get_gumbel_params.cache_clear()  # new process: loaded from the cache file
    assert gumbel.get_gumbel_params(custom) == GumbelParams(lambda_=0.3, K=0.1)
    assert len(calls) == 1


def test_cache_with_other_alp_settings_is_ignored(monkeypatch, cache_path):
    custom = NTSSWParams(match_score=2, mismatch_score=-3)
    key = gumbel.settings_key(custom)
    cache_path.write_text(json.dumps({"alp_settings": {"n_samples": 1}, "params": {key: {"lambda": 9.0, "K": 9.0}}}))
    monkeypatch.setattr(gumbel, "compute_gumbel_params", lambda _: GumbelParams(lambda_=0.5, K=0.2))
    assert gumbel.get_gumbel_params(custom) == GumbelParams(lambda_=0.5, K=0.2)


def test_alp_reproduces_precomputed_nt_defaults(cache_path):
    pytest.importorskip("alp")
    ssw_params = DEFAULT_ALIGNER_PARAMS.v_nt.aligner
    precomputed = gumbel._read(gumbel.PRECOMPUTED_PATH)[gumbel.settings_key(ssw_params)]
    computed = gumbel.compute_gumbel_params(ssw_params)
    assert computed.lambda_ == pytest.approx(precomputed.lambda_, rel=0.02)
    assert computed.K == pytest.approx(precomputed.K, rel=0.05)


def test_aligner_uses_gumbel_params_of_its_ssw_settings(cache_path, monkeypatch):
    monkeypatch.setattr(gumbel, "compute_gumbel_params", lambda _: GumbelParams(lambda_=0.31, K=0.07))
    custom = DEFAULT_ALIGNER_PARAMS.with_param("c_aa.aligner.gap_extend_penalty", 3)
    aligner = create_aa_c_gene_aligner(Organism.HOMO_SAPIENS, Locus.IGH, aa_genes_dir(GENE_DB_DIR), custom)
    assert aligner.gumbel_params == GumbelParams(lambda_=0.31, K=0.07)
    default = create_aa_c_gene_aligner(Organism.HOMO_SAPIENS, Locus.IGH, aa_genes_dir(GENE_DB_DIR))
    assert default.gumbel_params == gumbel.get_gumbel_params(DEFAULT_ALIGNER_PARAMS.c_aa.aligner)
