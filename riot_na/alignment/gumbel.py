"""
Gumbel parameters (lambda, K) of local alignment scores, used to convert raw Smith-Waterman scores to bit scores and
E-values. They depend on the scoring system (substitution matrix or match / mismatch scores and gap penalties) and are
computed by ALP (https://github.com/NaturalAntibody/alp), which takes from a fraction of a second to about a minute.

Lookup order of ``get_gumbel_params``:

1. memory,
2. ``riot_na/data/gumbel_params.json``: precomputed for the default settings (regenerate with
   ``python -m riot_na.alignment.gumbel``),
3. ``riot_na/data/gumbel_params_cache.json``: parameters computed earlier for settings overridden by the user,
4. computed with ALP and saved to the cache file (kept in memory only if the data directory is not writable).
"""

from __future__ import annotations

import json
import logging
import os
import tempfile
from functools import lru_cache
from pathlib import Path
from typing import Optional, Union

from riot_na.alignment.alignment_metrics import GumbelParams
from riot_na.alignment.params import DEFAULT_ALIGNER_PARAMS, AASSWParams, NTSSWParams

logger = logging.getLogger(__name__)

DATA_DIR = Path(__file__).parent.parent / "data"
PRECOMPUTED_PATH = DATA_DIR / "gumbel_params.json"
CACHE_PATH = DATA_DIR / "gumbel_params_cache.json"

# ALP settings; part of the files, parameters computed with other settings are not reused.
ALP_SETTINGS = {"n_samples": 100_000, "seed": 42, "gap_convention": "open_includes_extend"}

SSWParams = Union[NTSSWParams, AASSWParams]


def settings_key(ssw_params: SSWParams) -> str:
    """Key of the scoring system, e.g. ``nt|match=1|mismatch=-1|gap_open=4|gap_extend=1``."""
    if isinstance(ssw_params, NTSSWParams):
        return (
            f"nt|match={ssw_params.match_score}|mismatch={ssw_params.mismatch_score}"
            f"|gap_open={ssw_params.gap_open_penalty}|gap_extend={ssw_params.gap_extend_penalty}"
        )
    return f"aa|matrix={ssw_params.matrix}|gap_open={ssw_params.gap_open_penalty}|gap_extend={ssw_params.gap_extend_penalty}"


def _read(path: Path) -> dict[str, GumbelParams]:
    try:
        content = json.loads(path.read_text())
    except (OSError, ValueError):
        return {}
    if content.get("alp_settings") != ALP_SETTINGS:
        return {}
    return {
        key: GumbelParams(lambda_=value["lambda"], K=value["K"]) for key, value in content.get("params", {}).items()
    }


def _write(path: Path, params: dict[str, GumbelParams]) -> None:
    """Atomic write (several processes can compute parameters at once)."""
    content = {
        "alp_settings": ALP_SETTINGS,
        "params": {key: {"lambda": p.lambda_, "K": p.K} for key, p in sorted(params.items())},
    }
    fd, tmp_name = tempfile.mkstemp(dir=path.parent, prefix=path.name, suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as tmp:
            json.dump(content, tmp, indent=2)
            tmp.write("\n")
        os.replace(tmp_name, path)
    except OSError:
        if os.path.exists(tmp_name):
            os.remove(tmp_name)
        raise


def compute_gumbel_params(ssw_params: SSWParams) -> GumbelParams:
    """Run ALP for the scoring system (no caching)."""
    try:
        import alp  # type: ignore[import-untyped]  # pylint: disable=import-outside-toplevel
    except ImportError as error:
        raise ImportError(
            f"Gumbel parameters of {settings_key(ssw_params)} are not precomputed and the `alp` package is needed to "
            "compute them: pip install git+https://github.com/NaturalAntibody/alp.git"
        ) from error

    logger.info("Computing Gumbel parameters of %s with ALP", settings_key(ssw_params))
    if isinstance(ssw_params, NTSSWParams):
        result = alp.compute_gumbel_params(
            match=ssw_params.match_score,
            mismatch=-ssw_params.mismatch_score,
            gap_open=ssw_params.gap_open_penalty,
            gap_extend=ssw_params.gap_extend_penalty,
            **ALP_SETTINGS,
        )
    else:
        result = alp.compute_gumbel_params(
            ssw_params.substitution_matrix(),
            gap_open=ssw_params.gap_open_penalty,
            gap_extend=ssw_params.gap_extend_penalty,
            **ALP_SETTINGS,
        )
    return GumbelParams(lambda_=float(result.lambda_), K=float(result.K))


@lru_cache(maxsize=None)
def get_gumbel_params(ssw_params: SSWParams) -> GumbelParams:
    """Gumbel parameters of the scoring system: precomputed, cached or (and then cached) computed."""
    key = settings_key(ssw_params)
    for path in (PRECOMPUTED_PATH, CACHE_PATH):
        found = _read(path).get(key)
        if found is not None:
            return found

    computed = compute_gumbel_params(ssw_params)
    try:
        _write(CACHE_PATH, {**_read(CACHE_PATH), key: computed})
    except OSError as error:
        logger.warning("Could not save Gumbel parameters to %s: %s", CACHE_PATH, error)
    return computed


def precompute_defaults(output: Optional[Path] = None) -> dict[str, GumbelParams]:
    """Compute (ignoring files) and write the Gumbel parameters of all default scoring systems."""
    params = {settings_key(ssw): compute_gumbel_params(ssw) for ssw in DEFAULT_ALIGNER_PARAMS.ssw_params()}
    _write(output or PRECOMPUTED_PATH, params)
    return params


if __name__ == "__main__":
    for name, value in precompute_defaults().items():
        print(f"{name}: lambda={value.lambda_:.6f} K={value.K:.6f}")
