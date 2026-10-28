"""NGS AIRR baseline I/O and numbering helpers for hand-triggered regression.

Binary format: a single pickle blob with metadata + AIRR rearrangement entries.

Generate baselines (current ssw-aligner / riot stack)::

    poetry run python -m tests.regression.ngs_baseline generate
    poetry run python -m tests.regression.ngs_baseline generate --limit 100 --ncpu 8
    poetry run python -m tests.regression.ngs_baseline generate --suffix ssw  # -> baseline_airr_nt_ssw.pkl

Compare against baselines after swapping the aligner::

    poetry run python -m tests.regression.ngs_baseline compare --ncpu 8
    poetry run python -m tests.regression.ngs_baseline compare --suffix ssw

Compare two stored baseline sets by suffix, without re-numbering ('' = unsuffixed files)::

    poetry run python -m tests.regression.ngs_baseline diff ssw parasail
    poetry run python -m tests.regression.ngs_baseline diff '' parasail --no-aa

When mismatches are found, ``compare`` and ``diff`` write results to JSONL (default ``ngs_regression_result.jsonl``
next to the baselines, override with ``--output-path``): one record per differing entry
(with ``diffs``). Mismatch counts are printed to stdout.
"""

from __future__ import annotations

import functools
import itertools
import multiprocessing as mp
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import click
import tqdm
from psutil import cpu_count

from riot_na import CScheme
from riot_na.api.api_mp import _worker_call, _worker_initializer
from riot_na.common.fasta import parse_fasta
from riot_na.common.io import count_fasta_records
from riot_na.config import GENE_DB_DIR
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    AirrRearrangementEntryNT,
    InputType,
    Scheme,
)

from .common import (
    compare_entries,
    finish,
    limit_option,
    load_pickle,
    ncpu_option,
    output_path_option,
    report,
    save_pickle,
    with_suffix,
)

NGS_DIR = Path(__file__).resolve().parents[2] / "data" / "ngs_regression_data"
FASTA_NT = NGS_DIR / "ngs_sample_clean_nt.fasta"
FASTA_AA = NGS_DIR / "ngs_sample_clean_aa.fasta"
BASELINE_NT = NGS_DIR / "baseline_airr_nt.pkl"
BASELINE_AA = NGS_DIR / "baseline_airr_aa.pkl"
RESULT_JSONL = NGS_DIR / "ngs_regression_result.jsonl"

BASELINE_FORMAT_VERSION = 1
SCHEME = Scheme.IMGT


def number_fasta(
    fasta_path: Path,
    input_type: InputType,
    *,
    scheme: Scheme = SCHEME,
    limit: int | None = None,
    n_processes: int | None = None,
    db_dir: Path = GENE_DB_DIR,
) -> list[AirrRearrangementEntryNT] | list[AirrRearrangementEntryAA]:
    """Run riot NT/AA numbering over a FASTA file and return AIRR entries."""
    if n_processes is None:
        n_processes = cpu_count(logical=False) or 1

    total = count_fasta_records(fasta_path)
    if limit is not None:
        total = min(total, limit)

    records = itertools.islice(parse_fasta(fasta_path), limit)
    initializer = functools.partial(
        _worker_initializer,
        allowed_species=None,
        scheme=scheme,
        db_dir=db_dir,
        input_type=input_type,
        return_all_domains=False,
        extend_alignment=False,
        c_scheme=CScheme.IMGT,
    )
    with mp.Pool(processes=n_processes, initializer=initializer) as pool:
        results = list(
            tqdm.tqdm(
                pool.imap(_worker_call, records, chunksize=32),
                total=total,
                desc=f"numbering {fasta_path.name}",
            )
        )
    return results  # type: ignore[return-value]


def save_baseline(
    path: Path,
    entries: Sequence[AirrRearrangementEntryNT] | Sequence[AirrRearrangementEntryAA],
    *,
    fasta_path: Path,
    input_type: InputType,
    scheme: Scheme = SCHEME,
) -> None:
    """Serialize AIRR rearrangement entries to a pickle baseline file."""
    payload = {
        "format_version": BASELINE_FORMAT_VERSION,
        "fasta_name": fasta_path.name,
        "input_type": input_type.value,
        "scheme": scheme.value,
        "n_entries": len(entries),
        "entries": list(entries),
    }
    save_pickle(path, payload)


def load_baseline(path: Path) -> dict[str, Any]:
    """Load a baseline pickle produced by :func:`save_baseline`."""
    return load_pickle(path, BASELINE_FORMAT_VERSION)


def generate_baseline_nt(
    *,
    fasta_path: Path = FASTA_NT,
    baseline_path: Path = BASELINE_NT,
    limit: int | None = None,
    n_processes: int | None = None,
) -> Path:
    entries = number_fasta(fasta_path, InputType.NT, limit=limit, n_processes=n_processes)
    save_baseline(baseline_path, entries, fasta_path=fasta_path, input_type=InputType.NT)
    return baseline_path


def generate_baseline_aa(
    *,
    fasta_path: Path = FASTA_AA,
    baseline_path: Path = BASELINE_AA,
    limit: int | None = None,
    n_processes: int | None = None,
) -> Path:
    entries = number_fasta(fasta_path, InputType.AA, limit=limit, n_processes=n_processes)
    save_baseline(baseline_path, entries, fasta_path=fasta_path, input_type=InputType.AA)
    return baseline_path


def compare_baselines(baseline_a: Path, baseline_b: Path) -> dict[str, Any]:
    """Compare two stored baselines (a = expected, b = actual); see :func:`compare_entries`."""
    for path in (baseline_a, baseline_b):
        if not path.exists():
            raise FileNotFoundError(path)

    payload_a = load_baseline(baseline_a)
    payload_b = load_baseline(baseline_b)
    for key in ("input_type", "fasta_name", "scheme"):
        if payload_a[key] != payload_b[key]:
            raise ValueError(
                f"{key} differs: {baseline_a.name}={payload_a[key]!r} {baseline_b.name}={payload_b[key]!r}"
            )

    return compare_entries(payload_a["entries"], payload_b["entries"])


def compare_to_baseline(
    fasta_path: Path,
    baseline_path: Path,
    input_type: InputType,
    *,
    n_processes: int | None = None,
) -> dict[str, Any]:
    """Re-number FASTA and compare against baseline; see :func:`compare_entries`."""
    if not fasta_path.exists():
        raise FileNotFoundError(fasta_path)

    baseline = load_baseline(baseline_path)
    if baseline["input_type"] != input_type.value:
        raise ValueError(f"baseline input_type={baseline['input_type']!r} != {input_type.value!r}")
    if baseline["fasta_name"] != fasta_path.name:
        raise ValueError(f"baseline fasta_name={baseline['fasta_name']!r} != {fasta_path.name!r}")

    expected = baseline["entries"]
    actual = number_fasta(fasta_path, input_type, limit=len(expected), n_processes=n_processes)
    return compare_entries(expected, actual)


@click.group()
def cli() -> None:
    """NGS AIRR regression baseline tools."""


@cli.command("generate")
@limit_option
@ncpu_option
@click.option(
    "--suffix",
    type=str,
    default=None,
    help="Append '_<suffix>' to baseline file names, e.g. baseline_airr_nt_<suffix>.pkl.",
)
def generate_cmd(limit: int | None, ncpu: int | None, suffix: str = "") -> None:
    """Generate NT and AA baselines."""
    baseline_nt = with_suffix(BASELINE_NT, suffix)
    baseline_aa = with_suffix(BASELINE_AA, suffix)
    click.echo(f"Writing NT baseline -> {baseline_nt}")
    generate_baseline_nt(baseline_path=baseline_nt, limit=limit, n_processes=ncpu)
    click.echo(f"Writing AA baseline -> {baseline_aa}")
    generate_baseline_aa(baseline_path=baseline_aa, limit=limit, n_processes=ncpu)
    click.echo("Done.")


@cli.command("generate-nt")
@limit_option
@ncpu_option
@click.option(
    "--fasta",
    type=click.Path(path_type=Path, exists=True),
    default=FASTA_NT,
    show_default=True,
)
@click.option("--output", type=click.Path(path_type=Path), default=BASELINE_NT, show_default=True)
def generate_nt_cmd(limit: int | None, ncpu: int | None, fasta: Path, output: Path) -> None:
    """Generate NT baseline only."""
    click.echo(f"Writing NT baseline -> {output}")
    generate_baseline_nt(fasta_path=fasta, baseline_path=output, limit=limit, n_processes=ncpu)
    click.echo("Done.")


@cli.command("generate-aa")
@limit_option
@ncpu_option
@click.option(
    "--fasta",
    type=click.Path(path_type=Path, exists=True),
    default=FASTA_AA,
    show_default=True,
)
@click.option("--output", type=click.Path(path_type=Path), default=BASELINE_AA, show_default=True)
def generate_aa_cmd(limit: int | None, ncpu: int | None, fasta: Path, output: Path) -> None:
    """Generate AA baseline only."""
    click.echo(f"Writing AA baseline -> {output}")
    generate_baseline_aa(fasta_path=fasta, baseline_path=output, limit=limit, n_processes=ncpu)
    click.echo("Done.")


@cli.command("compare")
@ncpu_option
@click.option("--nt/--no-nt", default=True, show_default=True, help="Compare NT baseline.")
@click.option("--aa/--no-aa", default=True, show_default=True, help="Compare AA baseline.")
@click.option(
    "--suffix",
    type=str,
    default=None,
    help="Compare against baselines generated with the same '--suffix'.",
)
@output_path_option(RESULT_JSONL)
def compare_cmd(ncpu: int | None, nt: bool, aa: bool, output_path: Path, suffix: str = "") -> None:
    """Re-run numbering and compare against pickled baselines."""
    results = []
    for enabled, label, fasta_path, base_path, input_type in (
        (nt, "NT", FASTA_NT, BASELINE_NT, InputType.NT),
        (aa, "AA", FASTA_AA, BASELINE_AA, InputType.AA),
    ):
        if not enabled:
            continue
        baseline_path = with_suffix(base_path, suffix)
        click.echo(f"Comparing {label}: {fasta_path} vs {baseline_path}")
        result = compare_to_baseline(fasta_path, baseline_path, input_type, n_processes=ncpu)
        results.append({"dataset": label, **result})
        report(label, result)
    finish(results, output_path)


@cli.command("diff")
@click.argument("suffix_a")
@click.argument("suffix_b")
@click.option("--nt/--no-nt", default=True, show_default=True, help="Diff NT baselines.")
@click.option("--aa/--no-aa", default=True, show_default=True, help="Diff AA baselines.")
@output_path_option(RESULT_JSONL)
def diff_cmd(suffix_a: str, suffix_b: str, nt: bool, aa: bool, output_path: Path) -> None:
    """Compare two stored baseline sets identified by SUFFIX_A and SUFFIX_B ('' = unsuffixed).

    SUFFIX_A is reported as "expected", SUFFIX_B as "actual".
    """
    results = []
    for enabled, label, base_path in ((nt, "NT", BASELINE_NT), (aa, "AA", BASELINE_AA)):
        if not enabled:
            continue
        path_a = with_suffix(base_path, suffix_a)
        path_b = with_suffix(base_path, suffix_b)
        click.echo(f"Diffing {label}: {path_a} vs {path_b}")
        result = compare_baselines(path_a, path_b)
        results.append({"dataset": label, **result})
        report(label, result)
    finish(results, output_path)


if __name__ == "__main__":
    cli()
