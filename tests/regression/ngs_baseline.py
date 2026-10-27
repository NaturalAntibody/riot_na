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
"""

from __future__ import annotations

import itertools
import multiprocessing as mp
import pickle
import sys
from collections.abc import Sequence
from dataclasses import fields
from pathlib import Path
from typing import Any

import click
import tqdm
from psutil import cpu_count

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

NGS_DIR = Path(__file__).resolve().parents[2] / "data" / "ngs_regression_data"
FASTA_NT = NGS_DIR / "ngs_sample_clean_nt.fasta"
FASTA_AA = NGS_DIR / "ngs_sample_clean_aa.fasta"
BASELINE_NT = NGS_DIR / "baseline_airr_nt.pkl"
BASELINE_AA = NGS_DIR / "baseline_airr_aa.pkl"

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
    with mp.Pool(
        processes=n_processes,
        initializer=_worker_initializer,
        initargs=(None, scheme, db_dir, input_type, False, False),
    ) as pool:
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
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "format_version": BASELINE_FORMAT_VERSION,
        "fasta_name": fasta_path.name,
        "input_type": input_type.value,
        "scheme": scheme.value,
        "n_entries": len(entries),
        "entries": list(entries),
    }
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
    tmp_path.replace(path)


def load_baseline(path: Path) -> dict[str, Any]:
    """Load a baseline pickle produced by :func:`save_baseline`."""
    with open(path, "rb") as handle:
        payload = pickle.load(handle)
    if payload.get("format_version") != BASELINE_FORMAT_VERSION:
        raise ValueError(
            f"Unsupported baseline format_version={payload.get('format_version')!r} "
            f"(expected {BASELINE_FORMAT_VERSION})"
        )
    return payload


def _diff_entry(expected: Any, actual: Any) -> list[str]:
    diffs: list[str] = []
    for field in fields(expected):
        exp_val = getattr(expected, field.name)
        act_val = getattr(actual, field.name)
        if exp_val != act_val:
            diffs.append(f"{field.name}: expected={exp_val!r} actual={act_val!r}")
    return diffs


def compare_entries(
    expected: Sequence[Any],
    actual: Sequence[Any],
    *,
    max_report: int = 20,
) -> list[str]:
    """Return human-readable mismatch descriptions (empty if identical)."""
    problems: list[str] = []
    if len(expected) != len(actual):
        problems.append(f"length mismatch: baseline={len(expected)} actual={len(actual)}")

    for idx, (exp, act) in enumerate(zip(expected, actual)):
        if exp == act:
            continue
        header = getattr(exp, "sequence_header", None) or getattr(act, "sequence_header", idx)
        field_diffs = _diff_entry(exp, act) if type(exp) is type(act) else [f"type {type(exp)} vs {type(act)}"]
        problems.append(f"[{idx}] header={header!r}: " + "; ".join(field_diffs[:10]))
        if len(problems) >= max_report:
            problems.append(f"... truncated after {max_report} mismatches")
            break

    if len(expected) != len(actual):
        longer = "baseline" if len(expected) > len(actual) else "actual"
        problems.append(f"{longer} has {abs(len(expected) - len(actual))} extra entries")

    return problems


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


def compare_baselines(baseline_a: Path, baseline_b: Path) -> list[str]:
    """Compare two stored baselines and return mismatch descriptions (empty if identical)."""
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
) -> list[str]:
    """Re-number FASTA and return mismatch descriptions vs baseline (empty if OK)."""
    if not fasta_path.exists():
        raise FileNotFoundError(fasta_path)
    if not baseline_path.exists():
        raise FileNotFoundError(baseline_path)

    baseline = load_baseline(baseline_path)
    if baseline["input_type"] != input_type.value:
        raise ValueError(f"baseline input_type={baseline['input_type']!r} != {input_type.value!r}")
    if baseline["fasta_name"] != fasta_path.name:
        raise ValueError(f"baseline fasta_name={baseline['fasta_name']!r} != {fasta_path.name!r}")

    expected = baseline["entries"]
    actual = number_fasta(fasta_path, input_type, limit=len(expected), n_processes=n_processes)
    return compare_entries(expected, actual)


def _common_options(function):
    function = click.option(
        "--limit",
        type=int,
        default=None,
        show_default="all",
        help="Only process the first N sequences (debug).",
    )(function)
    function = click.option(
        "--ncpu",
        type=int,
        default=None,
        show_default="physical cores",
        help="Worker processes for numbering.",
    )(function)
    return function


@click.group()
def cli() -> None:
    """NGS AIRR regression baseline tools."""


def _with_suffix(path: Path, suffix: str) -> Path:
    """Return ``<stem>_<suffix><ext>`` for ``path`` (unchanged when suffix is empty)."""
    if not suffix:
        return path
    return path.with_name(f"{path.stem}_{suffix}{path.suffix}")


@cli.command("generate")
@_common_options
@click.option(
    "--suffix",
    type=str,
    default=None,
    help="Append '_<suffix>' to baseline file names, e.g. baseline_airr_nt_<suffix>.pkl.",
)
def generate_cmd(limit: int | None, ncpu: int | None, suffix: str = "") -> None:
    """Generate NT and AA baselines."""
    baseline_nt = _with_suffix(BASELINE_NT, suffix)
    baseline_aa = _with_suffix(BASELINE_AA, suffix)
    click.echo(f"Writing NT baseline -> {baseline_nt}")
    generate_baseline_nt(baseline_path=baseline_nt, limit=limit, n_processes=ncpu)
    click.echo(f"Writing AA baseline -> {baseline_aa}")
    generate_baseline_aa(baseline_path=baseline_aa, limit=limit, n_processes=ncpu)
    click.echo("Done.")


@cli.command("generate-nt")
@_common_options
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
@_common_options
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
@click.option(
    "--ncpu",
    type=int,
    default=None,
    show_default="physical cores",
    help="Worker processes for numbering.",
)
@click.option("--nt/--no-nt", default=True, show_default=True, help="Compare NT baseline.")
@click.option("--aa/--no-aa", default=True, show_default=True, help="Compare AA baseline.")
@click.option(
    "--suffix",
    type=str,
    default=None,
    help="Compare against baselines generated with the same '--suffix'.",
)
def compare_cmd(ncpu: int | None, nt: bool, aa: bool, suffix: str = "") -> None:
    """Re-run numbering and compare against pickled baselines."""
    baseline_nt = _with_suffix(BASELINE_NT, suffix)
    baseline_aa = _with_suffix(BASELINE_AA, suffix)
    failed = False
    if nt:
        click.echo(f"Comparing NT: {FASTA_NT} vs {baseline_nt}")
        problems = compare_to_baseline(FASTA_NT, baseline_nt, InputType.NT, n_processes=ncpu)
        if problems:
            failed = True
            click.echo("NT mismatches:", err=True)
            for problem in problems:
                click.echo(f"  {problem}", err=True)
        else:
            click.echo("NT OK")
    if aa:
        click.echo(f"Comparing AA: {FASTA_AA} vs {baseline_aa}")
        problems = compare_to_baseline(FASTA_AA, baseline_aa, InputType.AA, n_processes=ncpu)
        if problems:
            failed = True
            click.echo("AA mismatches:", err=True)
            for problem in problems:
                click.echo(f"  {problem}", err=True)
        else:
            click.echo("AA OK")
    if failed:
        sys.exit(1)


@cli.command("diff")
@click.argument("suffix_a")
@click.argument("suffix_b")
@click.option("--nt/--no-nt", default=True, show_default=True, help="Diff NT baselines.")
@click.option("--aa/--no-aa", default=True, show_default=True, help="Diff AA baselines.")
def diff_cmd(suffix_a: str, suffix_b: str, nt: bool, aa: bool) -> None:
    """Compare two stored baseline sets identified by SUFFIX_A and SUFFIX_B ('' = unsuffixed)."""
    failed = False
    for enabled, label, base_path in ((nt, "NT", BASELINE_NT), (aa, "AA", BASELINE_AA)):
        if not enabled:
            continue
        path_a = _with_suffix(base_path, suffix_a)
        path_b = _with_suffix(base_path, suffix_b)
        click.echo(f"Diffing {label}: {path_a} vs {path_b}")
        problems = compare_baselines(path_a, path_b)
        if problems:
            failed = True
            click.echo(f"{label} mismatches:", err=True)
            for problem in problems:
                click.echo(f"  {problem}", err=True)
        else:
            click.echo(f"{label} OK")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    cli()
