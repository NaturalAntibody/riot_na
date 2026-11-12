"""Helpers shared by the hand-triggered regression baselines (ngs_baseline, multidomain_baseline)."""

from __future__ import annotations

import json
import pickle
import sys
from collections import Counter
from collections.abc import Sequence
from dataclasses import fields
from pathlib import Path
from typing import Any

import click


def with_suffix(path: Path, suffix: str | None) -> Path:
    """Return ``<stem>_<suffix><ext>`` for ``path`` (unchanged when suffix is empty)."""
    if not suffix:
        return path
    return path.with_name(f"{path.stem}_{suffix}{path.suffix}")


def save_pickle(path: Path, payload: dict[str, Any]) -> None:
    """Atomically write ``payload`` to ``path`` as a pickle."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
    tmp_path.replace(path)


def load_pickle(path: Path, format_version: int) -> dict[str, Any]:
    """Load a baseline pickle and check its ``format_version``."""
    if not path.exists():
        raise FileNotFoundError(path)
    with open(path, "rb") as handle:
        payload = pickle.load(handle)
    if payload.get("format_version") != format_version:
        raise ValueError(
            f"Unsupported baseline format_version={payload.get('format_version')!r} (expected {format_version})"
        )
    return payload


def diff_entry(expected: Any, actual: Any) -> dict[str, dict[str, Any]]:
    if type(expected) is not type(actual):
        return {"__type__": {"expected": type(expected).__name__, "actual": type(actual).__name__}}
    diffs: dict[str, dict[str, Any]] = {}
    for field in fields(expected):
        exp_val = getattr(expected, field.name)
        act_val = getattr(actual, field.name)
        if exp_val != act_val:
            diffs[field.name] = {"expected": exp_val, "actual": act_val}
    return diffs


def query_sequence(*entries: Any) -> str | None:
    """Input sequence numbered into the first of ``entries`` that has one (full molecule for segmented entries)."""
    for entry in entries:
        for name in ("query_sequence", "sequence", "sequence_aa"):
            if value := getattr(entry, name, None):
                return value
    return None


def make_result(expected_count: int, actual_count: int, mismatches: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "ok": not mismatches and expected_count == actual_count,
        "n_expected": expected_count,
        "n_actual": actual_count,
        "n_mismatched_entries": len(mismatches),
        "mismatches": mismatches,
    }


def compare_entries(expected: Sequence[Any], actual: Sequence[Any]) -> dict[str, Any]:
    """Compare entry lists and return a result with counts and every per-entry mismatch."""
    mismatches: list[dict[str, Any]] = []
    for idx, (exp, act) in enumerate(zip(expected, actual)):
        if exp == act:
            continue
        header = getattr(exp, "sequence_header", None) or getattr(act, "sequence_header", None)
        mismatches.append(
            {
                "index": idx,
                "sequence_header": header,
                "query_sequence": query_sequence(exp, act),
                "diffs": diff_entry(exp, act),
            }
        )
    return make_result(len(expected), len(actual), mismatches)


def summarize_result(result: dict[str, Any]) -> list[str]:
    """Return summary lines for a result: mismatch count and mismatches per differing field."""
    field_counts = Counter(name for mismatch in result["mismatches"] for name in mismatch["diffs"])
    lines = [f"{result['n_mismatched_entries']} mismatches"]
    if result["n_expected"] != result["n_actual"]:
        lines.append(f"entry count differs: expected={result['n_expected']} actual={result['n_actual']}")
    lines.extend(f"  {name}: {count}" for name, count in field_counts.most_common())
    return lines


def write_results_jsonl(path: Path, results: Sequence[dict[str, Any]]) -> None:
    """Write one record per mismatched entry of every dataset result."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        for result in results:
            for mismatch in result["mismatches"]:
                record = {"dataset": result["dataset"], **mismatch}
                handle.write(json.dumps(record, default=str) + "\n")


def report(label: str, result: dict[str, Any]) -> None:
    if result["ok"]:
        click.echo(f"{label} OK")
        return
    first, *rest = summarize_result(result)
    click.echo(f"{label}: {first}")
    for line in rest:
        click.echo(line)


def finish(results: Sequence[dict[str, Any]], output_path: Path) -> None:
    """Write mismatches to ``output_path`` and exit with status 1, or report that there are none."""
    if all(result["ok"] for result in results):
        click.echo("No mismatches found.")
        return
    write_results_jsonl(output_path, results)
    click.echo(f"Mismatches found! See {output_path}")
    sys.exit(1)


def output_path_option(default: Path):
    return click.option(
        "--output-path",
        type=click.Path(path_type=Path, dir_okay=False),
        default=default,
        show_default=True,
        help="JSONL file to write comparison results to when mismatches are found (overwritten).",
    )


def ncpu_option(function):
    return click.option(
        "--ncpu",
        type=int,
        default=None,
        show_default="physical cores",
        help="Worker processes for numbering.",
    )(function)


def limit_option(function):
    return click.option(
        "--limit",
        type=int,
        default=None,
        show_default="all",
        help="Only process the first N sequences (debug).",
    )(function)
