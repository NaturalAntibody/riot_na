"""Therapeutics AIRR baseline I/O and numbering helpers for hand-triggered regression.

Numbers therapeutic molecules (data/therapeutics.csv, AA: heavy chains with constant regions, scFvs, bispecifics,
fusions) in both modes:

- single domain: one AirrRearrangementEntryAA per molecule,
- multiple domains: a list of domains (SegmentedAirrRearrangementEntryAA) per molecule.

Binary format: a single pickle blob with metadata + entries of both modes.

Generate baseline::

    poetry run python -m tests.regression.multidomain_baseline generate --ncpu 8
    poetry run python -m tests.regression.multidomain_baseline generate --suffix ssw  # -> baseline_airr_therapeutics_aa_ssw.pkl

Compare against baseline after changing the code::

    poetry run python -m tests.regression.multidomain_baseline compare --ncpu 8
    poetry run python -m tests.regression.multidomain_baseline compare --suffix ssw

Compare two stored baselines by suffix, without re-numbering ('' = unsuffixed file)::

    poetry run python -m tests.regression.multidomain_baseline diff ssw parasail
    poetry run python -m tests.regression.multidomain_baseline diff '' parasail --no-single

When mismatches are found, ``compare`` and ``diff`` write results to JSONL (default
``multidomain_regression_result.jsonl`` next to the baseline, override with ``--output-path``): one record
per differing domain (with ``diffs``). Mismatch counts are printed to stdout.
"""

from __future__ import annotations

import csv
import itertools
import multiprocessing as mp
from collections import Counter
from collections.abc import Sequence
from pathlib import Path
from typing import Any, cast

import click
from psutil import cpu_count
from tqdm import tqdm

from riot_na.api.riot_numbering import RiotNumberingAA, create_riot_aa
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    Scheme,
    SegmentedAirrRearrangementEntryAA,
)

from .common import (
    diff_entry,
    finish,
    load_pickle,
    make_result,
    ncpu_option,
    output_path_option,
    query_sequence,
    report,
    save_pickle,
    with_suffix,
)

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
THERAPEUTICS_CSV = DATA_DIR / "therapeutics.csv"
BASELINE = DATA_DIR / "baseline_airr_therapeutics_aa.pkl"
RESULT_JSONL = DATA_DIR / "multidomain_regression_result.jsonl"

BASELINE_FORMAT_VERSION = 2
SCHEME = Scheme.IMGT

Molecule = tuple[str, str]  # (header, sequence)
MoleculeDomains = list[SegmentedAirrRearrangementEntryAA]


def read_therapeutics(csv_path: Path = THERAPEUTICS_CSV, limit: int | None = None) -> list[Molecule]:
    """Read (header, sequence) of therapeutic molecules, header is `therapeutic|molecule_name`."""
    with open(csv_path, newline="") as handle:
        rows = csv.DictReader(handle)
        return [
            (f"{row['therapeutic']}|{row['molecule_name']}", row["molecule_sequence"])
            for row in itertools.islice(rows, limit)
        ]


_single_domain_numbering: RiotNumberingAA | None = None
_multidomain_numbering: RiotNumberingAA | None = None


def _init_single_domain_worker() -> None:
    global _single_domain_numbering  # pylint: disable=global-statement
    _single_domain_numbering = create_riot_aa()


def _init_multidomain_worker() -> None:
    global _multidomain_numbering  # pylint: disable=global-statement
    _multidomain_numbering = create_riot_aa(return_all_domains=True)


def _number_single_domain(molecule: Molecule) -> AirrRearrangementEntryAA:
    assert _single_domain_numbering is not None
    header, sequence = molecule
    return _single_domain_numbering.run_on_sequence(header, sequence, SCHEME)


def _number_multidomain(molecule: Molecule) -> MoleculeDomains:
    assert _multidomain_numbering is not None
    header, sequence = molecule
    domains = _multidomain_numbering.run_on_sequence(header, sequence, SCHEME, return_all_domains=True)
    return cast(MoleculeDomains, domains)  # the builder returns segmented entries, the API is typed more loosely


def number_single_domain(
    molecules: Sequence[Molecule], n_processes: int | None = None
) -> list[AirrRearrangementEntryAA]:
    """Number molecules in single domain mode, returns the best domain of every molecule."""
    with mp.Pool(processes=n_processes or cpu_count(logical=False), initializer=_init_single_domain_worker) as pool:
        return list(
            tqdm(
                pool.imap(_number_single_domain, molecules, chunksize=8),
                total=len(molecules),
                desc="numbering therapeutics (single domain)",
            )
        )


def number_multidomain(molecules: Sequence[Molecule], n_processes: int | None = None) -> list[MoleculeDomains]:
    """Number molecules in multiple domains mode, returns all domains of every molecule."""
    with mp.Pool(processes=n_processes or cpu_count(logical=False), initializer=_init_multidomain_worker) as pool:
        return list(
            tqdm(
                pool.imap(_number_multidomain, molecules, chunksize=8),
                total=len(molecules),
                desc="numbering therapeutics (multidomain)",
            )
        )


def save_baseline(
    path: Path,
    molecules: Sequence[Molecule],
    single_domain: Sequence[AirrRearrangementEntryAA],
    multidomain: Sequence[MoleculeDomains],
) -> None:
    payload = {
        "format_version": BASELINE_FORMAT_VERSION,
        "scheme": SCHEME.value,
        "headers": [header for header, _ in molecules],
        "single_domain": list(single_domain),
        "multidomain": list(multidomain),
    }
    save_pickle(path, payload)


def load_baseline(path: Path) -> dict[str, Any]:
    return load_pickle(path, BASELINE_FORMAT_VERSION)


def compare_molecules(
    headers: Sequence[str],
    expected: Sequence[Sequence[Any]],
    actual: Sequence[Sequence[Any]],
) -> dict[str, Any]:
    """Compare domains of every molecule, return a result in the format of common.compare_entries."""
    mismatches: list[dict[str, Any]] = []
    for idx, (header, exp_domains, act_domains) in enumerate(zip(headers, expected, actual)):
        if len(exp_domains) != len(act_domains):
            diffs = {"__domain_count__": {"expected": len(exp_domains), "actual": len(act_domains)}}
            mismatches.append(
                {
                    "index": idx,
                    "sequence_header": header,
                    "query_sequence": query_sequence(*exp_domains, *act_domains),
                    "diffs": diffs,
                }
            )
            continue
        for domain_idx, (exp, act) in enumerate(zip(exp_domains, act_domains)):
            if exp != act:
                mismatches.append(
                    {
                        "index": idx,
                        "sequence_header": header,
                        "query_sequence": query_sequence(exp, act),
                        "domain": domain_idx,
                        "diffs": diff_entry(exp, act),
                    }
                )

    return make_result(len(expected), len(actual), mismatches)


def _mode_results(
    headers: Sequence[str],
    expected: dict[str, Any],
    actual: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compare every mode present in ``actual`` (single_domain / multidomain payload keys) and report it."""
    results = []
    for key, label in (("single_domain", "Single domain"), ("multidomain", "Multidomain")):
        if key not in actual:
            continue
        if key == "single_domain":
            result = compare_molecules(headers, [[e] for e in expected[key]], [[a] for a in actual[key]])
            for mismatch in result["mismatches"]:
                mismatch.pop("domain", None)  # always 0 in single domain mode
        else:
            result = compare_molecules(headers, expected[key], actual[key])
        results.append({"dataset": label, **result})
        report(label, result)
    return results


@click.group()
def cli() -> None:
    """Therapeutics AIRR regression baseline tools."""


@cli.command("generate")
@click.option("--limit", type=int, default=None, show_default="all", help="Only process the first N molecules.")
@ncpu_option
@click.option(
    "--suffix",
    type=str,
    default=None,
    help="Append '_<suffix>' to baseline file name, e.g. baseline_airr_therapeutics_aa_<suffix>.pkl.",
)
def generate_cmd(limit: int | None, ncpu: int | None, suffix: str = "") -> None:
    """Generate single domain and multidomain AA baselines."""
    baseline_path = with_suffix(BASELINE, suffix)
    molecules = read_therapeutics(limit=limit)
    single_domain = number_single_domain(molecules, n_processes=ncpu)
    multidomain = number_multidomain(molecules, n_processes=ncpu)
    click.echo(f"Writing therapeutics baseline -> {baseline_path}")
    save_baseline(baseline_path, molecules, single_domain, multidomain)
    click.echo(f"Done: {len(molecules)} molecules, {sum(map(len, multidomain))} multidomain domains.")


@cli.command("compare")
@ncpu_option
@click.option("--single/--no-single", default=True, show_default=True, help="Compare single domain mode.")
@click.option("--multidomain/--no-multidomain", default=True, show_default=True, help="Compare multidomain mode.")
@click.option(
    "--suffix",
    type=str,
    default=None,
    help="Compare against baseline generated with the same '--suffix'.",
)
@output_path_option(RESULT_JSONL)
def compare_cmd(ncpu: int | None, single: bool, multidomain: bool, output_path: Path, suffix: str = "") -> None:
    """Re-run numbering and compare against pickled baseline."""
    baseline_path = with_suffix(BASELINE, suffix)
    baseline = load_baseline(baseline_path)
    headers = baseline["headers"]
    molecules = read_therapeutics(limit=len(headers))
    if [header for header, _ in molecules] != headers:
        raise ValueError(f"{THERAPEUTICS_CSV} molecules differ from the baseline")

    click.echo(f"Comparing therapeutics: {THERAPEUTICS_CSV} vs {baseline_path}")
    actual: dict[str, Any] = {}
    if single:
        actual["single_domain"] = number_single_domain(molecules, n_processes=ncpu)
    if multidomain:
        actual["multidomain"] = number_multidomain(molecules, n_processes=ncpu)
    finish(_mode_results(headers, baseline, actual), output_path)


@cli.command("diff")
@click.argument("suffix_a")
@click.argument("suffix_b")
@click.option("--single/--no-single", default=True, show_default=True, help="Diff single domain mode.")
@click.option("--multidomain/--no-multidomain", default=True, show_default=True, help="Diff multidomain mode.")
@output_path_option(RESULT_JSONL)
def diff_cmd(suffix_a: str, suffix_b: str, single: bool, multidomain: bool, output_path: Path) -> None:
    """Compare two stored baselines identified by SUFFIX_A and SUFFIX_B ('' = unsuffixed).

    SUFFIX_A is reported as "expected", SUFFIX_B as "actual".
    """
    path_a = with_suffix(BASELINE, suffix_a)
    path_b = with_suffix(BASELINE, suffix_b)
    payload_a = load_baseline(path_a)
    payload_b = load_baseline(path_b)
    for key in ("scheme", "headers"):
        if payload_a[key] != payload_b[key]:
            raise ValueError(f"{key} differs between {path_a.name} and {path_b.name}")

    click.echo(f"Diffing therapeutics: {path_a} vs {path_b}")
    modes = [key for key, enabled in (("single_domain", single), ("multidomain", multidomain)) if enabled]
    actual = {key: payload_b[key] for key in modes}
    finish(_mode_results(payload_a["headers"], payload_a, actual), output_path)


MODE_COMPARISON_FIELDS = (
    "locus",
    "v_call",
    "j_call",
    "c_call",
    "fwr1_aa",
    "cdr1_aa",
    "fwr2_aa",
    "cdr2_aa",
    "fwr3_aa",
    "cdr3_aa",
    "fwr4_aa",
    "scheme_residue_mapping",
)


def _short(value: Any, max_length: int = 60) -> str:
    text = repr(value)
    return text if len(text) <= max_length else text[: max_length - 3] + "..."


def _matching_domain(
    single_domain: AirrRearrangementEntryAA, domains: MoleculeDomains
) -> SegmentedAirrRearrangementEntryAA | None:
    """Domain of the same V gene as the single domain result, or the best scoring domain."""
    domains_with_v = [domain for domain in domains if domain.v_call]
    same_v = [domain for domain in domains_with_v if domain.v_call == single_domain.v_call]
    return max(same_v or domains_with_v, default=None)


def compare_modes(
    headers: Sequence[str],
    single_domain: Sequence[AirrRearrangementEntryAA],
    multidomain: Sequence[MoleculeDomains],
    *,
    max_examples: int = 10,
) -> list[str]:
    """
    Compare single domain results with the matching domain of multidomain results (domain-independent fields:
    gene calls, regions and numbering). Returns a human-readable summary.
    """
    outcomes: Counter[str] = Counter()
    field_differences: Counter[str] = Counter()
    examples: list[str] = []

    for header, single, domains in zip(headers, single_domain, multidomain):
        domain_count = "1 domain" if len(domains) == 1 else "2+ domains"
        domain = _matching_domain(single, domains)
        if not single.v_call or domain is None:
            outcome = {
                (False, True): "no V in both modes",
                (True, True): "V in single domain mode only",
                (False, False): "V in multidomain mode only",
            }[(bool(single.v_call), domain is None)]
            outcomes[outcome] += 1
            continue

        differing = [name for name in MODE_COMPARISON_FIELDS if getattr(single, name) != getattr(domain, name)]
        outcomes[f"{domain_count}: {'same' if not differing else 'different'}"] += 1
        field_differences.update(differing)
        if differing and len(examples) < max_examples:
            examples.append(
                f"{header!r} ({domain_count}): "
                + "; ".join(
                    f"{name}: {_short(getattr(single, name))} vs {_short(getattr(domain, name))}"
                    for name in differing[:3]
                )
            )

    return [
        f"{len(headers)} molecules",
        *(f"  {outcome}: {count}" for outcome, count in sorted(outcomes.items())),
        "differing fields (single domain vs matching domain):",
        *(f"  {name}: {count}" for name, count in field_differences.most_common()),
        "examples:",
        *(f"  {example}" for example in examples),
    ]


@cli.command("compare-modes")
@ncpu_option
@click.option("--baseline/--no-baseline", default=False, show_default=True, help="Use results stored in baseline.")
@click.option("--suffix", type=str, default=None, help="With --baseline, use baseline generated with this '--suffix'.")
def compare_modes_cmd(ncpu: int | None, baseline: bool, suffix: str = "") -> None:
    """Compare single domain and multidomain results of the same molecules."""
    if baseline:
        payload = load_baseline(with_suffix(BASELINE, suffix))
        headers, single_domain, multidomain = payload["headers"], payload["single_domain"], payload["multidomain"]
    else:
        molecules = read_therapeutics()
        headers = [header for header, _ in molecules]
        single_domain = number_single_domain(molecules, n_processes=ncpu)
        multidomain = number_multidomain(molecules, n_processes=ncpu)

    for line in compare_modes(headers, single_domain, multidomain):
        click.echo(line)


if __name__ == "__main__":
    cli()
