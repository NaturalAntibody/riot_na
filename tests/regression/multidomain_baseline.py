"""Therapeutics AIRR baseline I/O and numbering helpers for hand-triggered regression.

Numbers therapeutic molecules (data/therapeutics.csv, AA: heavy chains with constant regions, scFvs, bispecifics,
fusions) in both modes:

- single domain: one AirrRearrangementEntryAA per molecule,
- multiple domains: a list of domains (SegmentedAirrRearrangementEntryAA) per molecule.

Binary format: a single pickle blob with metadata + entries of both modes.

Generate baseline::

    poetry run python -m tests.regression.multidomain_baseline generate --ncpu 8

Compare against baseline after changing the code::

    poetry run python -m tests.regression.multidomain_baseline compare --ncpu 8
"""

from __future__ import annotations

import csv
import itertools
import multiprocessing as mp
import pickle
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Iterator, Optional, Sequence

import click
from psutil import cpu_count
from tqdm import tqdm

from riot_na.api.multidomain_riot_numbering import (
    MultidomainRiotNumberingAA,
    create_multidomain_riot_aa,
)
from riot_na.api.riot_numbering import RiotNumberingAA, create_riot_aa
from riot_na.data.model import (
    AirrRearrangementEntryAA,
    CScheme,
    Scheme,
    SegmentedAirrRearrangementEntryAA,
)
from tests.regression.ngs_baseline import _diff_entry

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
THERAPEUTICS_CSV = DATA_DIR / "therapeutics.csv"
BASELINE = DATA_DIR / "baseline_airr_therapeutics_aa.pkl"

BASELINE_FORMAT_VERSION = 2
SCHEME = Scheme.IMGT
C_SCHEME = CScheme.EU

Molecule = tuple[str, str]  # (header, sequence)
MoleculeDomains = list[SegmentedAirrRearrangementEntryAA]


def read_therapeutics(csv_path: Path = THERAPEUTICS_CSV, limit: Optional[int] = None) -> list[Molecule]:
    """Read (header, sequence) of therapeutic molecules, header is `therapeutic|molecule_name`."""
    with open(csv_path, newline="") as handle:
        rows = csv.DictReader(handle)
        return [
            (f"{row['therapeutic']}|{row['molecule_name']}", row["molecule_sequence"])
            for row in itertools.islice(rows, limit)
        ]


_single_domain_numbering: Optional[RiotNumberingAA] = None
_multidomain_numbering: Optional[MultidomainRiotNumberingAA] = None


def _init_single_domain_worker() -> None:
    global _single_domain_numbering  # pylint: disable=global-statement
    _single_domain_numbering = create_riot_aa()


def _init_multidomain_worker() -> None:
    global _multidomain_numbering  # pylint: disable=global-statement
    _multidomain_numbering = create_multidomain_riot_aa()


def _number_single_domain(molecule: Molecule) -> AirrRearrangementEntryAA:
    assert _single_domain_numbering is not None
    header, sequence = molecule
    return _single_domain_numbering.run_on_sequence(header, sequence, SCHEME, c_scheme=C_SCHEME)


def _number_multidomain(molecule: Molecule) -> MoleculeDomains:
    assert _multidomain_numbering is not None
    header, sequence = molecule
    return _multidomain_numbering.run_on_sequence(header, sequence, SCHEME, c_scheme=C_SCHEME)


def number_single_domain(
    molecules: Sequence[Molecule], n_processes: Optional[int] = None
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


def number_multidomain(molecules: Sequence[Molecule], n_processes: Optional[int] = None) -> list[MoleculeDomains]:
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
        "c_scheme": C_SCHEME.value,
        "headers": [header for header, _ in molecules],
        "single_domain": list(single_domain),
        "multidomain": list(multidomain),
    }
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with open(tmp_path, "wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
    tmp_path.replace(path)


def load_baseline(path: Path) -> dict[str, Any]:
    with open(path, "rb") as handle:
        payload = pickle.load(handle)
    if payload.get("format_version") != BASELINE_FORMAT_VERSION:
        raise ValueError(
            f"Unsupported baseline format_version={payload.get('format_version')!r} "
            f"(expected {BASELINE_FORMAT_VERSION})"
        )
    return payload


def _molecule_problems(header: str, expected: Sequence[Any], actual: Sequence[Any]) -> Iterator[str]:
    if len(expected) != len(actual):
        yield f"{header!r}: domain count mismatch: baseline={len(expected)} actual={len(actual)}"
        return
    for idx, (exp, act) in enumerate(zip(expected, actual)):
        if exp == act:
            continue
        field_diffs = _diff_entry(exp, act) if type(exp) is type(act) else [f"type {type(exp)} vs {type(act)}"]
        yield f"{header!r} domain {idx}: " + "; ".join(field_diffs[:10])


def compare_molecules(
    headers: Sequence[str],
    expected: Sequence[Sequence[Any]],
    actual: Sequence[Sequence[Any]],
    *,
    max_report: int = 20,
) -> list[str]:
    """Compare domains of every molecule, return human-readable mismatch descriptions (empty if identical)."""
    problems = list(
        itertools.chain.from_iterable(
            _molecule_problems(header, exp, act) for header, exp, act in zip(headers, expected, actual)
        )
    )
    if len(expected) != len(actual):
        problems.insert(0, f"molecule count mismatch: baseline={len(expected)} actual={len(actual)}")
    if len(problems) > max_report:
        problems = problems[:max_report] + [f"... {len(problems) - max_report} more mismatches"]
    return problems


def _report(mode: str, problems: list[str]) -> bool:
    if problems:
        click.echo(f"{mode} mismatches:", err=True)
        for problem in problems:
            click.echo(f"  {problem}", err=True)
        return False
    click.echo(f"{mode} OK")
    return True


@click.group()
def cli() -> None:
    """Therapeutics AIRR regression baseline tools."""


@cli.command("generate")
@click.option("--limit", type=int, default=None, show_default="all", help="Only process the first N molecules.")
@click.option("--ncpu", type=int, default=None, show_default="physical cores", help="Worker processes.")
def generate_cmd(limit: Optional[int], ncpu: Optional[int]) -> None:
    """Generate single domain and multidomain AA baselines."""
    molecules = read_therapeutics(limit=limit)
    single_domain = number_single_domain(molecules, n_processes=ncpu)
    multidomain = number_multidomain(molecules, n_processes=ncpu)
    click.echo(f"Writing therapeutics baseline -> {BASELINE}")
    save_baseline(BASELINE, molecules, single_domain, multidomain)
    click.echo(f"Done: {len(molecules)} molecules, {sum(map(len, multidomain))} multidomain domains.")


@cli.command("compare")
@click.option("--ncpu", type=int, default=None, show_default="physical cores", help="Worker processes.")
@click.option("--single/--no-single", default=True, show_default=True, help="Compare single domain mode.")
@click.option("--multidomain/--no-multidomain", default=True, show_default=True, help="Compare multidomain mode.")
def compare_cmd(ncpu: Optional[int], single: bool, multidomain: bool) -> None:
    """Re-run numbering and compare against pickled baseline."""
    baseline = load_baseline(BASELINE)
    headers = baseline["headers"]
    molecules = read_therapeutics(limit=len(headers))
    if [header for header, _ in molecules] != headers:
        raise ValueError(f"{THERAPEUTICS_CSV} molecules differ from the baseline")

    click.echo(f"Comparing therapeutics: {THERAPEUTICS_CSV} vs {BASELINE}")
    ok = True
    if single:
        actual = number_single_domain(molecules, n_processes=ncpu)
        expected = baseline["single_domain"]
        ok &= _report("Single domain", compare_molecules(headers, [[e] for e in expected], [[a] for a in actual]))
    if multidomain:
        actual_domains = number_multidomain(molecules, n_processes=ncpu)
        ok &= _report("Multidomain", compare_molecules(headers, baseline["multidomain"], actual_domains))
    if not ok:
        sys.exit(1)


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
    "c_scheme_residue_mapping",
)


def _short(value: Any, max_length: int = 60) -> str:
    text = repr(value)
    return text if len(text) <= max_length else text[: max_length - 3] + "..."


def _matching_domain(
    single_domain: AirrRearrangementEntryAA, domains: MoleculeDomains
) -> Optional[SegmentedAirrRearrangementEntryAA]:
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
@click.option("--ncpu", type=int, default=None, show_default="physical cores", help="Worker processes.")
@click.option("--baseline/--no-baseline", default=False, show_default=True, help="Use results stored in baseline.")
def compare_modes_cmd(ncpu: Optional[int], baseline: bool) -> None:
    """Compare single domain and multidomain results of the same molecules."""
    if baseline:
        payload = load_baseline(BASELINE)
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
