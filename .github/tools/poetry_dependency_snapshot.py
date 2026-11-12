"""Build a GitHub dependency submission snapshot from ``poetry.lock``.

GitHub's static ``poetry.lock`` parser reports every package as ``runtime``, so dev-only packages
(pre-commit -> virtualenv -> filelock, ...) fail ``actions/dependency-review-action`` on advisories that
never reach users. User submissions take priority over static parsing for the same manifest, so we
submit the ``poetry.lock`` manifest ourselves with accurate scopes:

* ``runtime``     - the package belongs to the ``main`` group (installed for users of the package),
* ``development`` - it belongs only to other groups (``dev``, ``release``, ...).

Usage (the snapshot is written to stdout)::

    python .github/tools/poetry_dependency_snapshot.py --sha SHA --ref REF --job-id ID | \\
        gh api -X POST repos/OWNER/REPO/dependency-graph/snapshots --input -

Docs: https://docs.github.com/en/rest/dependency-graph/dependency-submission
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    import tomli as tomllib  # type: ignore[no-redef,import-not-found]

MANIFEST = "poetry.lock"
DETECTOR = {
    "name": "riot-na-poetry-lock-scopes",
    "version": "1.0.0",
    "url": "https://github.com/NaturalAntibody/riot_na/blob/master/.github/tools/poetry_dependency_snapshot.py",
}
CORRELATOR = "riot-na-poetry-lock-scopes"


def normalize(name: str) -> str:
    """PEP 503 name normalization."""
    return re.sub(r"[-_.]+", "-", name).lower()


def requirement_name(requirement: str) -> str:
    """Extract the package name from a PEP 508 requirement string."""
    match = re.match(r"\s*([A-Za-z0-9][A-Za-z0-9._-]*)", requirement)
    if not match:
        raise ValueError(f"Cannot parse requirement: {requirement!r}")
    return normalize(match.group(1))


def direct_dependency_names(pyproject: dict) -> set[str]:
    """Names of packages declared directly in pyproject.toml (PEP 621 and Poetry groups)."""
    names: set[str] = set()
    project = pyproject.get("project", {})
    names.update(requirement_name(req) for req in project.get("dependencies", []))
    for requirements in project.get("optional-dependencies", {}).values():
        names.update(requirement_name(req) for req in requirements)
    poetry = pyproject.get("tool", {}).get("poetry", {})
    names.update(normalize(name) for name in poetry.get("dependencies", {}) if name != "python")
    for group in poetry.get("group", {}).values():
        names.update(normalize(name) for name in group.get("dependencies", {}))
    return names


def purl(name: str, version: str) -> str:
    return f"pkg:pypi/{normalize(name)}@{version}"


def build_manifest(lock: dict, pyproject: dict) -> dict:
    packages = lock["package"]
    direct = direct_dependency_names(pyproject)

    versions_by_name = defaultdict(list)
    for package in packages:
        versions_by_name[normalize(package["name"])].append(package["version"])

    resolved = {}
    for package in packages:
        name = normalize(package["name"])
        groups = package.get("groups", [])
        # A package is only a "development" dependency if no runtime group needs it.
        scope = "runtime" if "main" in groups else "development"
        # The lock does not say which version of a multi-version package (python_version markers)
        # a dependent uses, so all versions are listed as children.
        children = sorted(
            purl(dep, version)
            for dep in map(normalize, package.get("dependencies", {}))
            for version in versions_by_name.get(dep, [])
        )
        package_url = purl(name, package["version"])
        resolved[package_url] = {
            "package_url": package_url,
            "relationship": "direct" if name in direct else "indirect",
            "scope": scope,
            "dependencies": children,
        }
    return {"name": MANIFEST, "file": {"source_location": MANIFEST}, "resolved": resolved}


def build_snapshot(lock: dict, pyproject: dict, sha: str, ref: str, job_id: str) -> dict:
    return {
        "version": 0,
        "sha": sha,
        "ref": ref,
        "job": {"correlator": CORRELATOR, "id": job_id},
        "detector": DETECTOR,
        "scanned": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "manifests": {MANIFEST: build_manifest(lock, pyproject)},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sha", required=True, help="commit SHA the snapshot belongs to")
    parser.add_argument("--ref", required=True, help="full git ref, e.g. refs/heads/master")
    parser.add_argument("--job-id", required=True, help="CI run id (job.id)")
    parser.add_argument("--lock", type=Path, default=Path(MANIFEST))
    parser.add_argument("--pyproject", type=Path, default=Path("pyproject.toml"))
    args = parser.parse_args()

    with args.lock.open("rb") as file:
        lock = tomllib.load(file)
    with args.pyproject.open("rb") as file:
        pyproject = tomllib.load(file)

    json.dump(build_snapshot(lock, pyproject, args.sha, args.ref, args.job_id), sys.stdout, indent=2)
    return 0


if __name__ == "__main__":
    sys.exit(main())
