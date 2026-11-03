"""Print the notebooks this CI run checks, as a GitHub Actions output.

A pull request or a push checks the notebooks it changed; a scheduled or manual
run, or a change to the CI itself, checks every notebook. A notebook is
`<domain>/<slug>.py`.
"""

import json
import os
import subprocess
from pathlib import Path


def notebooks() -> list[str]:
    """Every notebook in the repo: a `.py` one level down, outside `ci/`."""
    return sorted(str(p) for p in Path().glob("*/*.py") if p.parts[0] != "ci")


def changed_paths(base: str) -> list[str] | None:
    """Paths changed between `base` and HEAD, or None when there's no base to diff."""
    if not base or set(base) == {"0"}:
        return None
    out = subprocess.run(
        ["git", "diff", "--name-only", base, "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    return out.split()


def notebooks_to_check(event: str, base: str) -> list[str]:
    """The notebooks the run checks for `event`, diffing against `base`."""
    every = notebooks()
    if event not in ("pull_request", "push"):
        return every
    paths = changed_paths(base)
    if paths is None or any(p.startswith(("ci/", ".github/")) for p in paths):
        return every
    return [n for n in every if n in paths]


if __name__ == "__main__":
    checked = notebooks_to_check(os.environ["EVENT"], os.environ.get("BASE", ""))
    print(f"notebooks={json.dumps(checked)}")
