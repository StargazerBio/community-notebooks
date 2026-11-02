"""Print the notebook folders this CI run checks, as a GitHub Actions output.

A pull request or a push checks the folders it changed; a scheduled or manual
run, or a change to the CI itself, checks every folder. A folder is
`<domain>/<slug>/`, holding `notebook.py`.
"""

import json
import os
import subprocess
from pathlib import Path


def notebook_folders() -> list[str]:
    """Every folder in the repo that holds a notebook."""
    return sorted(str(p.parent) for p in Path().glob("*/*/notebook.py"))


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


def folders_to_check(event: str, base: str) -> list[str]:
    """The folders the run checks for `event`, diffing against `base`."""
    every = notebook_folders()
    if event not in ("pull_request", "push"):
        return every
    paths = changed_paths(base)
    if paths is None or any(p.startswith(("ci/", ".github/")) for p in paths):
        return every
    touched = {"/".join(Path(p).parts[:2]) for p in paths if len(Path(p).parts) > 2}
    return [f for f in every if f in touched]


if __name__ == "__main__":
    folders = folders_to_check(os.environ["EVENT"], os.environ.get("BASE", ""))
    print(f"folders={json.dumps(folders)}")
