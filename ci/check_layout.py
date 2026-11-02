"""Check the repo's layout: every notebook folder and the catalog are well formed.

Each `<domain>/<slug>/` holds `notebook.py` and `meta.yaml`. The notebook's
setup block declares its dashboard tile as plain strings (`main_name`,
`main_description`), and `meta.yaml` names its `owners` and the `demo_bundle`
CI runs it on. Every catalog entry names an existing folder at a full commit.
Prints one line per problem and exits 1 if there are any.
"""

import ast
import re
import sys
from pathlib import Path

import yaml

SLUG = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")


def setup_strings(notebook: Path) -> dict[str, str]:
    """The plain-string assignments in the notebook's `with app.setup:` block."""
    tree = ast.parse(notebook.read_text())
    found = {}
    for node in ast.walk(tree):
        if not isinstance(node, ast.With):
            continue
        if ast.unparse(node.items[0].context_expr) != "app.setup":
            continue
        for stmt in node.body:
            if (
                isinstance(stmt, ast.Assign)
                and len(stmt.targets) == 1
                and isinstance(stmt.targets[0], ast.Name)
                and isinstance(stmt.value, ast.Constant)
                and isinstance(stmt.value.value, str)
            ):
                found[stmt.targets[0].id] = stmt.value.value
    return found


def folder_problems(folder: Path) -> list[str]:
    """What's wrong with one `<domain>/<slug>/` folder."""
    problems = []
    for part in folder.parts:
        if not SLUG.match(part):
            problems.append(f"{folder}: {part!r} isn't lowercase kebab-case")
    notebook, meta_path = folder / "notebook.py", folder / "meta.yaml"
    if not meta_path.exists():
        return [*problems, f"{folder}: no meta.yaml"]
    meta = yaml.safe_load(meta_path.read_text()) or {}
    if not meta.get("owners"):
        problems.append(f"{meta_path}: no owners")
    if not meta.get("demo_bundle"):
        problems.append(f"{meta_path}: no demo_bundle")
    strings = setup_strings(notebook)
    for name in ("main_name", "main_description"):
        if not strings.get(name):
            problems.append(f"{notebook}: setup block declares no {name}")
    return problems


def catalog_problems(catalog: Path) -> list[str]:
    """What's wrong with the catalog."""
    entries = (yaml.safe_load(catalog.read_text()) or {}).get("notebooks") or []
    problems = []
    for entry in entries:
        path, commit = entry.get("path", ""), entry.get("commit", "")
        if not (Path(path) / "notebook.py").exists():
            problems.append(f"{catalog}: {path!r} isn't a notebook folder")
        if not COMMIT.match(commit):
            problems.append(f"{catalog}: {path!r} needs a full commit SHA")
    return problems


def main() -> int:
    """Check every folder and the catalog; print the problems, return the exit code."""
    problems = []
    for notebook in sorted(Path().glob("*/*/notebook.py")):
        problems += folder_problems(notebook.parent)
    for stray in sorted(Path().glob("*/*/*.py")):
        if stray.name != "notebook.py" and not stray.is_relative_to("ci"):
            problems.append(f"{stray}: a folder's notebook is named notebook.py")
    problems += catalog_problems(Path("catalog.yaml"))
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
