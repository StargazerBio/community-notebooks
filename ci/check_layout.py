"""Check the repo's layout: every notebook and the catalog are well formed.

Each notebook is `<domain>/<slug>.py`. Its setup block declares its dashboard
tile as plain strings (`main_name`, `main_description`), and the
`[tool.stargazer]` table in its script header names its `owners`. Every
catalog entry names an existing notebook at a full commit. Prints one line per problem and exits 1 if there are any.
"""

import ast
import re
import sys
import tomllib
from pathlib import Path

import yaml

SLUG = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
# PEP 723's reference regex for an inline script metadata block.
SCRIPT_BLOCK = re.compile(
    r"(?m)^# /// (?P<type>[a-zA-Z0-9-]+)$\s(?P<content>(^#(| .*)$\s)+)^# ///$"
)


def notebooks() -> list[Path]:
    """Every notebook in the repo: a `.py` one level down, outside `ci/`."""
    return sorted(p for p in Path().glob("*/*.py") if p.parts[0] != "ci")


def stargazer_table(notebook: Path) -> dict:
    """The `[tool.stargazer]` table in the notebook's script header, or {}."""
    for match in SCRIPT_BLOCK.finditer(notebook.read_text()):
        if match["type"] != "script":
            continue
        toml = "".join(
            line[2:] if line.startswith("# ") else line[1:]
            for line in match["content"].splitlines(keepends=True)
        )
        return tomllib.loads(toml).get("tool", {}).get("stargazer", {})
    return {}


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


def notebook_problems(notebook: Path) -> list[str]:
    """What's wrong with one `<domain>/<slug>.py` notebook."""
    problems = []
    for part in (notebook.parent.name, notebook.stem):
        if not SLUG.match(part):
            problems.append(f"{notebook}: {part!r} isn't lowercase kebab-case")
    try:
        meta = stargazer_table(notebook)
    except tomllib.TOMLDecodeError as e:
        return [*problems, f"{notebook}: script header isn't valid TOML: {e}"]
    if not meta.get("owners"):
        problems.append(f"{notebook}: [tool.stargazer] names no owners")
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
        if Path(path) not in notebooks():
            problems.append(f"{catalog}: {path!r} isn't a notebook")
        if not COMMIT.match(commit):
            problems.append(f"{catalog}: {path!r} needs a full commit SHA")
    return problems


def main() -> int:
    """Check every notebook and the catalog; print the problems, return the exit code."""
    problems = []
    for notebook in notebooks():
        problems += notebook_problems(notebook)
    for stray in sorted(Path().glob("*/*/**/*.py")):
        problems.append(f"{stray}: a notebook is a file, <domain>/<slug>.py")
    problems += catalog_problems(Path("catalog.yaml"))
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
