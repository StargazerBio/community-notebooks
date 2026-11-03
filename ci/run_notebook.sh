#!/usr/bin/env bash
# Run one notebook headless the way the dashboard would, inside the Stargazer
# image: `docker run --entrypoint /repo/ci/run_notebook.sh <image> <notebook> <out>`.
#
# The run gets an asset store of its own, empty but for the notebook's demo
# bundle (`demo_bundle` in its script header's `[tool.stargazer]`), and no
# Pinata key, so the notebook sees exactly its demo inputs. `marimo export html --sandbox` builds the
# notebook's venv from its script header, as a hosted pod does, runs every cell
# and writes the page to <out>/<domain>--<slug>.html. It exits 1 if any cell
# raised; the log beside the page names the exception. The notebook gets
# `--run` (`mo.cli_args()`), so work it holds behind a button runs too.
set -euo pipefail

NOTEBOOK="$1"
OUT="$2"
NAME="$(basename "$(dirname "$NOTEBOOK")")--$(basename "$NOTEBOOK" .py)"

RUN_DIR="$(mktemp -d)"
export STARGAZER_LOCAL="$RUN_DIR/local"
export STARGAZER_STORE_ROOT="$RUN_DIR/store"
export STARGAZER_INDEX_URL="$RUN_DIR/index.db"
export PINATA_JWT=

# ci/check_layout.py checks the header is valid TOML, so one line is enough here.
BUNDLE="$(sed -n 's/^# demo_bundle[[:space:]]*=[[:space:]]*"\([A-Za-z0-9_-]*\)".*/\1/p' "$NOTEBOOK")"
echo "== $NOTEBOOK: fetching demo bundle $BUNDLE"
# The image's runtime venv holds no SDK, so the fetch runs in one uv builds from
# the image's SDK tree, the same tree each notebook's --sandbox venv builds from.
uv run --quiet --no-project --with /stargazer python -c "import asyncio, sys; from stargazer.bundles import fetch_bundle; asyncio.run(fetch_bundle(sys.argv[1]))" "$BUNDLE"

# marimo writes next to the notebook and the checkout is read-only.
mkdir "$RUN_DIR/notebook"
cp "$NOTEBOOK" "$RUN_DIR/notebook/"
cd "$RUN_DIR/notebook"

echo "== $NOTEBOOK: running every cell"
status=0
marimo -q export html --sandbox "$(basename "$NOTEBOOK")" -o "$OUT/$NAME.html" -- --run \
  > "$OUT/$NAME.log" 2>&1 || status=$?
cat "$OUT/$NAME.log"
exit "$status"
