#!/usr/bin/env bash
# Run one notebook headless the way the dashboard would, inside the Stargazer
# image: `docker run --entrypoint /repo/ci/run_notebook.sh <image> <folder> <out>`.
#
# The run gets an asset store of its own, empty but for the folder's demo
# bundle (`demo_bundle` in its meta.yaml), and no Pinata key, so the notebook
# sees exactly its demo inputs. `marimo export html --sandbox` builds the
# notebook's venv from its script header, as a hosted pod does, runs every cell
# and writes the page to <out>/<domain>--<slug>.html. It exits 1 if any cell
# raised; the log beside the page names the exception. The notebook gets
# `--run` (`mo.cli_args()`), so work it holds behind a button runs too.
set -euo pipefail

FOLDER="${1%/}"
OUT="$2"
NAME="$(basename "$(dirname "$FOLDER")")--$(basename "$FOLDER")"

RUN_DIR="$(mktemp -d)"
export STARGAZER_LOCAL="$RUN_DIR/local"
export STARGAZER_STORE_ROOT="$RUN_DIR/store"
export STARGAZER_INDEX_URL="$RUN_DIR/index.db"
export PINATA_JWT=

BUNDLE="$(sed -n 's/^demo_bundle:[[:space:]]*\([A-Za-z0-9_-]*\).*/\1/p' "$FOLDER/meta.yaml")"
echo "== $FOLDER: fetching demo bundle $BUNDLE"
# The image's MCP venv has the SDK installed; the notebook's own venv is built
# by --sandbox below.
/opt/mcp/bin/python -c "import asyncio, sys; from stargazer.bundles import fetch_bundle; asyncio.run(fetch_bundle(sys.argv[1]))" "$BUNDLE"

# marimo writes next to the notebook and the checkout is read-only.
cp -r "$FOLDER" "$RUN_DIR/notebook"
cd "$RUN_DIR/notebook"

echo "== $FOLDER: running every cell"
status=0
marimo -q export html --sandbox notebook.py -o "$OUT/$NAME.html" -- --run \
  > "$OUT/$NAME.log" 2>&1 || status=$?
cat "$OUT/$NAME.log"
exit "$status"
