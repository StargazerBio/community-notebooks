# Stargazer Community Notebooks

Pipeline notebooks for [Stargazer](https://github.com/StargazerBio/stargazer): parameterized [marimo](https://marimo.io/) notebooks that run a full analysis on data from your asset graph. Stargazer's own pipelines live here beside everyone else's.

| Notebook | What it does |
|---|---|
| [`scrna/clustering-pipeline`](scrna/clustering-pipeline/notebook.py) | scRNA-seq: parallel per-sample preprocessing, then interactive Leiden clustering with side-by-side UMAPs and marker genes |
| [`germline/short-variant-discovery`](germline/short-variant-discovery/notebook.py) | GATK best practices from paired reads to a joint-called cohort VCF |

## Running a notebook

On your Stargazer dashboard, paste a link to a notebook's `notebook.py` into **Open from GitHub**. A link to a commit (the **y** key on GitHub turns a page's link into one) always opens the same version. The notebook lands in your workspace as your own copy.

Locally, in the Stargazer image:

```bash
curl -LO https://raw.githubusercontent.com/StargazerBio/community-notebooks/main/germline/short-variant-discovery/notebook.py
docker run -p 8080:8080 -v "$PWD":/workspace ghcr.io/stargazerbio/stargazer
```

and open `notebook.py` from marimo's file explorer. Each notebook fetches its demo inputs when your graph has none.

## Layout

One folder per notebook, by topic: `<domain>/<slug>/`, both lowercase kebab-case.

- `notebook.py` is the notebook. Its setup block declares the dashboard tile as plain strings, `main_name` and `main_description`, and it imports only public `stargazer` APIs.
- Work a notebook holds behind a button runs anyway when the notebook is started with `--run` (`"run" in mo.cli_args()`), which is how CI starts it.
- `meta.yaml` names the bundle that is its demo input (`demo_bundle`) and who may update it (`owners`).

`catalog.yaml` lists the featured notebooks, each at the commit its maintainers reviewed.

## Contributing

Open a PR adding or changing a folder. CI runs every notebook the PR touches headless in the latest Stargazer image, on an empty asset store holding only its demo bundle, and checks the layout above. A PR whose checks pass merges. Every notebook is also run daily against the latest Stargazer release.

To check a notebook yourself before opening the PR, from the repo root:

```bash
mkdir -m 777 -p out
```

```bash
docker run --rm --platform linux/amd64 -v "$PWD":/repo:ro -v "$PWD/out":/out --entrypoint /repo/ci/run_notebook.sh ghcr.io/stargazerbio/stargazer /repo/<domain>/<slug> /out
```

The page it rendered and its log land in `out/`.

To remove a notebook, open a PR deleting its folder. Its history stays in the repo.

## Licence

[MIT](LICENSE). Contributions are published under it.
