# Stargazer Community Notebooks

Pipeline notebooks for [Stargazer](https://github.com/StargazerBio/stargazer): parameterized [marimo](https://marimo.io/) notebooks that run a full analysis on data from your asset graph. Stargazer's own pipelines live here beside everyone else's.

| Notebook | What it does |
|---|---|
| [`scrna/clustering-pipeline`](scrna/clustering-pipeline.py) | scRNA-seq: parallel per-sample preprocessing, then interactive Leiden clustering with side-by-side UMAPs and marker genes |
| [`germline/short-variant-discovery`](germline/short-variant-discovery.py) | GATK best practices from paired reads to a joint-called cohort VCF |

## Running a notebook

On your Stargazer dashboard, paste a link to a notebook into **Open from GitHub**. A link to a commit (the **y** key on GitHub turns a page's link into one) always opens the same version. The notebook lands in your workspace as your own copy.

Locally, in the Stargazer image:

```bash
curl -LO https://raw.githubusercontent.com/StargazerBio/community-notebooks/main/germline/short-variant-discovery.py
docker run -p 8080:8080 -v "$PWD":/workspace ghcr.io/stargazerbio/stargazer
```

and open `short-variant-discovery.py` from marimo's file explorer. Each notebook fetches its demo inputs when your graph has none.

## Layout

One file per notebook, by topic: `<domain>/<slug>.py`, both lowercase kebab-case. Everything about a notebook is in the notebook:

- Its setup block declares the dashboard tile as plain strings, `main_name` and `main_description`, and it imports only public `stargazer` APIs.
- Its script header, beside the dependencies, has a `[tool.stargazer]` table naming the bundle that is its demo input (`demo_bundle`) and who may update it (`owners`):

  ```python
  # [tool.stargazer]
  # owners = ["maintainers"]
  # demo_bundle = "variant_calling_demo"
  ```

- Work a notebook holds behind a button runs anyway when the notebook is started with `--run` (`"run" in mo.cli_args()`), which is how CI starts it.

`catalog.yaml` lists the featured notebooks, each at the commit its maintainers reviewed.

## Contributing

Open a PR adding or changing a notebook. CI runs every notebook the PR touches headless in the latest Stargazer image, on an empty asset store holding only its demo bundle, and checks the layout above. A PR whose checks pass merges. Every notebook is also run daily against the latest Stargazer release.

To check a notebook yourself before opening the PR, from the repo root:

```bash
mkdir -m 777 -p out
```

```bash
docker run --rm --platform linux/amd64 -v "$PWD":/repo:ro -v "$PWD/out":/out --entrypoint /repo/ci/run_notebook.sh ghcr.io/stargazerbio/stargazer /repo/<domain>/<slug>.py /out
```

The page it rendered and its log land in `out/`.

To remove a notebook, open a PR deleting it. Its history stays in the repo.

## Licence

[MIT](LICENSE). Contributions are published under it.
