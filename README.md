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
- Its script header, beside the dependencies, has a `[tool.stargazer]` table naming who may update it (`owners`, Union subjects, or `maintainers` for the org's admins):

  ```python
  # [tool.stargazer]
  # owners = ["maintainers"]
  ```

`catalog.yaml` lists the featured notebooks, each at the commit its maintainers reviewed.

## Publishing

Publish a notebook from its tile on your Stargazer dashboard: give it a name and a domain, and it's committed here as yours. Republishing commits the new version; unpublishing deletes it, though its history stays in the repo. You can also open a PR adding a notebook by hand.

CI checks the layout above on every push and PR.

## Licence

[MIT](LICENSE). Contributions are published under it.
