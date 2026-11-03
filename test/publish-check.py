# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "stargazer-bio",
# ]
#
# [tool.uv.sources]
# stargazer-bio = { path = "/stargazer", editable = true }
#
# [tool.stargazer]
# owners = ["387300641116005877"]
# ///
"""
### Blank workspace notebook.

An empty marimo notebook for authoring from scratch. `/workspace/create`
copies this seed under your chosen name and declares it as `main_name` in the
setup block. Edit `main_env` there to rightsize the pod for your workload.

spec: [docs/architecture/notebook.md](../architecture/notebook.md)
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium")

with app.setup:
    import flyte

    from stargazer.notebooks import notebook_env

    main_name = "publish-check"

    main_description = "Republished version."
    main_env = notebook_env.clone_with(
        name="notebook", resources=flyte.Resources(cpu=2, memory="4Gi")
    )


@app.cell
def _():
    """Imports."""
    import marimo as mo

    return (mo,)


if __name__ == "__main__":
    app.run()
