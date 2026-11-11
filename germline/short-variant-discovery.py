# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "marimo",
#   "stargazer-bio",
# ]
#
# [tool.pixi.workspace]
# channels = ["conda-forge", "bioconda"]
#
# [tool.pixi.dependencies]
# python = "3.14.*"
# gatk4 = "4.6.*"
# samtools = "1.*"
# bwa-mem2 = "2.*"
#
# [tool.pixi.pypi-dependencies]
# stargazer-bio = { path = "/stargazer" }
#
# [tool.stargazer]
# # Who may republish or unpublish this notebook from the dashboard, as Union
# # subjects; "maintainers" is the org's admins.
# owners = ["maintainers"]
# ///
"""
### Stargazer germline short variant discovery.

GATK best practices from raw reads to a joint-called cohort VCF, run on a
reference and samples picked from the asset graph. The workflow takes
questions (a build, sample IDs), not files: it assembles its own inputs, so
the same notebook runs on the demo slice of TP53, on your own data, or on
published data, with nothing to change.
"""

import marimo

__generated_with = "0.21.1"
app = marimo.App(width="medium")

with app.setup:
    import gzip
    from pathlib import Path

    import marimo as mo

    from stargazer.assets.asset import assemble
    from stargazer.assets.reads import R1
    from stargazer.assets.reference import Reference
    from stargazer.workflows import germline_short_variant_discovery

    main_name = "Germline variants"
    main_description = "GATK best practices from paired reads to a cohort VCF."

    DEMO_BUNDLE = "variant_calling_demo"


@app.function
def classify_variant(ref: str, alt: str) -> str:
    """The kind of one ALT allele against its REF: SNP, insertion, deletion, MNP."""
    if alt.startswith("<") or alt in ("*", "."):
        return "symbolic"
    if len(ref) == len(alt):
        return "SNP" if len(ref) == 1 else "MNP"
    return "insertion" if len(alt) > len(ref) else "deletion"


@app.function
def read_vcf(path: Path) -> tuple[list[str], list[list[str]]]:
    """The sample names and the data lines (split on tabs) of a VCF, gzipped or not."""
    with path.open("rb") as f:
        gzipped = f.read(2) == b"\x1f\x8b"
    opener = gzip.open if gzipped else open
    samples, records = [], []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = fields[9:]
            else:
                records.append(fields)
    return samples, records


@app.function
def count_variants(records: list[list[str]]) -> dict[str, int]:
    """Count ALT alleles by kind, a multiallelic site counting once per allele."""
    counts: dict[str, int] = {}
    for fields in records:
        for alt in fields[4].split(","):
            kind = classify_variant(fields[3], alt)
            counts[kind] = counts.get(kind, 0) + 1
    return counts


@app.cell
def _():
    """Title and Flyte init."""
    from stargazer.config import init_flyte

    init_flyte()

    mo.md(
        """
        # Germline Short Variant Discovery

        Raw paired reads to a joint-genotyped cohort VCF, following GATK's
        best practices:

        1. **Prepare the reference**: FASTA index, sequence dictionary and
           aligner index.
        2. **Preprocess each sample** in parallel: align with BWA-MEM2,
           sort, mark duplicates.
        3. **Call each sample** with HaplotypeCaller into a GVCF.
        4. **Joint-genotype the cohort**: GenomicsDBImport, then
           GenotypeGVCFs.

        The workflow takes a build and sample IDs and finds its inputs in
        your asset graph, so pick them below rather than pointing at files.
        """
    )


@app.cell
def _():
    """Section 1 — the demo inputs button."""
    add_demo = mo.ui.run_button(label="Add the demo inputs")
    mo.vstack(
        [
            mo.md(
                f"""
                ## 1. Inputs From Your Graph

                The pickers below list every reference and every sample with
                paired reads that your graph holds. If it holds none yet, the
                `{DEMO_BUNDLE}` bundle is added for you: a slice of GRCh38
                around TP53 (build `GRCh38_TP53`) and sample NA12829's reads.
                Add it to a graph that already has data with the button.
                """
            ),
            add_demo,
        ]
    )
    return (add_demo,)


@app.cell
async def _(add_demo):
    """The graph's references and samples, adding the demo bundle if asked or empty."""
    from stargazer.bundles import fetch_bundle

    async def _inputs():
        """The graph's references and the R1 reads of its samples."""
        _refs = [
            a for a in await assemble(asset="reference") if isinstance(a, Reference)
        ]
        _reads = [a for a in await assemble(asset="r1") if isinstance(a, R1)]
        return _refs, _reads

    references, reads = await _inputs()
    if add_demo.value or not references or not reads:
        with mo.status.spinner(title=f"Adding the {DEMO_BUNDLE} bundle..."):
            await fetch_bundle(DEMO_BUNDLE)
        references, reads = await _inputs()
    return reads, references


@app.cell
def _(reads, references):
    """Section 1 — build, sample and cohort pickers."""
    _builds = sorted({r.build for r in references if r.build})
    _samples = sorted({r.sample_id for r in reads if r.sample_id})
    mo.stop(
        not _builds or not _samples,
        mo.md("Your graph has no reference with a build, or no sample reads."),
    )
    build_picker = mo.ui.dropdown(
        options=_builds,
        value="GRCh38_TP53" if "GRCh38_TP53" in _builds else _builds[0],
        label="Reference build",
    )
    sample_picker = mo.ui.multiselect(
        options=_samples,
        value=["NA12829"] if "NA12829" in _samples else _samples[:1],
        label="Samples",
    )
    cohort_input = mo.ui.text(value="cohort", label="Cohort ID")
    mo.hstack([build_picker, sample_picker, cohort_input], justify="start")
    return build_picker, cohort_input, sample_picker


@app.cell
def _():
    """Section 2 — where the pipeline runs."""
    from stargazer.notebooks.ui import dispatch, run_location

    where = run_location()
    mo.vstack(
        [
            mo.md(
                """
                ## 2. Where It Runs

                **On the cluster** runs every step in its own pod, samples
                in parallel. **On this machine** runs it in the notebook's
                own process, with the tools this machine has; joint calling
                needs an x86-64 machine, so on Apple silicon run on the
                cluster.
                """
            ),
            where,
        ]
    )
    return dispatch, where


@app.cell
def _(build_picker, cohort_input, sample_picker):
    """Section 3 — the run button."""
    run_button = mo.ui.run_button(label="Run the pipeline")
    mo.vstack(
        [
            mo.md(
                f"""
                ## 3. Run It

                Calls **{", ".join(sample_picker.value) or "no samples"}**
                against **{build_picker.value}** as cohort
                **{cohort_input.value}**. A run reuses the outputs of
                earlier runs with the same inputs.
                """
            ),
            run_button,
        ]
    )
    return (run_button,)


@app.cell
async def _(build_picker, cohort_input, dispatch, run_button, sample_picker, where):
    """Run the workflow where the switch says, once asked or started with `--run`."""
    mo.stop(not sample_picker.value, mo.md("Pick at least one sample."))
    mo.stop(
        not run_button.value and "run" not in mo.cli_args(),
        mo.md("Press **Run the pipeline** to start."),
    )
    with mo.status.spinner(
        title=f"Calling {len(sample_picker.value)} sample(s) on {build_picker.value}..."
    ):
        cohort = await dispatch(
            where,
            germline_short_variant_discovery,
            build=build_picker.value,
            sample_ids=list(sample_picker.value),
            cohort_id=cohort_input.value,
        )
    mo.md(
        f"Cohort **{cohort.sample_id}**: {len(cohort.source_samples)} sample(s) "
        f"on {cohort.build}, VCF `{cohort.cid}`."
    )
    return (cohort,)


@app.cell
async def _(cohort):
    """Section 4 — what the graph now holds for the cohort and its samples."""
    _kinds = ["alignment", "alignment_index", "variants", "variants_index"]
    _found = await assemble(
        sample_id=[cohort.sample_id, *cohort.source_samples], asset=_kinds
    )
    _rows = sorted(
        (
            {
                "asset": _kv["asset"],
                "sample": _kv.get("sample_id", ""),
                "made by": _kv.get("caller") or _kv.get("tool", ""),
                "links to": ", ".join(
                    f"{_k}={_kv[_k][:16]}…"
                    for _k in ("reference_cid", "companion_of")
                    if _kv.get(_k)
                ),
                "cid": _a.cid,
            }
            for _a in _found
            for _kv in [_a.to_keyvalues()]
        ),
        key=lambda _r: (_r["sample"], _kinds.index(_r["asset"])),
    )
    mo.vstack(
        [
            mo.md(
                """
                ## 4. See the Graph

                Everything your graph holds for this cohort and its samples:
                the alignments (and their indices) with the CID of the
                reference they were aligned to, each sample's GVCF, and the
                cohort VCF. Indices link to the file they index. Rows from
                earlier runs on the same samples show here too; the cohort
                VCF above is this run's.
                """
            ),
            mo.ui.table(_rows, selection=None, page_size=20),
        ]
    )


@app.cell
async def _(cohort):
    """Section 5 — a summary of the cohort VCF."""
    vcf_samples, vcf_records = read_vcf(Path(await cohort.fetch()))
    _counts = count_variants(vcf_records)
    _first = [
        {
            "chrom": _f[0],
            "pos": int(_f[1]),
            "ref": _f[3],
            "alt": _f[4],
            "qual": _f[5],
            "filter": _f[6],
            **{_s: _g.split(":")[0] for _s, _g in zip(vcf_samples, _f[9:])},
        }
        for _f in vcf_records[:20]
    ]
    mo.vstack(
        [
            mo.md(
                f"""
                ## 5. The Result

                **{len(vcf_records):,}** sites across
                **{len(vcf_samples)}** sample(s). Alleles by kind:
                """
            ),
            mo.ui.table(
                [
                    {"kind": _k, "alleles": _n}
                    for _k, _n in sorted(_counts.items(), key=lambda _kv: -_kv[1])
                ],
                selection=None,
            ),
            mo.md("The first sites, with each sample's genotype:"),
            mo.ui.table(_first, selection=None),
        ]
    )


if __name__ == "__main__":
    app.run()
