# /// script
# requires-python = ">=3.14"
# dependencies = [
#   "flyte",
#   "marimo",
#   "numpy",
#   "pandas",
#   "scikit-learn",
#   "stargazer-bio",
# ]
#
# [tool.pixi.workspace]
# # bioconda's packages build against conda-forge's, which come first.
# channels = ["conda-forge", "bioconda"]
# platforms = ["linux-64"]
#
# [tool.pixi.dependencies]
# python = "3.14.*"
# gatk4 = "4.6.*"
# bwa = "0.7.*"
# bwa-mem2 = "2.*"
# samtools = "1.*"
# bcftools = "1.*"
# fastp = "*"
#
# [tool.pixi.pypi-dependencies]
# stargazer-bio = { path = "/stargazer" }
#
# [tool.stargazer]
# owners = ["387300641116005877"]
# ///
"""
### Crop GMM: unsupervised VQSR training sets for crop variant discovery.

Bergmann, MacNish, Batley & Edwards' framework
(github.com/tbergmann92/Advanced-Crop-Variant-Discovery) as one Stargazer
workflow. It combines the repo's three Nextflow scripts (GATK variant calling,
SNP VQSR and INDEL VQSR) with the R step that links them, which fits a Gaussian
mixture model to a random subset of the cohort's site annotations and turns
its best clusters into VQSR training resources. Crops don't need curated
truth sets this way. An optional SNP-array known-sites VCF serves for BQSR and
as VQSR's `known` resource.

`crop_variant_discovery` takes a build and sample IDs; each parameter of the
repo's scripts is one of its arguments, with the B. napus values as defaults.
"""

import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")

with app.setup:
    import asyncio
    from pathlib import Path

    import flyte
    import marimo as mo

    from stargazer import (
        R1,
        R2,
        Alignment,
        KnownSites,
        Reference,
        Variants,
    )
    from stargazer.assets import (
        KnownSitesIndex,
        ReferenceIndex,
        VariantsIndex,
        VQSRModel,
        VQSRModelIndex,
        VQSRTranches,
    )
    from stargazer.config import STARGAZER_ENV_VARS, STARGAZER_SECRETS
    from stargazer.images import notebook_image
    from stargazer.notebooks import notebook_env
    from stargazer.tasks import (
        apply_bqsr,
        apply_vqsr,
        base_recalibrator,
        bwa_index,
        bwa_mem,
        bwa_mem2_index,
        bwa_mem2_mem,
        create_sequence_dictionary,
        joint_call_gvcfs,
        mark_duplicates,
        samtools_faidx,
        sort_sam,
    )
    from stargazer.tasks.gatk.env import gatk_env
    from stargazer.tasks.general.bwa_mem2 import BWA_MEM2_INDEX_EXTENSIONS
    from stargazer.utils import get_client, run, temp_dir

    main_name = "Crop GMM"

    main_description = "Unsupervised Clustering Enables Cohort-Specific Variant Recalibration in Crops"
    main_env = notebook_env.clone_with(
        name="notebook", resources=flyte.Resources(cpu=2, memory="4Gi")
    )

    # The pods this notebook's tasks run in, on an image built from the header
    # above (GATK, BWA, samtools, bcftools, fastp, scikit-learn). The SDK's
    # tasks it calls run on `gatk_env`. GATK's JVM sizes its heap from the
    # memory limit; raise it for VariantRecalibrator on large cohorts.
    crop_env = flyte.TaskEnvironment(
        name="crop-gmm",
        description="GATK, BWA, bcftools, fastp and scikit-learn for crop VQSR",
        image=notebook_image(__file__, "stargazer-crop-gmm"),
        resources=flyte.Resources(cpu=2, memory=("4Gi", "16Gi"), disk="50Gi"),
        env_vars=STARGAZER_ENV_VARS,
        secrets=STARGAZER_SECRETS,
        depends_on=[gatk_env],
    )

    BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")

    # The annotations the repo's R scripts cluster on, in every case study.
    GMM_FEATURES = {
        "SNP": ["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"],
        "INDEL": ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"],
    }

    # GATK's generic hard filters as the R scripts apply them: (annotation,
    # fails when "<" or ">" the threshold). Clusters are ranked by how many a
    # site breaks on average.
    HARD_FILTERS = {
        "SNP": [
            ("QD", "<", 2.0),
            ("FS", ">", 60.0),
            ("SOR", ">", 3.0),
            ("MQ", "<", 40.0),
            ("MQRankSum", "<", -12.5),
            ("ReadPosRankSum", "<", -8.0),
        ],
        "INDEL": [
            ("QD", "<", 2.0),
            ("FS", ">", 200.0),
            ("SOR", ">", 3.0),
            ("MQRankSum", "<", -12.5),
            ("ReadPosRankSum", "<", -20.0),
        ],
    }

    # Training tiers, best first, named as the R scripts name their clusters.
    TIER_NAMES = ["very_good", "good", "medium"]

    # The case studies' parameters: the repo's R clustering scripts (k, the
    # posterior cutoffs, the subset size) and VQSR scripts (priors,
    # annotations, filter level). B. oleracea has no VQSR script; its priors
    # follow B. napus.
    PRESETS = {
        "B. napus (case study 1)": dict(
            subset=1_000_000,
            snp_k=3, snp_priors="15, 12", snp_truth="true, false",
            snp_probs="0.8, 0.7",
            snp_ann=["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"],
            snp_level=99.5,
            indel=True, indel_k=4, indel_priors="15, 12, 10",
            indel_truth="true, true, false", indel_probs="0.7, 0.7, 0.7",
            indel_ann=["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"],
            indel_level=99.5,
        ),
        "B. oleracea (case study 2)": dict(
            subset=1_000_000,
            snp_k=3, snp_priors="15, 12", snp_truth="true, false",
            snp_probs="0.9, 0.9",
            snp_ann=["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"],
            snp_level=99.5,
            indel=True, indel_k=3, indel_priors="15, 12",
            indel_truth="true, false", indel_probs="0.9, 0.9",
            indel_ann=["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"],
            indel_level=99.5,
        ),
        "T. aestivum (case study 3)": dict(
            subset=3_000_000,
            snp_k=3, snp_priors="15, 12", snp_truth="true, false",
            snp_probs="0.9, 0.9",
            snp_ann=["QD", "FS", "SOR", "MQ"],
            snp_level=99.9,
            indel=False, indel_k=4, indel_priors="15, 12, 10",
            indel_truth="true, true, false", indel_probs="0.7, 0.7, 0.7",
            indel_ann=["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"],
            indel_level=99.5,
        ),
    }


@app.function
@crop_env.task
async def prepare_crop_reference(build: str, aligner: str = "bwa") -> Reference:
    """Find the build's reference and build the indices it lacks: .fai, .dict, aligner index."""
    if aligner not in ("bwa", "bwa-mem2"):
        raise ValueError(f"aligner must be 'bwa' or 'bwa-mem2', got {aligner!r}")
    ref = await Reference.find(build=build)
    ref_path = await ref.fetch()
    if not ref_path.with_name(f"{ref_path.name}.fai").exists():
        await samtools_faidx(ref)
    if not ref_path.with_suffix(".dict").exists():
        await create_sequence_dictionary(ref)
    extensions = BWA_INDEX_EXTENSIONS if aligner == "bwa" else BWA_MEM2_INDEX_EXTENSIONS
    if not all(ref_path.with_name(ref_path.name + e).exists() for e in extensions):
        await (bwa_index(ref) if aligner == "bwa" else bwa_mem2_index(ref))
    return ref


@app.function
@crop_env.task
async def fastp_trim(r1: R1, r2: R2) -> tuple[R1, R2]:
    """Trim a read pair with fastp: adapters detected for PE, base correction on.

    The trimmed pair is stored with `trimmed_by=fastp`, so a later run still
    finds the raw reads.
    """
    r1_path, r2_path = await asyncio.gather(r1.fetch(), r2.fetch())
    output_dir = temp_dir()
    out1 = output_dir / f"{r1.sample_id}_f1.clean.fq.gz"
    out2 = output_dir / f"{r1.sample_id}_r2.clean.fq.gz"
    await run(
        [
            "fastp",
            "-i", r1_path, "-I", r2_path,
            "-o", out1, "-O", out2,
            "--detect_adapter_for_pe",
            "--correction",
            "--thread", 8,
            "--json", output_dir / "fastp.json",
            "--html", output_dir / "fastp.html",
        ],
        cwd=output_dir,
    )
    client = get_client()
    trimmed_r1 = R1(sample_id=r1.sample_id, sequencing_platform=r1.sequencing_platform)
    trimmed_r1.keyvalues["trimmed_by"] = "fastp"
    await client.upload(trimmed_r1, out1)
    trimmed_r2 = R2(
        sample_id=r2.sample_id,
        mate_cid=trimmed_r1.cid,
        sequencing_platform=r2.sequencing_platform,
    )
    trimmed_r2.keyvalues["trimmed_by"] = "fastp"
    await client.upload(trimmed_r2, out2)
    return trimmed_r1, trimmed_r2


@app.function
@crop_env.task
async def preprocess_crop_sample(
    ref: Reference,
    sample_id: str,
    run_fastp: bool = True,
    aligner: str = "bwa",
    known_sites: list[KnownSites] = [],
) -> Alignment:
    """One sample's reads to an analysis-ready BAM.

    fastp (paired reads only), BWA-MEM or BWA-MEM2 with a full read group,
    coordinate sort, MarkDuplicates, then BQSR against `known_sites` when
    there are any.
    """
    r1s = [r for r in await R1.find_all(sample_id=sample_id) if "trimmed_by" not in r.keyvalues]
    if not r1s:
        raise ValueError(f"No raw R1 reads found for sample {sample_id!r}")
    r1 = r1s[0]
    r2s = [r for r in await R2.find_all(sample_id=sample_id) if "trimmed_by" not in r.keyvalues]
    r2 = next((r for r in r2s if r.mate_cid == r1.cid), r2s[0] if r2s else None)

    if run_fastp and r2 is not None:
        r1, r2 = await fastp_trim(r1=r1, r2=r2)

    read_group = {"PL": "ILLUMINA", "LB": sample_id, "PU": sample_id}
    align = bwa_mem if aligner == "bwa" else bwa_mem2_mem
    alignment = await align(ref=ref, r1=r1, r2=r2, read_group=read_group)
    alignment = await sort_sam(alignment=alignment, sort_order="coordinate")
    alignment = await mark_duplicates(alignment=alignment)

    if known_sites:
        report = await base_recalibrator(alignment=alignment, ref=ref, known_sites=known_sites)
        alignment = await apply_bqsr(alignment=alignment, ref=ref, bqsr_report=report)
    return alignment


@app.function
@crop_env.task
async def haplotype_caller_shard(
    alignment: Alignment,
    ref: Reference,
    intervals: list[str],
    shard: str,
) -> Variants:
    """HaplotypeCaller in GVCF mode over one shard of the reference's contigs."""
    bam_path, ref_path = await asyncio.gather(alignment.fetch(), ref.fetch())
    output_dir = temp_dir()
    interval_list = output_dir / f"{shard}.list"
    interval_list.write_text("\n".join(intervals) + "\n")
    output_gvcf = output_dir / f"{alignment.sample_id}.{shard}.g.vcf.gz"
    await run(
        [
            "gatk", "HaplotypeCaller",
            "-R", ref_path,
            "-I", bam_path,
            "-L", interval_list,
            "-O", output_gvcf,
            "-ERC", "GVCF",
            "--native-pair-hmm-threads", 4,
        ],
        cwd=output_dir,
    )
    gvcf = Variants(
        sample_id=alignment.sample_id,
        caller="haplotype_caller",
        variant_type="gvcf",
        build=ref.build,
        source_samples=[alignment.sample_id],
        source_cids=[alignment.cid],
    )
    gvcf.keyvalues["shard"] = shard
    await get_client().upload(gvcf, output_gvcf)
    await VariantsIndex.store(
        output_dir / f"{output_gvcf.name}.tbi",
        sample_id=alignment.sample_id,
        companion_of=gvcf.cid,
    )
    return gvcf


@app.function
@crop_env.task
async def concat_vcfs(parts: list[Variants], ref: Reference, cohort_id: str) -> Variants:
    """Concatenate the shards' joint-called VCFs, in reference order, into one bgzipped VCF."""
    part_paths = [await p.fetch() for p in parts]
    output_dir = temp_dir()
    output_vcf = output_dir / f"{cohort_id}.merged.raw.vcf.gz"
    await run(["bcftools", "concat", "-Oz", "-o", output_vcf, *part_paths], cwd=output_dir)
    await run(["bcftools", "index", "-t", output_vcf], cwd=output_dir)
    samples = parts[0].source_samples or []
    vcf = await Variants.store(
        output_vcf,
        sample_id=cohort_id,
        caller="crop_joint_call",
        variant_type="vcf",
        build=ref.build,
        source_samples=samples,
        source_cids=[p.cid for p in parts],
    )
    await VariantsIndex.store(
        output_dir / f"{output_vcf.name}.tbi", sample_id=cohort_id, companion_of=vcf.cid
    )
    return vcf


@app.function
async def store_vcf_resource(
    vcf_path: Path, build: str, resource_name: str, known: bool, training: bool,
    truth: bool, prior: float, mode: str, extras: dict[str, str],
) -> KnownSites:
    """Index a bgzipped VCF and store it as a VQSR resource, the index as its companion."""
    await run(["bcftools", "index", "-f", "-t", vcf_path])
    resource = KnownSites(
        build=build,
        resource_name=resource_name,
        known=str(known).lower(),
        training=str(training).lower(),
        truth=str(truth).lower(),
        prior=f"{prior:g}",
        vqsr_mode=mode,
    )
    resource.keyvalues.update(extras)
    await get_client().upload(resource, vcf_path)
    await KnownSitesIndex.store(
        vcf_path.with_name(f"{vcf_path.name}.tbi"), companion_of=resource.cid
    )
    return resource


@app.function
@crop_env.task
async def gmm_training_sets(
    vcf: Variants,
    mode: str = "SNP",
    n_components: int = 3,
    priors: list[float] = [15.0, 12.0],
    truth: list[bool] = [True, False],
    min_probs: list[float] = [0.8, 0.7],
    subset_size: int = 1_000_000,
    seed: int = 123,
) -> list[KnownSites]:
    """Cluster a cohort's site annotations with a GMM and store its best clusters as VQSR training sets.

    The repo's R step, unattended: VariantsToTable, a random subset of
    `subset_size` sites of this type, z-scored annotations, a
    diagonal-covariance GMM with `n_components` components (ClusterR's
    `GMM`). The components are ranked by the mean number of GATK hard filters
    their sites break, ties broken by higher QD, in place of the scripts' hand
    labels. The best `len(priors)` become training tiers (very_good, good,
    medium); the worst always stays out. A tier keeps the sites assigned to it
    with a posterior of at least its `min_probs`, and is stored as a
    `KnownSites` with its `prior` and `truth` flag.
    """
    import numpy as np
    import pandas as pd
    from sklearn.mixture import GaussianMixture

    if mode not in GMM_FEATURES:
        raise ValueError(f"mode must be 'SNP' or 'INDEL', got {mode!r}")
    n_tiers = len(priors)
    if not (len(truth) == len(min_probs) == n_tiers):
        raise ValueError("priors, truth and min_probs need one entry per training tier")
    if not 0 < n_tiers < n_components:
        raise ValueError(f"Need 1 to {n_components - 1} tiers for k={n_components}, got {n_tiers}")
    names = TIER_NAMES[-n_tiers:] if n_tiers <= len(TIER_NAMES) else [f"tier{i + 1}" for i in range(n_tiers)]

    features = GMM_FEATURES[mode]
    vcf_path = await vcf.fetch()
    output_dir = temp_dir()
    table = output_dir / "annotations.tsv"
    fields = [a for f in ["CHROM", "POS", "TYPE", *features] for a in ("-F", f)]
    await run(["gatk", "VariantsToTable", "-V", vcf_path, *fields, "-O", table], cwd=output_dir)

    df = pd.read_csv(table, sep="\t", dtype={"CHROM": str})
    df = df[df["TYPE"] == mode]
    if len(df) > subset_size:
        df = df.sample(n=subset_size, random_state=seed)
    df = df.dropna(subset=features)
    if len(df) < n_components * 10:
        raise ValueError(f"Only {len(df)} {mode} sites with every annotation; too few to cluster")

    values = df[features].to_numpy(dtype=float)
    scaled = (values - values.mean(axis=0)) / values.std(axis=0, ddof=1)
    gmm = GaussianMixture(
        n_components=n_components, covariance_type="diag", max_iter=100,
        init_params="kmeans", random_state=seed,
    ).fit(scaled)
    proba = gmm.predict_proba(scaled)
    labels = proba.argmax(axis=1)

    violations = np.zeros(len(df))
    for annotation, op, threshold in HARD_FILTERS[mode]:
        column = df[annotation].to_numpy(dtype=float)
        violations += column < threshold if op == "<" else column > threshold
    qd = df["QD"].to_numpy(dtype=float)
    ranked = sorted(
        range(n_components),
        key=lambda c: (violations[labels == c].mean() if (labels == c).any() else np.inf,
                       -qd[labels == c].mean() if (labels == c).any() else 0.0),
    )

    select_type = "snps" if mode == "SNP" else "indels"
    resources = []
    for tier, (name, prior, is_truth, min_prob) in enumerate(zip(names, priors, truth, min_probs)):
        component = ranked[tier]
        keep = (labels == component) & (proba[:, component] >= min_prob)
        sites = df.loc[keep, ["CHROM", "POS"]]
        positions = output_dir / f"{name}.positions.txt"
        sites.to_csv(positions, sep="\t", header=False, index=False)
        tier_vcf = output_dir / f"{vcf.sample_id}.{mode}.GMM_{name}.vcf.gz"
        await run(
            ["bcftools", "view", "-T", positions, "-v", select_type, "-G",
             "-Oz", "-o", tier_vcf, vcf_path],
            cwd=output_dir,
        )
        resources.append(
            await store_vcf_resource(
                tier_vcf, vcf.build, f"GMM_{name}", known=False, training=True,
                truth=is_truth, prior=prior, mode=mode,
                extras={
                    "cohort_id": vcf.sample_id,
                    "variants_cid": vcf.cid,
                    "gmm_k": str(n_components),
                    "gmm_component": str(component),
                    "min_prob": f"{min_prob:g}",
                    "n_sites": str(int(keep.sum())),
                    "clustered_sites": str(len(df)),
                    "mean_hf_violations": f"{violations[labels == component].mean():.3f}",
                    "median_qd": f"{np.median(qd[labels == component]):.2f}",
                },
            )
        )
    return resources


@app.function
@crop_env.task
async def array_known_resource(vcf: Variants, array_sites: KnownSites, prior: float = 2.0) -> KnownSites:
    """The cohort's SNPs at SNP-array sites, as VQSR's `known` resource (the repo's TRUTH_SET VCF)."""
    vcf_path, array_path = await asyncio.gather(vcf.fetch(), array_sites.fetch())
    output_vcf = temp_dir() / f"{vcf.sample_id}.TRUTH_SET.vcf.gz"
    await run(
        ["bcftools", "view", "-T", array_path, "-v", "snps", "-G", "-Oz", "-o", output_vcf, vcf_path]
    )
    return await store_vcf_resource(
        output_vcf, vcf.build, "truth_set", known=True, training=False, truth=False,
        prior=prior, mode="SNP",
        extras={"cohort_id": vcf.sample_id, "variants_cid": vcf.cid, "array_cid": array_sites.cid},
    )


@app.function
@crop_env.task
async def crop_variant_recalibrator(
    vcf: Variants,
    ref: Reference,
    resources: list[KnownSites],
    mode: str = "SNP",
    annotations: list[str] = ["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"],
    tranches: list[float] = [100.0, 99.99, 99.9, 99.7, 99.5, 99.0, 95.0, 90.0],
    trust_all_polymorphic: bool = False,
    max_gaussians: int = 0,
) -> VQSRModel:
    """GATK VariantRecalibrator with the repo's annotations, tranches and options.

    The SDK's `variant_recalibrator` fixes its annotations; this one takes
    them, so INDELs can use MQRankSum and ReadPosRankSum as the paper does.
    Stores the model with its index and tranches as companions, so the SDK's
    `apply_vqsr` takes it.
    """
    vcf_path, ref_path = await asyncio.gather(vcf.fetch(), ref.fetch())
    resource_paths = [await r.fetch() for r in resources]
    output_dir = temp_dir()
    stem = f"{vcf.sample_id or 'cohort'}.{mode}"
    recal = output_dir / f"{stem}.recal"
    tranches_file = output_dir / f"{stem}.tranches"
    cmd = [
        "gatk", "VariantRecalibrator",
        "-R", ref_path, "-V", vcf_path, "-mode", mode,
        "-O", recal, "--tranches-file", tranches_file,
    ]
    if trust_all_polymorphic:
        cmd.append("--trust-all-polymorphic")
    if max_gaussians:
        cmd.extend(["--max-gaussians", max_gaussians])
    for ann in annotations:
        cmd.extend(["-an", ann])
    for r, path in zip(resources, resource_paths, strict=True):
        cmd.extend([
            f"--resource:{r.resource_name},known={r.known},training={r.training},"
            f"truth={r.truth},prior={r.prior}",
            path,
        ])
    for level in tranches:
        cmd.extend(["-tranche", f"{level:g}"])
    await run(cmd, cwd=output_dir)

    model = await VQSRModel.store(
        recal, sample_id=vcf.sample_id, mode=mode, build=ref.build, variants_cid=vcf.cid
    )
    await VQSRModelIndex.store(
        output_dir / f"{recal.name}.idx", sample_id=vcf.sample_id, companion_of=model.cid
    )
    await VQSRTranches.store(tranches_file, sample_id=vcf.sample_id, companion_of=model.cid)
    return model


@app.function
def contig_shards(fai_text: str, scatter_count: int) -> list[list[str]]:
    """Split a .fai's contigs, in order, into at most `scatter_count` runs of about equal length."""
    contigs = [(f[0], int(f[1])) for line in fai_text.splitlines() if line.strip() for f in [line.split("\t")]]
    target = sum(length for _, length in contigs) / max(scatter_count, 1)
    shards, current, size = [], [], 0
    for name, length in contigs:
        if current and size + length / 2 > target and len(shards) < scatter_count - 1:
            shards.append(current)
            current, size = [], 0
        current.append(name)
        size += length
    shards.append(current)
    return shards


@app.function
@crop_env.task(cache="disable")
async def crop_variant_discovery(
    build: str,
    sample_ids: list[str],
    cohort_id: str = "cohort",
    call_variants: bool = True,
    run_fastp: bool = True,
    aligner: str = "bwa",
    known_sites_name: str = "snp_array",
    run_bqsr: bool = True,
    scatter_count: int = 24,
    gmm_subset_size: int = 1_000_000,
    gmm_seed: int = 123,
    snp_components: int = 3,
    snp_priors: list[float] = [15.0, 12.0],
    snp_truth: list[bool] = [True, False],
    snp_min_probs: list[float] = [0.8, 0.7],
    snp_annotations: list[str] = ["QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum"],
    snp_filter_level: float = 99.5,
    array_prior: float = 2.0,
    run_indel_vqsr: bool = True,
    indel_components: int = 4,
    indel_priors: list[float] = [15.0, 12.0, 10.0],
    indel_truth: list[bool] = [True, True, False],
    indel_min_probs: list[float] = [0.7, 0.7, 0.7],
    indel_annotations: list[str] = ["QD", "FS", "SOR", "MQRankSum", "ReadPosRankSum"],
    indel_filter_level: float = 99.5,
    tranches: list[float] = [100.0, 99.99, 99.9, 99.7, 99.5, 99.0, 95.0, 90.0],
    trust_all_polymorphic: bool = True,
    vqsr_max_gaussians: int = 0,
) -> Variants:
    """Reads to a GMM-trained, VQSR-filtered cohort VCF for a crop.

    1. Calling (the repo's GATK script), unless `call_variants` is off, which
       reuses the newest raw VCF stored for `cohort_id`: fastp, BWA, sort,
       MarkDuplicates, BQSR on the `known_sites_name` sites, HaplotypeCaller
       and GenomicsDBImport + GenotypeGVCFs on `scatter_count` contig shards,
       then bcftools concat.
    2. GMM training sets (the R scripts), SNPs and INDELs in parallel, plus
       the array sites in the cohort as SNP VQSR's `known` resource.
    3. SNP VQSR, then INDEL VQSR on its output (the two VQSR scripts).

    Returns the filtered VCF; each tier, the models and the intermediate VCFs
    are stored too.
    """
    ref = await prepare_crop_reference(build=build, aligner=aligner)
    known = (
        await KnownSites.find(build=build, resource_name=known_sites_name, required=False)
        if known_sites_name
        else None
    )

    if call_variants:
        bqsr_sites = [known] if (known and run_bqsr) else []
        alignments = await asyncio.gather(
            *[
                preprocess_crop_sample(
                    ref=ref, sample_id=sid, run_fastp=run_fastp,
                    aligner=aligner, known_sites=bqsr_sites,
                )
                for sid in sample_ids
            ]
        )
        fai = await ReferenceIndex.find(companion_of=ref.cid)
        shards = contig_shards((await fai.fetch()).read_text(), scatter_count)
        names = [f"shard{i:03d}" for i in range(len(shards))]
        gvcfs = await asyncio.gather(
            *[
                haplotype_caller_shard(alignment=aln, ref=ref, intervals=contigs, shard=name)
                for name, contigs in zip(names, shards)
                for aln in alignments
            ]
        )
        n = len(alignments)
        parts = await asyncio.gather(
            *[
                joint_call_gvcfs(
                    gvcfs=list(gvcfs[i * n : (i + 1) * n]), ref=ref,
                    intervals=contigs, cohort_id=f"{cohort_id}.{name}",
                )
                for i, (name, contigs) in enumerate(zip(names, shards))
            ]
        )
        raw_vcf = await concat_vcfs(parts=list(parts), ref=ref, cohort_id=cohort_id)
    else:
        found = [
            v for v in await Variants.find_all(sample_id=cohort_id, build=build)
            if v.variant_type == "vcf" and not v.vqsr_mode
        ]
        if not found:
            raise ValueError(f"No raw cohort VCF stored for cohort {cohort_id!r} on {build!r}")
        raw_vcf = found[0]

    jobs = [
        gmm_training_sets(
            vcf=raw_vcf, mode="SNP", n_components=snp_components, priors=snp_priors,
            truth=snp_truth, min_probs=snp_min_probs, subset_size=gmm_subset_size, seed=gmm_seed,
        )
    ]
    if run_indel_vqsr:
        jobs.append(
            gmm_training_sets(
                vcf=raw_vcf, mode="INDEL", n_components=indel_components,
                priors=indel_priors, truth=indel_truth, min_probs=indel_min_probs,
                subset_size=gmm_subset_size, seed=gmm_seed,
            )
        )
    if known:
        jobs.append(array_known_resource(vcf=raw_vcf, array_sites=known, prior=array_prior))
    results = await asyncio.gather(*jobs)
    snp_resources = list(results[0]) + ([results[-1]] if known else [])

    snp_model = await crop_variant_recalibrator(
        vcf=raw_vcf, ref=ref, resources=snp_resources, mode="SNP",
        annotations=snp_annotations, tranches=tranches,
        trust_all_polymorphic=trust_all_polymorphic, max_gaussians=vqsr_max_gaussians,
    )
    filtered = await apply_vqsr(
        vcf=raw_vcf, ref=ref, vqsr_model=snp_model,
        truth_sensitivity_filter_level=snp_filter_level,
    )
    if not run_indel_vqsr:
        return filtered

    indel_model = await crop_variant_recalibrator(
        vcf=filtered, ref=ref, resources=list(results[1]), mode="INDEL",
        annotations=indel_annotations, tranches=tranches,
        max_gaussians=vqsr_max_gaussians,
    )
    return await apply_vqsr(
        vcf=filtered, ref=ref, vqsr_model=indel_model,
        truth_sensitivity_filter_level=indel_filter_level,
    )


@app.function
def parse_list(text: str, kind: type) -> list:
    """A comma-separated text field as a list of `kind`; "true"/"yes"/"1" read as True for bool."""
    items = [t.strip() for t in text.split(",") if t.strip()]
    if kind is bool:
        return [t.lower() in ("true", "yes", "1") for t in items]
    return [kind(t) for t in items]


@app.cell
def _():
    """Title."""
    mo.md(
        """
        # Crop GMM: Variant Recalibration Without Truth Sets

        GATK's VQSR needs curated truth and training resources, which human
        genomics has (HapMap, Omni, 1000G) and most crops don't. Bergmann et
        al. build them from the cohort: a Gaussian mixture model clusters a
        random subset of the raw calls by their site annotations (QD, FS,
        SOR, MQ, MQRankSum, ReadPosRankSum), and its high- and
        intermediate-quality clusters become VQSR's training sets. An SNP
        array, where one exists, gives BQSR its known sites and VQSR a
        `known` resource.

        `crop_variant_discovery` runs the repo's three scripts and the
        clustering step between them as one workflow:

        1. **Call**: fastp → BWA → sort → MarkDuplicates → BQSR →
           HaplotypeCaller and joint genotyping on contig shards →
           concatenate.
        2. **Cluster**: SNPs and INDELs separately. Components are ranked by
           how many GATK hard filters their sites break; the best become
           *very good / good / medium* training tiers, each tier keeping
           sites above a posterior cutoff.
        3. **Recalibrate**: SNP VQSR on the raw VCF, then INDEL VQSR on its
           output.
        """
    )
    return


@app.cell
async def _():
    """The graph's references, raw reads and known-sites VCFs."""
    references = [a for a in await Reference.find_all() if a.build]
    reads = [r for r in await R1.find_all() if "trimmed_by" not in r.keyvalues]
    known_sites = await KnownSites.find_all()
    return known_sites, reads, references


@app.cell
def _(reads, references):
    """Section 1 — build, samples, cohort."""
    _builds = sorted({r.build for r in references})
    _samples = sorted({r.sample_id for r in reads if r.sample_id})
    mo.stop(
        not _builds,
        mo.md("Your graph holds no reference with a `build`. Store one as a `Reference` asset first."),
    )
    build_picker = mo.ui.dropdown(options=_builds, value=_builds[0], label="Reference build")
    sample_picker = mo.ui.multiselect(options=_samples, label="Samples")
    cohort_input = mo.ui.text(value="cohort", label="Cohort ID")
    mo.vstack(
        [
            mo.md("## 1. Inputs"),
            mo.hstack([build_picker, cohort_input], justify="start"),
            sample_picker,
        ]
    )
    return build_picker, cohort_input, sample_picker


@app.cell
def _(build_picker, known_sites):
    """Section 2 — calling options."""
    # The GMM tiers and array resources earlier runs stored are KnownSites too.
    _names = sorted(
        {
            k.resource_name
            for k in known_sites
            if k.build == build_picker.value
            and k.resource_name
            and k.resource_name != "truth_set"
            and not k.resource_name.startswith("GMM_")
        }
    )
    known_picker = mo.ui.dropdown(
        options=["(none)", *_names],
        value="snp_array" if "snp_array" in _names else "(none)",
        label="SNP-array known sites",
    )
    call_toggle = mo.ui.checkbox(value=True, label="Call variants (off: reuse the cohort's stored raw VCF)")
    fastp_toggle = mo.ui.checkbox(value=True, label="Trim with fastp")
    bqsr_toggle = mo.ui.checkbox(value=True, label="BQSR on the known sites")
    aligner_picker = mo.ui.dropdown(options=["bwa", "bwa-mem2"], value="bwa", label="Aligner")
    scatter_input = mo.ui.number(start=1, stop=500, value=24, label="Contig shards")
    mo.vstack(
        [
            mo.md(
                f"""
                ## 2. Calling

                Known sites are a `KnownSites` asset on **{build_picker.value}**,
                found by its `resource_name`, such as the repo's
                `Resources/known_sites_*.vcf.gz`, the 60k Brassica or 660k wheat
                array mapped to the reference. Without one, BQSR is skipped and
                SNP VQSR trains on the GMM tiers alone. BWA's index needs far
                less memory than BWA-MEM2's on large genomes like wheat.
                """
            ),
            mo.hstack([call_toggle, fastp_toggle, bqsr_toggle], justify="start"),
            mo.hstack([known_picker, aligner_picker, scatter_input], justify="start"),
        ]
    )
    return (
        aligner_picker,
        bqsr_toggle,
        call_toggle,
        fastp_toggle,
        known_picker,
        scatter_input,
    )


@app.cell
def _():
    """Section 3 — the case-study preset."""
    preset_picker = mo.ui.dropdown(
        options=list(PRESETS), value="B. napus (case study 1)", label="Start from"
    )
    mo.vstack(
        [
            mo.md(
                """
                ## 3. GMM and VQSR

                Each case study's parameters from the repo: the number of
                Gaussians (k), each tier's posterior cutoff and VQSR prior,
                the annotations VQSR uses and the truth-sensitivity filter
                level. Pick the closest one and edit it.
                """
            ),
            preset_picker,
        ]
    )
    return (preset_picker,)


@app.cell
def _(preset_picker):
    """Section 3 — the GMM and VQSR parameters, from the preset."""
    _p = PRESETS[preset_picker.value]
    _snp_options = GMM_FEATURES["SNP"]
    gmm_form = mo.ui.dictionary(
        {
            "subset": mo.ui.number(start=1000, stop=50_000_000, step=1000, value=_p["subset"], label="Sites clustered"),
            "seed": mo.ui.number(start=0, stop=2**31 - 1, value=123, label="Seed"),
            "snp_k": mo.ui.number(start=2, stop=10, value=_p["snp_k"], label="SNP k"),
            "snp_priors": mo.ui.text(value=_p["snp_priors"], label="SNP tier priors"),
            "snp_truth": mo.ui.text(value=_p["snp_truth"], label="SNP tier truth"),
            "snp_probs": mo.ui.text(value=_p["snp_probs"], label="SNP tier posterior ≥"),
            "snp_ann": mo.ui.multiselect(options=_snp_options, value=_p["snp_ann"], label="SNP VQSR annotations"),
            "snp_level": mo.ui.number(start=90, stop=100, step=0.01, value=_p["snp_level"], label="SNP filter level"),
            "array_prior": mo.ui.number(start=0, stop=20, step=0.5, value=2.0, label="Array known prior"),
            "indel": mo.ui.checkbox(value=_p["indel"], label="INDEL VQSR"),
            "indel_k": mo.ui.number(start=2, stop=10, value=_p["indel_k"], label="INDEL k"),
            "indel_priors": mo.ui.text(value=_p["indel_priors"], label="INDEL tier priors"),
            "indel_truth": mo.ui.text(value=_p["indel_truth"], label="INDEL tier truth"),
            "indel_probs": mo.ui.text(value=_p["indel_probs"], label="INDEL tier posterior ≥"),
            "indel_ann": mo.ui.multiselect(options=_snp_options, value=_p["indel_ann"], label="INDEL VQSR annotations"),
            "indel_level": mo.ui.number(start=90, stop=100, step=0.01, value=_p["indel_level"], label="INDEL filter level"),
            "tranches": mo.ui.text(value="100.0, 99.99, 99.9, 99.7, 99.5, 99.0, 95.0, 90.0", label="Tranches"),
            "max_gaussians": mo.ui.number(start=0, stop=10, value=0, label="VQSR max Gaussians (0 = GATK's)"),
        }
    )
    _f = gmm_form
    mo.vstack(
        [
            mo.hstack([_f["subset"], _f["seed"]], justify="start"),
            mo.md("**SNPs**: tiers are best first, one entry each, separated by commas."),
            mo.hstack([_f["snp_k"], _f["snp_priors"], _f["snp_truth"], _f["snp_probs"]], justify="start", wrap=True),
            mo.hstack([_f["snp_ann"], _f["snp_level"], _f["array_prior"]], justify="start", wrap=True),
            mo.md("**INDELs**"),
            _f["indel"],
            mo.hstack([_f["indel_k"], _f["indel_priors"], _f["indel_truth"], _f["indel_probs"]], justify="start", wrap=True),
            mo.hstack([_f["indel_ann"], _f["indel_level"]], justify="start", wrap=True),
            mo.md("**VQSR**"),
            mo.hstack([_f["tranches"], _f["max_gaussians"]], justify="start", wrap=True),
        ]
    )
    return (gmm_form,)


@app.cell
def _():
    """Section 4 — where it runs."""
    from stargazer.notebooks.ui import dispatch, run_location

    where = run_location()
    run_button = mo.ui.run_button(label="Run crop_variant_discovery")
    mo.vstack(
        [
            mo.md(
                """
                ## 4. Run It

                **On the cluster** gives each sample, shard and step its own
                pod. **On this machine** runs everything in the notebook's
                process with the tools from its header; that only suits a
                small test region.
                """
            ),
            where,
            run_button,
        ]
    )
    return dispatch, run_button, where


@app.cell
def _(
    aligner_picker,
    bqsr_toggle,
    build_picker,
    call_toggle,
    cohort_input,
    fastp_toggle,
    gmm_form,
    known_picker,
    sample_picker,
    scatter_input,
):
    """The workflow's inputs, from every control above."""
    _g = gmm_form.value
    workflow_inputs = dict(
        build=build_picker.value,
        sample_ids=list(sample_picker.value),
        cohort_id=cohort_input.value,
        call_variants=call_toggle.value,
        run_fastp=fastp_toggle.value,
        aligner=aligner_picker.value,
        known_sites_name="" if known_picker.value == "(none)" else known_picker.value,
        run_bqsr=bqsr_toggle.value,
        scatter_count=int(scatter_input.value),
        gmm_subset_size=int(_g["subset"]),
        gmm_seed=int(_g["seed"]),
        snp_components=int(_g["snp_k"]),
        snp_priors=parse_list(_g["snp_priors"], float),
        snp_truth=parse_list(_g["snp_truth"], bool),
        snp_min_probs=parse_list(_g["snp_probs"], float),
        snp_annotations=list(_g["snp_ann"]),
        snp_filter_level=float(_g["snp_level"]),
        array_prior=float(_g["array_prior"]),
        run_indel_vqsr=_g["indel"],
        indel_components=int(_g["indel_k"]),
        indel_priors=parse_list(_g["indel_priors"], float),
        indel_truth=parse_list(_g["indel_truth"], bool),
        indel_min_probs=parse_list(_g["indel_probs"], float),
        indel_annotations=list(_g["indel_ann"]),
        indel_filter_level=float(_g["indel_level"]),
        tranches=parse_list(_g["tranches"], float),
        vqsr_max_gaussians=int(_g["max_gaussians"]),
    )
    mo.accordion({"Inputs the run will get": mo.json(workflow_inputs)})
    return (workflow_inputs,)


@app.cell
async def _(dispatch, run_button, where, workflow_inputs):
    """Run the workflow where the switch says, once asked."""
    mo.stop(
        workflow_inputs["call_variants"] and not workflow_inputs["sample_ids"],
        mo.md("Pick at least one sample, or turn off calling to reuse a stored cohort VCF."),
    )
    mo.stop(not run_button.value, mo.md("Press **Run crop_variant_discovery** to start."))
    with mo.status.spinner(title=f"Running cohort {workflow_inputs['cohort_id']}..."):
        final_vcf = await dispatch(where, crop_variant_discovery, **workflow_inputs)
    mo.md(f"Filtered VCF for **{final_vcf.sample_id}** on {final_vcf.build}: `{final_vcf.cid}`")
    return (final_vcf,)


@app.cell
async def _(final_vcf):
    """Section 5 — the training resources the GMM built for this cohort."""
    _resources = [
        k for k in await KnownSites.find_all(build=final_vcf.build)
        if k.keyvalues.get("cohort_id") == final_vcf.sample_id
    ]
    _rows = [
        {
            "resource": k.resource_name,
            "mode": k.vqsr_mode,
            "truth": k.truth,
            "known": k.known,
            "prior": k.prior,
            "sites": k.keyvalues.get("n_sites", ""),
            "of clustered": k.keyvalues.get("clustered_sites", ""),
            "posterior ≥": k.keyvalues.get("min_prob", ""),
            "mean HF violations": k.keyvalues.get("mean_hf_violations", ""),
            "median QD": k.keyvalues.get("median_qd", ""),
            "raw VCF": k.keyvalues.get("variants_cid", "")[:16],
            "cid": k.cid,
        }
        for k in _resources
    ]
    mo.vstack(
        [
            mo.md(
                """
                ## 5. Training Resources

                Each GMM tier and the array resource stored for this cohort,
                newest first. Earlier runs' rows show too; the `raw VCF`
                column tells them apart.
                """
            ),
            mo.ui.table(_rows, selection=None),
        ]
    )
    return


@app.cell
async def _(final_vcf):
    """Section 6 — PASS and filtered sites by type in the recalibrated VCF."""
    import gzip
    from collections import Counter

    _path = Path(await final_vcf.fetch())
    _opener = gzip.open if _path.suffix == ".gz" else open
    _counts = Counter()
    with _opener(_path, "rt") as _f:
        for _line in _f:
            if _line.startswith("#"):
                continue
            _fields = _line.split("\t", 8)
            _ref, _alts = _fields[3], _fields[4].split(",")
            _kind = "SNP" if all(len(a) == len(_ref) == 1 for a in _alts) else "INDEL/other"
            _counts[(_kind, _fields[6])] += 1
    mo.vstack(
        [
            mo.md("## 6. The Result"),
            mo.ui.table(
                [
                    {"type": _k, "filter": _flt, "sites": _n}
                    for (_k, _flt), _n in sorted(_counts.items())
                ],
                selection=None,
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()
