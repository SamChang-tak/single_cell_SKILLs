#!/usr/bin/env python3
"""Summarize patient-level spatial niche prevalence and condition differences."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path, help="H5AD containing spot-level niche labels")
    p.add_argument("--patient-key", default="patient_id")
    p.add_argument("--condition-key", default="condition")
    p.add_argument("--niche-key", default="leiden")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    a = args()
    if not a.input.is_file():
        raise FileNotFoundError(a.input)
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(config, indent=2))
        return

    import anndata as ad
    import matplotlib.pyplot as plt
    import pandas as pd
    import scanpy as sc
    from scipy.stats import mannwhitneyu

    adata = sc.read_h5ad(a.input, backed="r")
    keys = [a.patient_key, a.condition_key, a.niche_key]
    missing = [key for key in keys if key not in adata.obs]
    if missing:
        raise KeyError(f"Missing obs columns: {missing}")
    obs = adata.obs[keys].copy()
    adata.file.close()
    if obs.isna().any().any():
        raise ValueError(f"Missing metadata values: {obs.isna().sum().to_dict()}")
    patient_conditions = obs.groupby(a.patient_key, observed=True)[a.condition_key].nunique()
    if (patient_conditions != 1).any():
        raise ValueError("Each patient must map to exactly one condition")

    counts = obs.groupby(keys, observed=True).size().rename("n_spots").reset_index()
    totals = obs.groupby(a.patient_key, observed=True).size().rename("patient_total_spots")
    counts = counts.join(totals, on=a.patient_key)
    counts["proportion"] = counts["n_spots"] / counts["patient_total_spots"]
    all_pairs = pd.MultiIndex.from_product(
        [obs[a.patient_key].unique(), obs[a.niche_key].unique()], names=[a.patient_key, a.niche_key]
    ).to_frame(index=False)
    condition_map = obs.drop_duplicates(a.patient_key).set_index(a.patient_key)[a.condition_key]
    prevalence = all_pairs.merge(counts, on=[a.patient_key, a.niche_key], how="left")
    prevalence[a.condition_key] = prevalence[a.condition_key].fillna(prevalence[a.patient_key].map(condition_map))
    prevalence["n_spots"] = prevalence["n_spots"].fillna(0).astype(int)
    prevalence["patient_total_spots"] = prevalence["patient_total_spots"].fillna(prevalence[a.patient_key].map(totals))
    prevalence["proportion"] = prevalence["proportion"].fillna(0.0)

    a.output_dir.mkdir(parents=True, exist_ok=True)
    prevalence.to_csv(a.output_dir / "patient_niche_prevalence.tsv", sep="\t", index=False)
    conditions = list(prevalence[a.condition_key].unique())
    tests = []
    if len(conditions) == 2:
        for niche, frame in prevalence.groupby(a.niche_key, observed=True):
            x = frame.loc[frame[a.condition_key].eq(conditions[0]), "proportion"]
            y = frame.loc[frame[a.condition_key].eq(conditions[1]), "proportion"]
            stat, pvalue = mannwhitneyu(x, y, alternative="two-sided") if len(x) >= 2 and len(y) >= 2 else (float("nan"), float("nan"))
            tests.append({"niche": niche, "condition_1": conditions[0], "condition_2": conditions[1], "n_1": len(x), "n_2": len(y), "median_difference_2_minus_1": y.median() - x.median(), "mannwhitney_u": stat, "pvalue_unadjusted": pvalue})
    pd.DataFrame(tests).to_csv(a.output_dir / "niche_condition_tests.tsv", sep="\t", index=False)

    niches = list(prevalence[a.niche_key].astype(str).unique())
    fig, axes = plt.subplots(max(1, len(niches)), 1, figsize=(7, max(3, 2.5 * len(niches))), squeeze=False)
    for ax, niche in zip(axes.ravel(), niches):
        frame = prevalence.loc[prevalence[a.niche_key].astype(str).eq(niche)]
        for i, condition in enumerate(conditions):
            values = frame.loc[frame[a.condition_key].eq(condition), "proportion"]
            ax.scatter([i] * len(values), values, label=condition, alpha=0.8)
        ax.set(title=f"Niche {niche}", ylabel="Patient spot proportion", xticks=range(len(conditions)), xticklabels=conditions)
    fig.tight_layout()
    fig.savefig(a.output_dir / "patient_niche_prevalence.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    config.update({"n_patients": int(obs[a.patient_key].nunique()), "n_niches": int(obs[a.niche_key].nunique()), "conditions": conditions, "versions": {"python": platform.python_version(), "anndata": ad.__version__, "scanpy": sc.__version__}})
    (a.output_dir / "run_metadata.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
