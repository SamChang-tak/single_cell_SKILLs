#!/usr/bin/env python3
"""Validate and summarize an inferred spatial activity AnnData matrix by group."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("--group-key", required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--top-n", type=int, default=10)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def main() -> None:
    a = args()
    if not a.input.is_file():
        raise FileNotFoundError(a.input)
    if a.top_n < 1:
        raise ValueError("--top-n must be positive")
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(config, indent=2))
        return

    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    adata = sc.read_h5ad(a.input)
    if a.group_key not in adata.obs:
        raise KeyError(f"Missing obs[{a.group_key!r}]")
    if adata.obs[a.group_key].isna().any() or adata.obs[a.group_key].nunique() < 2:
        raise ValueError("Grouping variable must be complete and contain at least two groups")
    matrix = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    if matrix.shape != adata.shape or not np.isfinite(matrix).all():
        raise ValueError(f"Invalid activity matrix: {matrix.shape}")
    frame = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
    frame[a.group_key] = adata.obs[a.group_key].astype(str).to_numpy()
    means = frame.groupby(a.group_key, observed=True).mean()
    counts = adata.obs[a.group_key].value_counts().sort_index()
    ranked = []
    for group, values in means.iterrows():
        for activity, value in values.abs().sort_values(ascending=False).head(a.top_n).items():
            ranked.append({"group": group, "activity": activity, "mean_activity": values[activity], "absolute_mean": abs(values[activity])})

    a.output_dir.mkdir(parents=True, exist_ok=True)
    means.to_csv(a.output_dir / "mean_activity_by_group.tsv", sep="\t")
    counts.to_csv(a.output_dir / "group_counts.tsv", sep="\t", header=["n_observations"])
    pd.DataFrame(ranked).to_csv(a.output_dir / "top_activities_by_group.tsv", sep="\t", index=False)
    selected = list(dict.fromkeys(row["activity"] for row in ranked))
    fig, ax = plt.subplots(figsize=(max(7, 0.35 * len(selected)), max(4, 0.45 * len(means))))
    image = ax.imshow(means[selected].to_numpy(), aspect="auto", cmap="coolwarm")
    ax.set_xticks(range(len(selected)), selected, rotation=90)
    ax.set_yticks(range(len(means)), means.index)
    ax.set(title=f"Mean inferred activity by {a.group_key}")
    fig.colorbar(image, ax=ax, label="Mean inferred activity")
    fig.tight_layout()
    fig.savefig(a.output_dir / "activity_heatmap.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    config.update({"shape": list(adata.shape), "n_groups": int(means.shape[0]), "finite": True, "warning": "Activity scores are model-inferred; interpret with resource coverage and target-gene support.", "versions": {"python": platform.python_version(), "anndata": ad.__version__, "scanpy": sc.__version__}})
    (a.output_dir / "run_metadata.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
