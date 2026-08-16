#!/usr/bin/env python3
"""Compute and validate per-section Squidpy neighborhood enrichment."""

from __future__ import annotations

import argparse
import inspect
import json
import platform
from pathlib import Path

SEED = 1729


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("--cluster-key", required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--spatial-key", default="spatial")
    p.add_argument("--n-perms", type=int, default=1000)
    p.add_argument("--coord-type", choices=["generic", "grid"], default="generic")
    p.add_argument("--delaunay", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def validate(a: argparse.Namespace) -> None:
    if not a.input.is_file():
        raise FileNotFoundError(a.input)
    if a.n_perms < 100:
        raise ValueError("Use at least 100 permutations; 1000 or more is recommended")


def main() -> None:
    a = args()
    validate(a)
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(config, indent=2))
        return

    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import squidpy as sq
    from scipy.sparse.csgraph import connected_components

    np.random.seed(SEED)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_h5ad(a.input)
    if a.cluster_key not in adata.obs or a.spatial_key not in adata.obsm:
        raise KeyError(f"Missing obs[{a.cluster_key!r}] or obsm[{a.spatial_key!r}]")
    coords = np.asarray(adata.obsm[a.spatial_key])
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError(f"Invalid spatial coordinates: {coords.shape}")
    labels = adata.obs[a.cluster_key]
    if labels.isna().any() or labels.nunique() < 2:
        raise ValueError("Cluster labels must be complete and contain at least two groups")
    counts = labels.value_counts().sort_index()
    counts.to_csv(a.output_dir / "cluster_counts.tsv", sep="\t", header=["n_spots"])

    sq.gr.spatial_neighbors(
        adata, spatial_key=a.spatial_key, coord_type=a.coord_type, delaunay=a.delaunay
    )
    graph = adata.obsp["spatial_connectivities"]
    n_components, component = connected_components(graph, directed=False)
    component_sizes = np.bincount(component).tolist()
    kwargs = {"n_perms": a.n_perms}
    signature = inspect.signature(sq.gr.nhood_enrichment)
    if "seed" in signature.parameters:
        kwargs["seed"] = SEED
    sq.gr.nhood_enrichment(adata, cluster_key=a.cluster_key, **kwargs)
    result = adata.uns[f"{a.cluster_key}_nhood_enrichment"]
    categories = list(adata.obs[a.cluster_key].astype("category").cat.categories)
    pd.DataFrame(result["zscore"], index=categories, columns=categories).to_csv(
        a.output_dir / "neighborhood_enrichment_zscore.tsv", sep="\t"
    )
    if "count" in result:
        pd.DataFrame(result["count"], index=categories, columns=categories).to_csv(
            a.output_dir / "neighborhood_adjacency_counts.tsv", sep="\t"
        )
    sq.pl.nhood_enrichment(adata, cluster_key=a.cluster_key, show=False)
    plt.gcf().savefig(a.output_dir / "neighborhood_enrichment.png", dpi=180, bbox_inches="tight")
    plt.close("all")
    adata.write_h5ad(a.output_dir / "neighborhood_results.h5ad", compression="gzip")
    config.update({
        "seed": SEED, "shape": list(adata.shape), "n_clusters": int(labels.nunique()),
        "graph_components": int(n_components), "graph_component_sizes": component_sizes,
        "versions": {"python": platform.python_version(), "anndata": ad.__version__, "scanpy": sc.__version__, "squidpy": sq.__version__},
    })
    (a.output_dir / "run_metadata.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
