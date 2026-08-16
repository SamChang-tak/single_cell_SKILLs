#!/usr/bin/env python3
"""Reproduce the GSE138669 nc_2021_Tabib Scanpy baseline workflow."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

SEED = 1729
PUBLISHED_CLUSTER_MAP = {
    **{str(x): "Keratinocytes" for x in [0, 2, 5, 7, 17, 19, 20, 23]},
    **{str(x): "Endothelial" for x in [1, 14, 26]},
    **{str(x): "Fibroblasts" for x in [3, 4, 6, 12, 21]},
    "8": "Macrophages/cDC",
    **{str(x): "Pericytes" for x in [9, 11, 13]},
    "10": "T cells", "15": "Neural", "16": "NK cells", "18": "Smooth Muscle",
    "22": "Secretory cells", "24": "Mast cells", "25": "Melanocytes", "27": "cDC",
}


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sra-run-table", type=Path, required=True)
    p.add_argument("--h5-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--sample-column", default="Sample Name")
    p.add_argument("--condition-column", default="CONDITION")
    p.add_argument("--leiden-resolution", type=float, default=1.0)
    p.add_argument("--max-pct-mt", type=float, default=5.0)
    p.add_argument("--skip-published-annotation", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def validate(a: argparse.Namespace) -> None:
    if not a.sra_run_table.is_file():
        raise FileNotFoundError(a.sra_run_table)
    if not a.h5_dir.is_dir():
        raise NotADirectoryError(a.h5_dir)
    if a.output_dir.resolve() in {a.h5_dir.resolve(), a.h5_dir.resolve().parent}:
        raise ValueError("Use a separate output directory")
    if a.leiden_resolution <= 0 or not 0 < a.max_pct_mt <= 100:
        raise ValueError("Invalid resolution or mitochondrial threshold")


def save_figure(path, plt) -> None:
    plt.gcf().savefig(path, dpi=180, bbox_inches="tight")
    plt.close("all")


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
    from scipy import sparse

    np.random.seed(SEED)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    plots, tables = a.output_dir / "plots", a.output_dir / "tables"
    plots.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)
    run_info = pd.read_csv(a.sra_run_table)
    for key in [a.sample_column, a.condition_column]:
        if key not in run_info.columns:
            raise KeyError(f"Missing SRA column: {key}")

    datasets = []
    for path in sorted(a.h5_dir.glob("*.h5")):
        sample_name = path.name.split("_")[0]
        sample = run_info.loc[run_info[a.sample_column].astype(str).eq(sample_name)]
        if len(sample) != 1:
            raise ValueError(f"{path.name}: expected one metadata row for {sample_name}, found {len(sample)}")
        sample_adata = sc.read_10x_h5(path)
        sample_adata.var_names_make_unique()
        if not sparse.issparse(sample_adata.X) or sample_adata.X.data.size == 0:
            raise ValueError(f"{path.name}: missing sparse counts")
        if (sample_adata.X.data < 0).any() or not np.allclose(sample_adata.X.data, np.round(sample_adata.X.data)):
            raise ValueError(f"{path.name}: X is not non-negative integer-like counts")
        sc.pp.filter_cells(sample_adata, min_genes=200)
        sample_adata.obs["sample"] = sample_name
        sample_adata.obs["condition"] = str(sample.iloc[0][a.condition_column])
        datasets.append(sample_adata)
    if not datasets:
        raise RuntimeError("No .h5 count matrices found")

    adata = ad.concat(datasets, label="batch", index_unique="_")
    adata.write_h5ad(a.output_dir / "01_scrna_raw.h5ad", compression="gzip")
    cells_loaded = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=300)
    sc.pp.filter_cells(adata, min_counts=500)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs["n_genes_by_counts"] < 5000) & (adata.obs["pct_counts_mt"] < a.max_pct_mt)].copy()
    adata.layers["counts"] = adata.X.copy()
    pd.crosstab([adata.obs["sample"], adata.obs["condition"]], adata.obs["batch"]).to_csv(
        tables / "retained_cells.tsv", sep="\t"
    )
    sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], multi_panel=True, jitter=0.3, show=False)
    save_figure(plots / "01_qc_violin.png", plt)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.tl.pca(adata, mask_var="highly_variable", svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(adata, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.pl.umap(adata, color=["sample", "batch", "condition", "n_genes_by_counts"], ncols=2, show=False)
    save_figure(plots / "02_umap_batch_condition.png", plt)
    adata.write_h5ad(a.output_dir / "02_scrna_processed.h5ad", compression="gzip")

    sc.tl.leiden(adata, resolution=a.leiden_resolution, random_state=SEED)
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    cluster_markers = sc.get.rank_genes_groups_df(adata, group=None)
    cluster_markers.to_csv(tables / "cluster_markers_complete.tsv", sep="\t", index=False)
    cluster_markers.groupby("group", observed=True).head(20).to_csv(tables / "cluster_markers_top20.tsv", sep="\t", index=False)
    pd.crosstab([adata.obs["leiden"], adata.obs["condition"]], adata.obs["sample"]).to_csv(
        tables / "cluster_sample_condition_counts.tsv", sep="\t"
    )
    sc.pl.umap(adata, color=["leiden", "condition"], legend_loc="on data", show=False, frameon=False)
    save_figure(plots / "03_umap_leiden_condition.png", plt)

    unknown_clusters = sorted(set(adata.obs["leiden"].astype(str)) - set(PUBLISHED_CLUSTER_MAP))
    if a.skip_published_annotation:
        adata.obs["cell_type"] = adata.obs["leiden"].astype(str).astype("category")
    else:
        mapped = adata.obs["leiden"].astype(str).map(PUBLISHED_CLUSTER_MAP)
        adata.obs["cell_type"] = mapped.fillna("Unmapped").astype("category")
        parent = adata.obs["leiden"].astype(str).eq("8")
        if parent.sum() >= 20:
            myeloid = adata[parent].copy()
            sc.tl.leiden(myeloid, resolution=0.1, key_added="leiden_0.1", random_state=SEED)
            refined = myeloid.obs["leiden_0.1"].astype(str).map({"0": "cDC", "1": "Macrophages"}).fillna("Macrophages/cDC")
            adata.obs["cell_type"] = adata.obs["cell_type"].astype(str)
            adata.obs.loc[myeloid.obs_names, "cell_type"] = refined
            adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")

    pd.Series(unknown_clusters, name="unmapped_cluster").to_csv(tables / "unmapped_clusters.tsv", sep="\t", index=False)
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False, frameon=False)
    save_figure(plots / "04_umap_cell_type.png", plt)
    sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")
    type_markers = sc.get.rank_genes_groups_df(adata, group=None)
    type_markers.to_csv(tables / "cell_type_markers_complete.tsv", sep="\t", index=False)
    type_markers.groupby("group", observed=True).head(20).to_csv(tables / "cell_type_markers_top20.tsv", sep="\t", index=False)
    adata.write_h5ad(a.output_dir / "03_scrna_clustered_annotated.h5ad", compression="gzip")

    config.update({
        "seed": SEED, "cells_loaded": cells_loaded, "cells_final": adata.n_obs,
        "genes_final": adata.n_vars, "unknown_clusters": unknown_clusters,
        "batch_correction": "none (matches upstream executable code)",
        "annotation_warning": "Published numeric cluster labels require marker validation before transfer.",
        "versions": {"python": platform.python_version(), "scanpy": sc.__version__, "anndata": ad.__version__},
    })
    (a.output_dir / "run_metadata.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
