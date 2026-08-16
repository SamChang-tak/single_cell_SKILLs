#!/usr/bin/env python3
"""Integrate annotated Gur 2022 and Tabib 2021 SSc skin AnnData objects."""

from __future__ import annotations

import argparse
import inspect
import json
import platform
from pathlib import Path

SEED = 1729
MARKERS = {
    "B cells": ["CD79A"], "Endothelial": ["VWF", "PECAM1", "CDH5", "ERG", "KDR", "TEK"],
    "Fibroblasts": ["DCN", "COL1A1", "PDGFRA", "PCOLCE", "THY1", "ACKR3"],
    "Keratinocytes": ["KRT1", "KRT5"], "Langerhans cells": ["CD1A"],
    "Macrophages": ["CD68", "CD163"], "Mast cells": ["TPSAB1"], "Melanocytes": ["PMEL"],
    "Monocytes": ["FCGR3A", "C5AR1", "FCN1"], "NK cells": ["NKG7"], "Neural": ["PLP1"],
    "Pericytes": ["RGS5", "PDGFRB"], "Plasma cells": ["SDC1"],
    "Secretory cells": ["DCD", "SCGB2A2", "MUCL1", "CALML5"],
    "Smooth Muscle": ["ACTA2", "DES", "SMTN", "TAGLN", "TPM2"],
    "T cells": ["CD2", "CD3D"], "cDC": ["CD1C", "CD1B", "CD1E"], "pDC": ["NRP1", "IRF8"],
}


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gur2022", type=Path, required=True)
    p.add_argument("--tabib2021", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--harmony-key", default="data")
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def validate(a: argparse.Namespace) -> None:
    for path in [a.gur2022, a.tabib2021]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if a.gur2022.resolve() == a.tabib2021.resolve():
        raise ValueError("The two inputs must be different files")


def save_figure(path: Path, plt) -> None:
    plt.gcf().savefig(path, dpi=180, bbox_inches="tight")
    plt.close("all")


def condition_map(series):
    mapping = {"Control": "Healthy", "CONTROL": "Healthy", "Healthy": "Healthy", "SSC": "SSc", "SSc": "SSc"}
    result = series.astype(str).map(mapping)
    if result.isna().any():
        raise ValueError(f"Unrecognized condition labels: {sorted(series[result.isna()].astype(str).unique())}")
    return result


def main() -> None:
    a = args()
    validate(a)
    config = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    if a.dry_run:
        print(json.dumps(config, indent=2))
        return

    import anndata as ad
    import harmonypy as hm
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
    gur, tabib = sc.read_h5ad(a.gur2022), sc.read_h5ad(a.tabib2021)
    for name, obj, required in [
        ("Gur2022", gur, {"patient_id", "condition", "cell_type"}),
        ("Tabib2021", tabib, {"sample", "condition", "cell_type"}),
    ]:
        missing = required.difference(obj.obs.columns)
        if missing or "counts" not in obj.layers:
            raise KeyError(f"{name}: missing obs={sorted(missing)}, counts layer={('counts' in obj.layers)}")
        counts = obj.layers["counts"]
        values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
        if values.size == 0 or (values < 0).any() or not np.allclose(values, np.round(values)):
            raise ValueError(f"{name}: layers['counts'] is not non-negative integer-like counts")

    gene_summary = {
        "gur_genes": gur.n_vars, "tabib_genes": tabib.n_vars,
        "common_genes": len(gur.var_names.intersection(tabib.var_names)),
        "gur_only_genes": len(gur.var_names.difference(tabib.var_names)),
        "tabib_only_genes": len(tabib.var_names.difference(gur.var_names)),
    }
    if gene_summary["common_genes"] == 0:
        raise ValueError("Inputs have no common gene identifiers")
    pd.DataFrame([gene_summary]).to_csv(tables / "feature_overlap.tsv", sep="\t", index=False)
    gur.obs = gur.obs[["patient_id", "condition", "cell_type"]].copy()
    tabib.obs = tabib.obs[["sample", "condition", "cell_type"]].rename(columns={"sample": "patient_id"}).copy()
    gur.obs["condition"], tabib.obs["condition"] = condition_map(gur.obs["condition"]), condition_map(tabib.obs["condition"])
    type_summary = pd.DataFrame({
        "Gur2022": pd.Series(sorted(gur.obs["cell_type"].astype(str).unique())),
        "Tabib2021": pd.Series(sorted(tabib.obs["cell_type"].astype(str).unique())),
    })
    type_summary.to_csv(tables / "source_cell_types.tsv", sep="\t", index=False)

    adata = ad.concat(
        {"Gur2022": gur, "Tabib2021": tabib}, join="inner", merge="same",
        label="data", index_unique="-",
    )
    del gur, tabib
    adata.X = adata.layers["counts"].copy()
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    pd.crosstab([adata.obs["cell_type"], adata.obs["condition"]], adata.obs["data"]).to_csv(
        tables / "cell_type_condition_source_counts.tsv", sep="\t"
    )
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    sc.tl.pca(adata, mask_var="highly_variable", random_state=SEED)
    sc.pp.neighbors(adata, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    adata.obsm["X_umap_uncorrected"] = adata.obsm["X_umap"].copy()
    sc.pl.umap(adata, color=["data", "condition", "cell_type"], ncols=3, show=False, frameon=False)
    save_figure(plots / "01_umap_pre_harmony.png", plt)

    harmony_kwargs = {}
    if "random_state" in inspect.signature(hm.run_harmony).parameters:
        harmony_kwargs["random_state"] = SEED
    harmony = hm.run_harmony(
        np.asarray(adata.obsm["X_pca"], dtype=np.float64), adata.obs, [a.harmony_key], **harmony_kwargs
    )
    harmony_matrix = np.asarray(harmony.Z_corr)
    corrected = harmony_matrix.T if harmony_matrix.shape[1] == adata.n_obs else harmony_matrix
    if corrected.shape[0] != adata.n_obs or not np.isfinite(corrected).all():
        raise RuntimeError(f"Invalid Harmony coordinates: {corrected.shape}")
    adata.obsm["X_pca_harmony"] = corrected
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.pl.umap(adata, color=["data", "condition", "cell_type"], ncols=3, show=False, frameon=False)
    save_figure(plots / "02_umap_post_harmony.png", plt)

    sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(tables / "cell_type_markers_complete.tsv", sep="\t", index=False)
    markers.query("logfoldchanges > 0").groupby("group", observed=True).head(5).to_csv(
        tables / "cell_type_markers_top5.tsv", sep="\t", index=False
    )
    available = {group: [g for g in genes if g in adata.var_names] for group, genes in MARKERS.items()}
    available = {group: genes for group, genes in available.items() if genes}
    missing_markers = sorted({g for genes in MARKERS.values() for g in genes if g not in adata.var_names})
    pd.Series(missing_markers, name="missing_marker").to_csv(tables / "missing_markers.tsv", sep="\t", index=False)
    if available:
        sc.pl.dotplot(adata, available, groupby="cell_type", standard_scale="var", show=False)
        save_figure(plots / "03_canonical_marker_dotplot.png", plt)
    metabolic = [g for g in ["NAMPT", "NNMT"] if g in adata.var_names]
    if metabolic:
        adata.obs["condition_cell_type"] = adata.obs["cell_type"].astype(str) + "_" + adata.obs["condition"].astype(str)
        sc.pl.dotplot(adata, metabolic, groupby="condition_cell_type", swap_axes=True, show=False)
        save_figure(plots / "04_nampt_nnmt_dotplot.png", plt)

    adata.write_h5ad(a.output_dir / "integrated.h5ad", compression="gzip")
    config.update({
        "seed": SEED, "shape": list(adata.shape), "feature_overlap": gene_summary,
        "harmony_components": corrected.shape[1], "missing_markers": missing_markers,
        "warning": "Harmony corrects study label only; assess study/condition/patient confounding. Use patient-aware inference.",
        "versions": {"python": platform.python_version(), "scanpy": sc.__version__, "anndata": ad.__version__, "harmonypy": getattr(hm, "__version__", "unknown")},
    })
    (a.output_dir / "run_metadata.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
