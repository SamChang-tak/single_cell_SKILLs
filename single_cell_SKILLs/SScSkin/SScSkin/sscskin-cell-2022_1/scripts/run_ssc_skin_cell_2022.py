#!/usr/bin/env python3
"""Reproduce the Spatial_SSc_project 01_cell_2022 Scanpy workflow."""

from __future__ import annotations

import argparse
import json
import platform
from pathlib import Path

SEED = 1729

ANNOTATION_MAP = {
    **{x: "T cells" for x in ["sT_Effector", "T_Effector", "T_GD", "T_Naive", "sTreg_CXCR4", "sTreg", "Treg", "sT"]},
    **{x: "NK cells" for x in ["NK", "NK_XCL1", "NK_XCL1_CXCR4"]},
    "pDC_CXCR4": "pDC",
    **{x: "cDC" for x in ["DC", "DC_CCL22", "DC_CXCL10", "DC_XCR1"]},
    **{x: "Macrophages" for x in ["Mf_TREM2", "Mf"]},
    **{x: "Monocytes" for x in ["Mo_CD16", "M_CD16_IL1B", "M_IL1B", "Mo"]},
    "LC": "Langerhans cells",
    **{x: "Mast cells" for x in ["Mast", "Mast_CLC"]},
    **{x: "B cells" for x in ["B", "B_CXCR4"]},
    "Plasma": "Plasma cells",
    **{x: "Fibroblasts" for x in ["Fibro_ACTA2", "Fibro_Bad", "Fibro_COCH", "Fibro_COMP", "Fibro_IGFBP2", "Fibro_LGR5", "Fibro_MYOC1", "Fibro_MYOC2", "Fibro_POSTN", "Fibro_PTGDS", "Fibro_POSTN_PTGDS"]},
    **{x: "Pericytes" for x in ["Peri_RGS5", "Peri_TGFBI"]},
    **{x: "Endothelial" for x in ["Vascular_ACKR1", "Vascular_RBP7", "Lymphatic_TFF3"]},
    "Melanocytes_MLANA": "Melanocytes",
    **{x: "Keratinocytes" for x in ["KRT1_KRT10", "KRT14_ACTA2", "KRT14_S100A2_GJA1"]},
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--annotation", type=Path, required=True)
    p.add_argument("--sra-run-table", type=Path, required=True)
    p.add_argument("--rna-dir", type=Path, required=True)
    p.add_argument("--patient-workbook", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--annotation-index-col", type=int, default=0)
    p.add_argument("--min-genes", type=int, default=300)
    p.add_argument("--min-counts", type=int, default=500)
    p.add_argument("--max-genes", type=int, default=5000)
    p.add_argument("--max-pct-mt", type=float, default=30.0)
    p.add_argument("--dry-run", action="store_true")
    return p.parse_args()


def validate(a: argparse.Namespace) -> None:
    for path in [a.annotation, a.sra_run_table, a.patient_workbook]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if not a.rna_dir.is_dir():
        raise NotADirectoryError(a.rna_dir)
    if a.output_dir.resolve() in {a.rna_dir.resolve(), a.rna_dir.resolve().parent}:
        raise ValueError("Use a separate output directory; do not overwrite source data")


def load_patient_table(path, pd):
    df = pd.read_excel(path, header=1)
    if not {"PID", "Disease"}.issubset(df.columns):
        df.columns = df.iloc[0]
        df = df.iloc[1:].copy()
    if not {"PID", "Disease"}.issubset(df.columns):
        raise KeyError("Patient workbook must contain PID and Disease")
    out = df[["PID", "Disease"]].rename(columns={"PID": "patient_id", "Disease": "condition"})
    out["patient_id"] = out["patient_id"].astype(str)
    out["condition"] = out["condition"].astype(str)
    return out.drop_duplicates("patient_id")


def save_current_figure(path, plt):
    plt.gcf().savefig(path, dpi=180, bbox_inches="tight")
    plt.close("all")


def main() -> None:
    a = parse_args()
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
    plots = a.output_dir / "plots"
    tables = a.output_dir / "tables"
    plots.mkdir(exist_ok=True)
    tables.mkdir(exist_ok=True)

    annotation = pd.read_csv(a.annotation, sep="\t", index_col=a.annotation_index_col)
    annotation.index = annotation.index.astype(str)
    if not annotation.index.is_unique or "annotation" not in annotation.columns:
        raise ValueError("Annotation index must be unique and include an 'annotation' column")
    run_info = pd.read_csv(a.sra_run_table)
    required = {"Sample Name", "selection_marker", "source_name", "Tissue", "PATIENT_ID", "Organism"}
    missing = required.difference(run_info.columns)
    if missing:
        raise KeyError(f"Missing SRA columns: {sorted(missing)}")

    datasets = []
    for path in sorted(x for x in a.rna_dir.iterdir() if x.is_file()):
        frame = pd.read_csv(path, sep="\t", index_col=0).transpose()
        frame.index = frame.index.astype(str)
        absent = frame.index.difference(annotation.index)
        if len(absent):
            raise KeyError(f"{path.name}: {len(absent)} barcodes absent from annotation table")
        sample_name = path.name.split("_")[0]
        sample = run_info.loc[run_info["Sample Name"].astype(str).eq(sample_name)]
        if len(sample) != 1:
            raise ValueError(f"{path.name}: expected one SRA row for {sample_name}, found {len(sample)}")
        row = sample.iloc[0]
        x = sparse.csr_matrix(frame.to_numpy())
        sample_adata = ad.AnnData(x, obs=annotation.loc[frame.index].copy(), var=pd.DataFrame(index=frame.columns.astype(str)))
        sample_adata.obs["sample_name"] = sample_name
        for target, source in [("selection_marker", "selection_marker"), ("source_name", "source_name"), ("tissue", "Tissue"), ("patient_id", "PATIENT_ID"), ("organism", "Organism")]:
            sample_adata.obs[target] = str(row[source])
        datasets.append(sample_adata)
    if not datasets:
        raise RuntimeError("No count files found")

    adata = ad.concat(datasets, label="batch", index_unique="_")
    adata = adata[adata.obs["annotation"].astype(str).ne("_")].copy()
    adata.write_h5ad(a.output_dir / "01_scrna_raw.h5ad", compression="gzip")
    patients = load_patient_table(a.patient_workbook, pd)
    adata.obs = adata.obs.reset_index(names="cell_id").merge(patients, on="patient_id", how="left").set_index("cell_id")
    if adata.obs["condition"].isna().any():
        raise ValueError(f"Missing condition for {adata.obs['condition'].isna().sum()} cells")

    before = adata.n_obs
    keep = (
        adata.obs["condition"].isin(["SSC", "Control"])
        & adata.obs["tissue"].isin(["Skin", "skin"])
        & adata.obs["selection_marker"].isin(["CD45+", "CD90+"])
        & ~adata.obs["annotation"].isin(["UN", "GBP1", "NRXN1"])
    )
    adata = adata[keep].copy()
    sc.pp.filter_cells(adata, min_genes=a.min_genes)
    sc.pp.filter_cells(adata, min_counts=a.min_counts)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[(adata.obs["n_genes_by_counts"] < a.max_genes) & (adata.obs["pct_counts_mt"] < a.max_pct_mt)].copy()
    sc.pp.filter_genes(adata, min_cells=1)
    adata.layers["counts"] = adata.X.copy()
    pd.crosstab([adata.obs["patient_id"], adata.obs["condition"]], adata.obs["selection_marker"]).to_csv(tables / "retained_cells.tsv", sep="\t")
    sc.pl.violin(adata, ["pct_counts_mt", "n_genes_by_counts"], jitter=0.4, multi_panel=True, show=False)
    save_current_figure(plots / "01_qc_violin.png", plt)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata)
    sc.tl.pca(adata, mask_var="highly_variable", svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(adata, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    sc.pl.umap(adata, color=["annotation", "condition", "selection_marker", "pct_counts_mt"], ncols=2, show=False, frameon=False)
    save_current_figure(plots / "02_umap_overview.png", plt)
    adata.write_h5ad(a.output_dir / "02_scrna_skin.h5ad", compression="gzip")

    detailed = adata.obs["annotation"].astype(str)
    mapped = detailed.map(ANNOTATION_MAP)
    adata.obs["cell_type"] = mapped.fillna(detailed).astype("category")
    unmatched = sorted(detailed[mapped.isna()].unique())
    pd.Series(unmatched, name="unmatched_annotation").to_csv(tables / "unmatched_annotations.tsv", sep="\t", index=False)
    sc.pl.umap(adata, color="cell_type", legend_loc="on data", show=False, frameon=False)
    save_current_figure(plots / "03_umap_cell_type.png", plt)
    sc.tl.rank_genes_groups(adata, "cell_type", method="wilcoxon")
    markers = sc.get.rank_genes_groups_df(adata, group=None)
    markers.to_csv(tables / "cell_type_markers_complete.tsv", sep="\t", index=False)
    markers.groupby("group", observed=True).head(20).to_csv(tables / "cell_type_markers_top20.tsv", sep="\t", index=False)
    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=5, standard_scale="var", cmap="Blues", show=False)
    save_current_figure(plots / "04_marker_matrix.png", plt)
    adata.write_h5ad(a.output_dir / "03_scrna_skin_annotated.h5ad", compression="gzip")

    config.update({
        "seed": SEED, "cells_before_cohort_filter": before, "cells_final": adata.n_obs,
        "genes_final": adata.n_vars, "unmatched_annotations": unmatched,
        "versions": {"python": platform.python_version(), "scanpy": sc.__version__, "anndata": ad.__version__},
    })
    (a.output_dir / "run_metadata.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(config, indent=2))


if __name__ == "__main__":
    main()
