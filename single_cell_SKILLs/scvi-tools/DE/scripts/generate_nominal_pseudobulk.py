#!/usr/bin/env python3
"""Generate exploratory nominal-P pseudobulk tables and volcano plots."""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")


def draw(ax, table: pd.DataFrame, cell_type: str, annotate: bool = False) -> tuple[int, int]:
    call = table.PValue.lt(.05) & table.logFC.abs().ge(1)
    up = call & table.logFC.gt(0)
    down = call & table.logFC.lt(0)
    score = -np.log10(np.maximum(table.PValue.astype(float), np.nextafter(0, 1)))
    ax.scatter(table.loc[~call, "logFC"], score[~call], s=7 if annotate else 3,
               c="#c8c8c8", alpha=.45, linewidths=0)
    ax.scatter(table.loc[up, "logFC"], score[up], s=11 if annotate else 5,
               c="#d62728", alpha=.75, linewidths=0, label=f"Nominal up ({up.sum():,})")
    ax.scatter(table.loc[down, "logFC"], score[down], s=11 if annotate else 5,
               c="#1f77b4", alpha=.75, linewidths=0, label=f"Nominal down ({down.sum():,})")
    ax.axvline(1, ls="--", c="black", lw=1 if annotate else .7)
    ax.axvline(-1, ls="--", c="black", lw=1 if annotate else .7)
    ax.axhline(-np.log10(.05), ls=":", c="#555555", lw=.9 if annotate else .7,
               label="Raw P = 0.05" if annotate else None)
    ax.set_xlabel("edgeR log2FC: Fibrotic liver minus Healthy" if annotate else "edgeR log2FC")
    ax.set_ylabel("−log10(edgeR raw QL P-value)")
    ax.set_title(f"{cell_type} — exploratory nominal P" if annotate else
                 f"{cell_type}\nup {up.sum():,} / down {down.sum():,}")
    if annotate:
        labels = table.loc[call].nsmallest(8, "PValue")
        for row in labels.itertuples(index=False):
            y = -np.log10(max(float(row.PValue), np.nextafter(0, 1)))
            ax.annotate(row.gene, (row.logFC, y), xytext=(3, 3), textcoords="offset points", fontsize=7)
        ax.legend(frameon=False, fontsize=9)
    return int(up.sum()), int(down.sum())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    args = parser.parse_args()
    root = args.results_dir
    primary = root / "pseudobulk"
    out = root / "pseudobulk_nominal"
    (out / "celltypes").mkdir(parents=True, exist_ok=True)
    primary_summary = pd.read_csv(primary / "pseudobulk_summary_by_celltype.csv")
    summaries = []
    combined = []
    combined_sig = []
    plot_tables = []

    for meta in primary_summary.itertuples(index=False):
        cell_type = meta.cell_type
        table = pd.read_csv(primary / "celltypes" / slug(cell_type) / "edger_results.csv")
        call = table.PValue.lt(.05) & table.logFC.abs().ge(1)
        table["nominal_p05_abs_log2fc1"] = call
        folder = out / "celltypes" / slug(cell_type)
        folder.mkdir(parents=True, exist_ok=True)
        table.to_csv(folder / "edger_results_nominal_annotated.csv", index=False)
        table.loc[call].sort_values("PValue").to_csv(folder / "nominal_significant.csv", index=False)
        fig, ax = plt.subplots(figsize=(8.5, 6.2))
        up, down = draw(ax, table, cell_type, annotate=True)
        fig.tight_layout()
        fig.savefig(folder / "volcano_nominal.png", dpi=170, facecolor="white")
        plt.close(fig)
        summaries.append({"cell_type": cell_type, "genes_tested": len(table), "nominal_up": up,
                          "nominal_down": down, "nominal_total": int(call.sum()),
                          "disease_donors": meta.disease_donors, "control_donors": meta.control_donors,
                          "low_replication": meta.low_replication})
        combined.append(table.assign(cell_type=cell_type))
        combined_sig.append(table.loc[call].assign(cell_type=cell_type))
        plot_tables.append((cell_type, table))

    summary = pd.DataFrame(summaries).sort_values("nominal_total", ascending=False)
    summary.to_csv(out / "nominal_summary_by_celltype.csv", index=False)
    pd.concat(combined, ignore_index=True).to_csv(out / "all_celltypes_nominal_annotated.csv.gz",
                                                  index=False, compression="gzip")
    pd.concat(combined_sig, ignore_index=True).to_csv(out / "all_celltypes_nominal_significant.csv.gz",
                                                      index=False, compression="gzip")
    cols = 3
    rows = math.ceil(len(plot_tables) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(15, 4.5 * rows), squeeze=False)
    for ax, (cell_type, table) in zip(axes.flat, plot_tables):
        draw(ax, table, cell_type, annotate=False)
    for ax in axes.flat[len(plot_tables):]:
        ax.axis("off")
    fig.suptitle("Exploratory nominal P<0.05: Fibrotic liver vs Healthy", fontsize=16)
    fig.tight_layout()
    fig.savefig(out / "combined_nominal_volcano.png", dpi=180, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
