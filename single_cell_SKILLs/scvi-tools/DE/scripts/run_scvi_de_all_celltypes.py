#!/usr/bin/env python3
"""Fit one SCVI model and run condition DE independently within each cell type."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import math
import platform
import re
from html import escape
from pathlib import Path


def arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--condition-key", default="Condition")
    p.add_argument("--disease", required=True)
    p.add_argument("--control", required=True)
    p.add_argument("--celltype-key", default="Cell type (standardized)")
    p.add_argument("--tissue-key", default="Tissue")
    p.add_argument("--tissue", default="Liver")
    p.add_argument("--donor-key", default="Subject ID")
    p.add_argument("--count-layer")
    p.add_argument("--n-top-genes", type=int, default=4000)
    p.add_argument("--min-gene-cells", type=int, default=20)
    p.add_argument("--min-cells-per-group", type=int, default=25)
    p.add_argument("--max-epochs", type=int, default=150)
    p.add_argument("--early-stopping-patience", type=int, default=20)
    p.add_argument("--weights", choices=("importance", "uniform"), default="importance")
    p.add_argument("--delta", type=float, default=.25)
    p.add_argument("--fdr-target", type=float, default=.05)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def slug(x): return re.sub(r"[^a-z0-9]+", "_", str(x).lower()).strip("_")


def main():
    a=arguments(); out=a.output_dir
    if out.exists() and any(out.iterdir()): raise FileExistsError(f"Output directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True); (out/"celltypes").mkdir()
    import anndata as ad
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scvi
    from scipy import sparse

    data=ad.read_h5ad(a.input)
    for key in (a.condition_key,a.celltype_key,a.tissue_key,a.donor_key):
        if key not in data.obs: raise KeyError(f"Missing obs key: {key}")
    counts=data.layers[a.count_layer] if a.count_layer else data.X
    vals=counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    probe=vals[:min(len(vals),500_000)]
    if len(probe) and (np.any(probe<0) or not np.allclose(probe,np.round(probe))):
        raise ValueError("Selected matrix is not raw nonnegative integer-like counts")
    keep=(data.obs[a.tissue_key].astype(str)==a.tissue) & data.obs[a.condition_key].astype(str).isin([a.disease,a.control])
    data=data[keep].copy(); data.X=(counts[keep.to_numpy()].copy() if a.count_layer else data.X)
    data.obs[a.condition_key]=data.obs[a.condition_key].astype(str).astype("category")
    data.obs[a.celltype_key]=data.obs[a.celltype_key].astype(str).astype("category")
    sc.pp.filter_genes(data,min_cells=a.min_gene_cells)
    data.layers["counts"]=sparse.csr_matrix(data.X.copy())
    scvi.settings.seed=a.seed
    scvi.data.poisson_gene_selection(data,n_top_genes=min(a.n_top_genes,data.n_vars),inplace=True)
    data=data[:,data.var["highly_variable"]].copy(); data.layers["counts"]=sparse.csr_matrix(data.layers["counts"])

    ct_counts=pd.crosstab(data.obs[a.celltype_key].astype(str),data.obs[a.condition_key].astype(str)).reindex(columns=[a.disease,a.control],fill_value=0)
    donor_counts=data.obs.groupby([a.celltype_key,a.condition_key],observed=True)[a.donor_key].nunique().unstack(fill_value=0).reindex(columns=[a.disease,a.control],fill_value=0)
    audit=ct_counts.rename(columns={a.disease:"disease_cells",a.control:"control_cells"}).join(
        donor_counts.rename(columns={a.disease:"disease_donors",a.control:"control_donors"}))
    audit["status"]=np.where((audit.disease_cells>=a.min_cells_per_group)&(audit.control_cells>=a.min_cells_per_group)&(audit.index!="Unassigned"),"run","skip")
    audit.index.name="cell_type"; audit.to_csv(out/"celltype_audit.csv")

    scvi.model.SCVI.setup_anndata(data,layer="counts")
    model=scvi.model.SCVI(data,gene_likelihood="nb")
    model.train(max_epochs=a.max_epochs,check_val_every_n_epoch=1,early_stopping=True,
                early_stopping_patience=a.early_stopping_patience,early_stopping_monitor="elbo_validation")
    model.save(out/"model",save_anndata=False)
    hist={}
    for k,v in model.history.items():
        if isinstance(v,pd.DataFrame): hist[k]=v.iloc[:,0].reset_index(drop=True)
        elif isinstance(v,pd.Series): hist[k]=v.reset_index(drop=True)
    history=pd.DataFrame(hist); history.to_csv(out/"training_history.csv",index_label="epoch")
    fig,ax=plt.subplots(figsize=(8,5))
    for k in ("elbo_train","elbo_validation"):
        if k in history: ax.plot(history.index,history[k],label=k)
    ax.set_yscale("log"); ax.set_xlabel("Epoch"); ax.set_ylabel("ELBO loss (log scale)"); ax.set_title("SCVI convergence"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(out/"convergence.png",dpi=180); plt.close(fig)
    data.obsm["X_scVI"]=model.get_latent_representation(); sc.pp.neighbors(data,use_rep="X_scVI"); sc.tl.umap(data,random_state=a.seed)
    fig,axes=plt.subplots(1,2,figsize=(16,7))
    sc.pl.umap(data,color=a.condition_key,ax=axes[0],show=False,title="Disease vs control")
    sc.pl.umap(data,color=a.celltype_key,ax=axes[1],show=False,title="Cell types")
    fig.tight_layout(); fig.savefig(out/"latent_umap.png",dpi=180); plt.close(fig)

    summaries=[]; combined=[]; plot_data=[]
    cond=data.obs[a.condition_key].astype(str).to_numpy(); ctype=data.obs[a.celltype_key].astype(str).to_numpy()
    for ct,row in audit.iterrows():
        if row.status!="run": continue
        idx1=(ctype==ct)&(cond==a.disease); idx2=(ctype==ct)&(cond==a.control)
        de=model.differential_expression(idx1=idx1,idx2=idx2,mode="change",weights=a.weights,
            filter_outlier_cells=(a.weights=="importance"),batch_correction=False,
            delta=a.delta,fdr_target=a.fdr_target)
        de.index.name="gene"; fdr=[c for c in de if c.startswith("is_de_fdr_")]
        sig=de[fdr[0]].astype(bool) if fdr else de.proba_de.ge(.95)
        up=sig & de.lfc_mean.gt(0); down=sig & de.lfc_mean.lt(0)
        folder=out/"celltypes"/slug(ct); folder.mkdir()
        de.to_csv(folder/"de_results.csv"); de.loc[sig].to_csv(folder/"de_genes_fdr.csv")
        score=-np.log10(np.clip(de.proba_not_de.astype(float),1e-300,1))
        fig,ax=plt.subplots(figsize=(8,6)); ax.scatter(de.loc[~sig,"lfc_mean"],score[~sig],s=6,c="#bdbdbd",alpha=.5,linewidths=0)
        ax.scatter(de.loc[up,"lfc_mean"],score[up],s=8,c="#d62728",alpha=.7,linewidths=0,label=f"Up ({up.sum()})")
        ax.scatter(de.loc[down,"lfc_mean"],score[down],s=8,c="#1f77b4",alpha=.7,linewidths=0,label=f"Down ({down.sum()})")
        ax.axvline(a.delta,ls="--",c="black",lw=1); ax.axvline(-a.delta,ls="--",c="black",lw=1)
        ax.set_xlabel(f"Posterior mean LFC: {a.disease} minus {a.control}"); ax.set_ylabel("-log10 posterior non-DE probability")
        ax.set_title(ct); ax.legend(frameon=False); fig.tight_layout(); fig.savefig(folder/"volcano.png",dpi=160); plt.close(fig)
        de2=de.reset_index(); de2.insert(0,"cell_type",ct); combined.append(de2)
        summaries.append({"cell_type":ct,"disease_cells":int(row.disease_cells),"control_cells":int(row.control_cells),
                          "disease_donors":int(row.disease_donors),"control_donors":int(row.control_donors),
                          "de_up":int(up.sum()),"de_down":int(down.sum()),"de_total":int(sig.sum())})
        plot_data.append((ct,de.lfc_mean.to_numpy(),score.to_numpy(),sig.to_numpy(),up.to_numpy(),down.to_numpy()))
    summary=pd.DataFrame(summaries).sort_values("de_total",ascending=False); summary.to_csv(out/"de_summary_by_celltype.csv",index=False)
    pd.concat(combined,ignore_index=True).to_csv(out/"all_celltypes_de_results.csv.gz",index=False,compression="gzip")

    n=len(plot_data); cols=3; rows=math.ceil(n/cols); fig,axes=plt.subplots(rows,cols,figsize=(15,4.5*rows),squeeze=False)
    for ax,(ct,lfc,score,sig,up,down) in zip(axes.flat,plot_data):
        ax.scatter(lfc[~sig],score[~sig],s=3,c="#c8c8c8",alpha=.45,linewidths=0)
        ax.scatter(lfc[up],score[up],s=4,c="#d62728",alpha=.65,linewidths=0)
        ax.scatter(lfc[down],score[down],s=4,c="#1f77b4",alpha=.65,linewidths=0)
        ax.axvline(a.delta,ls="--",c="black",lw=.7); ax.axvline(-a.delta,ls="--",c="black",lw=.7); ax.set_title(ct)
        ax.set_xlabel("Posterior mean LFC"); ax.set_ylabel("-log10 posterior non-DE probability")
    for ax in axes.flat[n:]: ax.axis("off")
    fig.suptitle(f"{a.disease} vs {a.control} within liver cell types",fontsize=16); fig.tight_layout()
    fig.savefig(out/"combined_volcano.png",dpi=180,facecolor="white"); plt.close(fig)
    data.write_h5ad(out/"processed.h5ad",compression="gzip")
    config={k:(str(v) if isinstance(v,Path) else v) for k,v in vars(a).items()}
    config.update(cells=data.n_obs,genes_tested=data.n_vars,celltypes_run=len(summary),
                  batch_correction=False,batch_note="Donor is nested within condition; no identifiable shared donor batch correction",
                  versions={"python":platform.python_version(),"scvi_tools":importlib.metadata.version("scvi-tools"),"anndata":importlib.metadata.version("anndata")})
    (out/"run_config.json").write_text(json.dumps(config,indent=2)+"\n")
    html=f"""<!doctype html><html><head><meta charset='utf-8'><title>GSE136103 scVI DE</title><style>body{{font:15px/1.5 -apple-system,BlinkMacSystemFont,Segoe UI,sans-serif;max-width:1180px;margin:30px auto;padding:0 20px}}img{{max-width:100%}}table{{border-collapse:collapse;font-size:13px}}th,td{{padding:6px 9px;border:1px solid #ddd}}th{{background:#f3f5f7}}code{{background:#eee;padding:2px 4px}}</style></head><body><h1>GSE136103 disease-versus-control scVI DE</h1><p><b>Contrast:</b> {escape(a.disease)} minus {escape(a.control)} within each annotated liver cell type. Positive LFC is higher in disease.</p><h2>Methodology</h2><p>Raw counts were restricted to {escape(a.tissue)} cells, genes detected in at least {a.min_gene_cells} cells, then {data.n_vars:,} genes were selected with Poisson gene selection. One negative-binomial SCVI model was trained and reused for all within-cell-type change-mode contrasts. DE used weights={a.weights}, delta={a.delta}, and Bayesian FDR={a.fdr_target}. Donor was not entered as a batch covariate because donors are fully nested within disease/control; such adjustment is not identifiable.</p><img src='convergence.png'><img src='latent_umap.png'><h2>DE results</h2>{summary.to_html(index=False)}<img src='combined_volcano.png'><h2>Cell-type audit</h2>{audit.reset_index().to_html(index=False)}<h2>Interpretation limits</h2><ul><li>Cells are not independent biological replicates. Validate condition effects with donor-aware pseudobulk.</li><li>Results apply only to the {data.n_vars:,} genes used to train scVI.</li><li>Posterior DE probabilities are not frequentist p-values.</li><li>Mesothelial and any other underpowered groups are skipped according to the audit.</li><li>PBMC cells were excluded to avoid tissue/source confounding.</li></ul></body></html>"""
    (out/"index.html").write_text(html)


if __name__=="__main__": main()
