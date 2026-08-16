#!/usr/bin/env python3
"""Finalize an all-cell-type scVI DE output directory after computation."""
from __future__ import annotations

import argparse
import json
import re
from html import escape
from pathlib import Path

import pandas as pd
from PIL import Image


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    out = args.output_dir

    # RGB PNGs are more portable across Preview, browsers, and Office software.
    for path in out.rglob("*.png"):
        with Image.open(path) as image:
            if image.mode != "RGB":
                background = Image.new("RGB", image.size, "white")
                if "A" in image.getbands():
                    background.paste(image, mask=image.getchannel("A"))
                else:
                    background.paste(image.convert("RGB"))
                background.save(path, format="PNG", optimize=True)

    summary = pd.read_csv(out / "de_summary_by_celltype.csv")
    audit = pd.read_csv(out / "celltype_audit.csv")
    config = json.loads((out / "run_config.json").read_text())
    edge_rel = "pseudobulk_edgeR" if (out / "pseudobulk_edgeR").is_dir() else "pseudobulk"
    nominal_rel = (f"{edge_rel}/pseudobulk_nominal"
                   if (out / edge_rel / "pseudobulk_nominal").is_dir()
                   else "pseudobulk_nominal")
    rows = []
    for row in summary.itertuples(index=False):
        folder = f"celltypes/{slug(row.cell_type)}"
        rows.append(
            {
                "cell_type": row.cell_type,
                "disease_cells": row.disease_cells,
                "control_cells": row.control_cells,
                "disease_donors": row.disease_donors,
                "control_donors": row.control_donors,
                "DE_up": row.de_up,
                "DE_down": row.de_down,
                "DE_total": row.de_total,
                "outputs": (
                    f'<a href="{folder}/de_results.csv">all genes CSV</a> · '
                    f'<a href="{folder}/de_genes_fdr.csv">FDR genes CSV</a> · '
                    f'<a href="{folder}/volcano.png">volcano</a>'
                ),
            }
        )
    result_table = pd.DataFrame(rows).to_html(index=False, escape=False)
    pseudobulk = pd.read_csv(out / edge_rel / "pseudobulk_summary_by_celltype.csv")
    pseudobulk_rows = []
    for row in pseudobulk.itertuples(index=False):
        folder = f"{edge_rel}/celltypes/{slug(row.cell_type)}"
        pseudobulk_rows.append(
            {
                "cell_type": row.cell_type,
                "genes_tested": row.genes_tested,
                "disease_donors": row.disease_donors,
                "control_donors": row.control_donors,
                "DE_up": row.up,
                "DE_down": row.down,
                "DE_total": row.total,
                "replication_flag": "LOW: <3 donors in one group" if row.low_replication else "OK",
                "outputs": (
                    f'<a href="{folder}/edger_results.csv">full edgeR CSV</a> · '
                    f'<a href="{folder}/edger_significant.csv">significant CSV</a> · '
                    f'<a href="{folder}/sample_metadata.csv">sample audit</a> · '
                    f'<a href="{folder}/pseudobulk_counts.csv.gz">counts</a> · '
                    f'<a href="{folder}/volcano.png">volcano</a>'
                ),
            }
        )
    pseudobulk_table = pd.DataFrame(pseudobulk_rows).to_html(index=False, escape=False)
    nominal = pd.read_csv(out / nominal_rel / "nominal_summary_by_celltype.csv")
    nominal_rows = []
    for row in nominal.itertuples(index=False):
        folder = f"{nominal_rel}/celltypes/{slug(row.cell_type)}"
        nominal_rows.append(
            {
                "cell_type": row.cell_type,
                "genes_tested": row.genes_tested,
                "nominal_up": row.nominal_up,
                "nominal_down": row.nominal_down,
                "nominal_total": row.nominal_total,
                "disease_donors": row.disease_donors,
                "control_donors": row.control_donors,
                "outputs": (
                    f'<a href="{folder}/edger_results_nominal_annotated.csv">annotated full CSV</a> · '
                    f'<a href="{folder}/nominal_significant.csv">nominal genes CSV</a> · '
                    f'<a href="{folder}/volcano_nominal.png">nominal volcano</a>'
                ),
            }
        )
    nominal_table = pd.DataFrame(nominal_rows).to_html(index=False, escape=False)
    limma = pd.read_csv(out / "pseudobulk_limma" / "limma_summary_by_celltype.csv")
    limma_rows = []
    for row in limma.itertuples(index=False):
        folder = f"pseudobulk_limma/celltypes/{slug(row.cell_type)}"
        limma_rows.append(
            {
                "cell_type": row.cell_type,
                "genes_tested": row.genes_tested,
                "FDR_up": row.fdr_up,
                "FDR_down": row.fdr_down,
                "FDR_total": row.fdr_total,
                "nominal_up": row.nominal_up,
                "nominal_down": row.nominal_down,
                "nominal_total": row.nominal_total,
                "donors": f"{row.disease_donors} vs {row.control_donors}",
                "replication_flag": "LOW" if row.low_replication else "OK",
                "outputs": (
                    f'<a href="{folder}/limma_results.csv">full CSV</a> · '
                    f'<a href="{folder}/limma_fdr_significant.csv">FDR CSV</a> · '
                    f'<a href="{folder}/limma_nominal_significant.csv">nominal CSV</a> · '
                    f'<a href="{folder}/volcano_fdr.png">FDR volcano</a> · '
                    f'<a href="{folder}/volcano_nominal.png">nominal volcano</a>'
                ),
            }
        )
    limma_table = pd.DataFrame(limma_rows).to_html(index=False, escape=False)
    edge_top = pd.read_csv(out / "top20_edgeR_B_and_plasma_up_down.csv")
    edge_top["method"] = "edgeR"
    other_top = pd.read_csv(out / "top20_limma_scvi_B_and_plasma_up_down.csv")
    top = pd.concat(
        [edge_top[["method", "cell_type", "direction", "rank", "gene"]],
         other_top[["method", "cell_type", "direction", "rank", "gene"]]],
        ignore_index=True,
    )
    method_names = {"edgeR": "edgeR pseudobulk", "limma-voom": "limma-voom pseudobulk",
                    "scVI": "scVI sensitivity analysis"}
    top_blocks = []
    for method in ("edgeR", "limma-voom", "scVI"):
        top_blocks.append(f"<h3>{method_names[method]}</h3>")
        for cell_type in ("B cell", "plasma cell"):
            for direction in ("Up", "Down"):
                genes = top.loc[
                    top.method.eq(method) & top.cell_type.eq(cell_type) & top.direction.eq(direction)
                ].sort_values("rank").gene.tolist()
                top_blocks.append(
                    f"<p><b>{escape(cell_type)} {direction.lower()}:</b> "
                    f"{', '.join(escape(str(g)) for g in genes)}.</p>"
                )
    top20_html = "\n".join(top_blocks)
    overlap_rows = []
    for cell_type in ("B cell", "plasma cell"):
        for direction in ("Up", "Down"):
            sets = {
                method: set(top.loc[top.method.eq(method) & top.cell_type.eq(cell_type) &
                                    top.direction.eq(direction), "gene"])
                for method in ("edgeR", "limma-voom", "scVI")
            }
            shared = sorted(sets["edgeR"] & sets["limma-voom"])
            overlap_rows.append({"cell_type": cell_type, "direction": direction,
                                 "edgeR–limma overlap": len(shared),
                                 "shared pseudobulk genes": ", ".join(shared) or "None",
                                 "edgeR–scVI overlap": len(sets["edgeR"] & sets["scVI"]),
                                 "limma–scVI overlap": len(sets["limma-voom"] & sets["scVI"])})
    overlap_table = pd.DataFrame(overlap_rows).to_html(index=False)
    broad = summary.loc[summary.de_total >= 0.8 * config["genes_tested"], "cell_type"].tolist()
    broad_note = ""
    if broad:
        names = ", ".join(escape(x) for x in broad)
        broad_note = (
            f"<p class='warning'><b>Broad DE flag:</b> {names} had at least 80% of tested genes called DE. "
            "This can reflect strong global/compositional or technical shifts; interpret gene-level calls cautiously "
            "and confirm them with donor-aware pseudobulk.</p>"
        )

    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>GSE136103 donor-aware DE</title>
<style>
body{{font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;max-width:1200px;margin:30px auto;padding:0 20px;color:#202124}}
img{{max-width:100%;height:auto}} table{{border-collapse:collapse;font-size:13px;display:block;overflow-x:auto}}
th,td{{padding:6px 9px;border:1px solid #ddd}} th{{background:#f3f5f7}} code{{background:#eee;padding:2px 4px}}
.warning{{background:#fff4ce;border-left:5px solid #d79b00;padding:10px 12px}} .plot{{margin:12px 0 28px}}
</style></head><body>
<h1>GSE136103 disease-versus-control differential expression</h1>
<p><b>Contrast:</b> {escape(config['disease'])} minus {escape(config['control'])} within each annotated {escape(config['tissue'])} cell type. Positive LFC is higher in disease.</p>
<h2>Primary methodology: donor-aware edgeR pseudobulk</h2>
<p>Raw UMI counts were restricted to liver cells and aggregated by <b>cell type × biological donor</b>. Each donor therefore contributes one independent pseudobulk sample per cell type. Within each cell type, low-expression genes were removed with edgeR <code>filterByExpr</code>; library composition was normalized with TMM; negative-binomial dispersions were estimated with robust estimation; and the disease coefficient was tested using edgeR's quasi-likelihood framework (<code>glmQLFit</code>/<code>glmQLFTest</code>). Positive <code>logFC</code> denotes higher expression in fibrotic liver. Primary DE calls require edgeR FDR &lt; 0.05 and |log2FC| ≥ 1. Volcano y-axes display <code>−log10(FDR)</code>, so the horizontal line and colored calls use the same multiple-testing criterion; nominal raw P-values remain available in the CSV tables. No donor blocking term is added because every donor occurs in only one condition—the donor-level count matrix itself supplies biological replication.</p>
<p><b>Replication:</b> most contrasts contain four fibrotic and five healthy donors. Hepatocytes contain four versus four. Mast cells contain only two fibrotic versus five healthy donors and are explicitly flagged as low-replication.</p>
<p><a href="{edge_rel}/pseudobulk_summary_by_celltype.csv">Pseudobulk summary CSV</a> · <a href="{edge_rel}/all_celltypes_edger_results.csv.gz">Combined edgeR results (gzip)</a></p>
{pseudobulk_table}
<div class="plot"><a href="{edge_rel}/combined_pseudobulk_volcano.png"><img src="{edge_rel}/combined_pseudobulk_volcano.png" alt="Combined donor-level pseudobulk volcano plots"></a></div>
<h2>Complementary donor-aware methodology: limma-voom pseudobulk</h2>
<p>The same cell type × donor raw-count pseudobulks and condition design were analyzed with limma-voom. Genes were filtered with edgeR <code>filterByExpr</code>, library composition was normalized with TMM, and <code>voom</code> estimated the mean–variance relationship and observation-level precision weights. Disease effects were fitted with <code>lmFit</code> and moderated using robust <code>eBayes</code>. Positive log2FC denotes fibrotic liver minus healthy. FDR-controlled limma calls require adjusted P &lt; 0.05 and |log2FC| ≥ 1; exploratory nominal calls require raw P &lt; 0.05 and the same effect threshold. FDR volcano plots use <code>−log10(adjusted P)</code>; nominal plots use <code>−log10(raw P)</code>.</p>
<p><a href="pseudobulk_limma/limma_summary_by_celltype.csv">limma summary CSV</a> · <a href="pseudobulk_limma/all_celltypes_limma_results.csv.gz">combined full results</a> · <a href="pseudobulk_limma/all_celltypes_limma_fdr_significant.csv.gz">combined FDR genes</a> · <a href="pseudobulk_limma/all_celltypes_limma_nominal_significant.csv.gz">combined nominal genes</a></p>
{limma_table}
<div class="plot"><a href="pseudobulk_limma/combined_limma_fdr_volcano.png"><img src="pseudobulk_limma/combined_limma_fdr_volcano.png" alt="Combined limma FDR volcano plots"></a></div>
<div class="plot"><a href="pseudobulk_limma/combined_limma_nominal_volcano.png"><img src="pseudobulk_limma/combined_limma_nominal_volcano.png" alt="Combined limma nominal volcano plots"></a></div>
<h2>Exploratory pseudobulk results: nominal significance</h2>
<p>This additional result set uses the same donor-level edgeR quasi-likelihood models but defines exploratory calls as <b>raw P &lt; 0.05 and |log2FC| ≥ 1</b>, without multiple-testing correction. Volcano y-axes show <code>−log10(raw P)</code>, and the horizontal line is raw P = 0.05. These candidates may be useful for ranking and hypothesis generation, but they must not be described as FDR-controlled discoveries. The FDR column is retained in every annotated CSV.</p>
<p><a href="{nominal_rel}/nominal_summary_by_celltype.csv">Nominal summary CSV</a> · <a href="{nominal_rel}/all_celltypes_nominal_annotated.csv.gz">Combined annotated results (gzip)</a> · <a href="{nominal_rel}/all_celltypes_nominal_significant.csv.gz">Combined nominal genes (gzip)</a></p>
{nominal_table}
<div class="plot"><a href="{nominal_rel}/combined_nominal_volcano.png"><img src="{nominal_rel}/combined_nominal_volcano.png" alt="Combined nominal pseudobulk volcano plots"></a></div>
<h2>Secondary methodology: scVI cell-level sensitivity analysis</h2>
<p>For sensitivity analysis only, genes detected in at least {config['min_gene_cells']} liver cells were retained, followed by Poisson selection of {config['genes_tested']:,} genes. A shared negative-binomial scVI model was trained for {config['max_epochs']} maximum epochs and reused for all within-cell-type change-mode contrasts. Differential expression used <code>weights={escape(config['weights'])}</code>, <code>delta={config['delta']}</code>, and Bayesian FDR {config['fdr_target']}. Because the default 5,000 posterior samples yield exact <code>proba_not_de=0</code> for saturated genes, visualization uses <code>−log10(max(proba_not_de, 10⁻⁴))</code>, capped at 4. Hollow downward triangles mark saturated genes at the cap. The cap changes only visualization, not the saved posterior results or Bayesian FDR calls.</p>
<p><a href="run_config.json">scVI run configuration</a> · <a href="training_history.csv">Training history</a> · <a href="celltype_audit.csv">Cell-type audit CSV</a> · <a href="de_summary_by_celltype.csv">scVI summary CSV</a> · <a href="all_celltypes_de_results.csv.gz">Combined scVI results (gzip)</a> · <a href="scvi_plot_saturation_audit.csv">Saturation audit</a></p>
<h2>Model diagnostics and embedding</h2>
<div class="plot"><a href="convergence.png"><img src="convergence.png" alt="scVI convergence"></a></div>
<div class="plot"><a href="latent_umap.png"><img src="latent_umap.png" alt="scVI latent UMAP"></a></div>
<h2>Secondary scVI results and downloads</h2>
{broad_note}{result_table}
<div class="plot"><a href="combined_volcano.png"><img src="combined_volcano.png" alt="Combined volcano plots"></a></div>
<h2>Cell-type audit</h2>
{audit.to_html(index=False)}
<h2>Interpretation limits</h2>
<ul><li>The donor-aware pseudobulk analysis is primary; scVI cell-level DE is retained only as a sensitivity analysis.</li>
<li>Condition remains inseparable from donor identity because donors are nested within condition; results represent between-donor condition differences without repeated measures.</li>
<li>Only two fibrotic donors contribute mast-cell pseudobulks, so that contrast has limited power and unstable dispersion estimation.</li>
<li>Results cover only the {config['genes_tested']:,} genes selected for scVI training.</li>
<li>Posterior DE probabilities are Bayesian quantities, not frequentist p-values.</li>
<li>Mesothelial cells and other underpowered or unassigned groups were skipped as recorded in the audit.</li>
<li>PBMC cells were excluded to avoid tissue/source confounding.</li></ul>
<h2>Appendix: top-20 B-cell and plasma-cell genes across methods</h2>
<p>Genes are ranked within direction by raw quasi-likelihood P-value for edgeR, raw moderated P-value for limma-voom, and posterior DE probability for scVI. None of the B-cell or plasma-cell pseudobulk genes pass FDR &lt; 0.05, so those lists are exploratory rankings rather than discoveries. scVI is a cell-level sensitivity analysis and is not directly equivalent to donor-level pseudobulk inference.</p>
<p><a href="top20_edgeR_B_and_plasma_up_down.csv">Download edgeR rankings</a> · <a href="top20_limma_scvi_B_and_plasma_up_down.csv">Download limma/scVI rankings</a></p>
{top20_html}
<h3>Brief comparison</h3>
<p>edgeR and limma-voom show substantial concordance: their top-20 lists share 8 B-cell up, 7 B-cell down, 12 plasma-cell up, and 11 plasma-cell down genes. Shared pseudobulk candidates are the most reproducible exploratory signals. scVI shows essentially no agreement with the donor-level methods—only one edgeR–scVI overlap and one limma–scVI overlap across all four lists—while returning cross-cell-type markers. This reinforces using edgeR and limma as the primary evidence and retaining scVI only as a sensitivity analysis.</p>
{overlap_table}
<h2>Appendix: biological interpretation of shared plasma-cell candidates</h2>
<p>These genes suggest several biological themes, but not a formally significant pathway enrichment—the plasma-cell comparison has no FDR-significant genes.</p>
<ul>
<li><b>ECM remodeling/fibrosis:</b> <i>SPP1</i> encodes osteopontin, which signals through integrins and CD44 and participates in extracellular-matrix interactions. <a href="https://reactome.org/content/detail/R-HSA-2752125">Reactome</a></li>
<li><b>Lipid handling:</b> <i>FABP5</i> supports intracellular fatty-acid transport and metabolic adaptation. <a href="https://pubmed.ncbi.nlm.nih.gov/28219080/">PubMed</a></li>
<li><b>Inflammatory stress:</b> <i>MAPKAPK2</i> is the downstream MK2 component of p38-MAPK signaling, regulating inflammatory cytokine production, RNA stability, and stress responses. <a href="https://pubmed.ncbi.nlm.nih.gov/34559922/">PubMed</a></li>
<li><b>ROS and oxidative damage:</b> <i>CYBA</i> is part of the NOX2 NADPH-oxidase complex; <i>NUDT1</i> sanitizes oxidized nucleotide pools. <a href="https://pubmed.ncbi.nlm.nih.gov/40227295/">NOX2 review</a>; <a href="https://pubmed.ncbi.nlm.nih.gov/28035004/">NUDT1 study</a></li>
<li><b>Immune adhesion/signaling:</b> <i>ADGRG5/GPR114</i> is an adhesion GPCR expressed in leukocyte lineages. <a href="https://pubmed.ncbi.nlm.nih.gov/27832495/">PubMed</a></li>
<li><b>Antibody biology:</b> <i>IGHV3-48</i> represents immunoglobulin heavy-chain rearrangement and is consistent with plasma/B-cell identity.</li>
<li><b>Chromatin/cellular activation:</b> <i>H2AFZ</i> is a histone variant; <i>C17ORF62</i> and <i>NUDT1</i> may reflect cellular activation or proliferation.</li>
<li><b>Poorly characterized:</b> <i>CHID1</i>, <i>LINC01353</i>, and <i>RP11-16E12</i> cannot currently support a confident pathway assignment.</li>
</ul>
<p><b>Overall interpretation:</b> an SPP1–integrin/ECM, lipid-metabolic, p38–MK2 inflammatory, and oxidative-stress program.</p>
<p class="warning"><b>Cell-identity caution:</b> <i>SPP1</i>, <i>FABP5</i>, <i>CYBA</i>, <i>ADGRG5</i>, and <i>CHID1</i> also look relatively myeloid-like. In a plasma-cell cluster, this could indicate ambient RNA, doublets, or annotation heterogeneity. Confirm that these genes co-express within individual cells alongside <i>JCHAIN</i>, <i>MZB1</i>, <i>SDC1</i>, <i>XBP1</i>, and immunoglobulin genes before claiming a plasma-cell fibrotic program.</p>
</body></html>
"""
    (out / "index.html").write_text(html)


if __name__ == "__main__":
    main()
