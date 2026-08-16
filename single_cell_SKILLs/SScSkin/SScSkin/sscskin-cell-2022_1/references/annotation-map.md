# Detailed-to-broad annotation map

Use exact label matching. Preserve and report any label not listed here.

| Broad cell type | Detailed annotations |
|---|---|
| T cells | `sT_Effector`, `T_Effector`, `T_GD`, `T_Naive`, `sTreg_CXCR4`, `sTreg`, `Treg`, `sT` |
| NK cells | `NK`, `NK_XCL1`, `NK_XCL1_CXCR4` |
| pDC | `pDC_CXCR4` |
| cDC | `DC`, `DC_CCL22`, `DC_CXCL10`, `DC_XCR1` |
| Macrophages | `Mf_TREM2`, `Mf` |
| Monocytes | `Mo_CD16`, `M_CD16_IL1B`, `M_IL1B`, `Mo` |
| Langerhans cells | `LC` |
| Mast cells | `Mast`, `Mast_CLC` |
| B cells | `B`, `B_CXCR4` |
| Plasma cells | `Plasma` |
| Fibroblasts | `Fibro_ACTA2`, `Fibro_Bad`, `Fibro_COCH`, `Fibro_COMP`, `Fibro_IGFBP2`, `Fibro_LGR5`, `Fibro_MYOC1`, `Fibro_MYOC2`, `Fibro_POSTN`, `Fibro_PTGDS`, `Fibro_POSTN_PTGDS` |
| Pericytes | `Peri_RGS5`, `Peri_TGFBI` |
| Endothelial | `Vascular_ACKR1`, `Vascular_RBP7`, `Lymphatic_TFF3` |
| Melanocytes | `Melanocytes_MLANA` |
| Keratinocytes | `KRT1_KRT10`, `KRT14_ACTA2`, `KRT14_S100A2_GJA1` |

The upstream marker dot plot also groups selected genes under labels such as
`T_NK_Cells`, `Macrophage/DC`, `Secretory`, and `Neural`; those are visualization groups,
not values assigned to `adata.obs['cell_type']` by the annotation-collapse notebook.
