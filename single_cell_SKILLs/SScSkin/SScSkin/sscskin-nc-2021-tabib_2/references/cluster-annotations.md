# Published cluster annotations and markers

Apply these mappings only after marker validation. They reflect the later
`04_annotation.ipynb` followed by `05_annotation.ipynb`.

| Broad label | Parent Leiden clusters |
|---|---|
| Keratinocytes | 0, 2, 5, 7, 17, 19, 20, 23 |
| Endothelial | 1, 14, 26 |
| Fibroblasts | 3, 4, 6, 12, 21 |
| Macrophages/cDC parent | 8 |
| Pericytes | 9, 11, 13 |
| T cells | 10 |
| Neural | 15 |
| NK cells | 16 |
| Smooth Muscle | 18 |
| Secretory cells | 22 |
| Mast cells | 24 |
| Melanocytes | 25 |
| cDC | 27 |

Within parent cluster 8, the published resolution-0.1 subcluster map is:

- subcluster 0 → cDC
- subcluster 1 → Macrophages

Useful marker panels from the notebooks include:

- Keratinocytes: `KRT1`, `KRT10`, `KRT14`, `KRT5`, `KRT17`, `KRT6A`, `KRT6B`, `KRT2`
- Endothelial: `AQP1`, `IFI27`, `CCL21`, `CLDN5`, `RAMP2`, `PECAM1`, `VWF`
- Fibroblasts: `SFRP2`, `COMP`, `DCN`, `APOD`, `CCL19`, `APOE`, `POSTN`, `COCH`,
  `ASPN`, `ANGPTL7`, `COL1A1`, `PDGFRA`
- Pericytes: `RGS5`, `GEM`, `TAGLN`, `RERGL`, `PDGFRB`
- T/NK: `CD3D`, `CD4`, `CXCR4`, `IL32`
- Macrophage/myeloid: `CCL2`, `CD14`, `CD163`, `CD68`, `FCGR3A`, `AIF1`
- cDC: `CD1C`, `CD74`, `IRF8`
- Mast: `TPSAB1`
- Melanocyte: `MLANA`, `PMEL`, `TYR`
- Lymphatic endothelial: `PDPN`, `PROX1`, `LYVE1`
- Smooth muscle: `ACTA2`, `DES`, `TAGLN`
- Secretory: `DCD`, `KRT19`, `SCGB2A2`, `MUCL1`
- Neural/Schwann: `SOX10`, `PLP1`, `PMP22`, `MPZ`

The notebooks report that canonical pDC, B-cell, and plasma-cell marker panels were not
detected clearly. Treat absence cautiously because dropout, gene naming, and low abundance
can all contribute.
