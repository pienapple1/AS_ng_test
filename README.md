This repository contains the code and analysis pipelines for the manuscript:
Spatiotemporal Control of Vascular Cell Plasticity during Atherosclerosis
<img width="4961" height="2026" alt="figure5_v7" src="https://github.com/user-attachments/assets/5a70dac0-2efa-4593-8776-bf6200d1d105" />


## Software and algorithms

| Software | Reference | Version | URL |
| --- | --- | --- | --- |
| GraphPad Prism | GraphPad Software | V8.0 | <https://www.graphpad.com/> |
| ImageJ | NIH | N/A | <https://imagej.net/ij/> |
| R | R Project | V4.3.1 | <https://www.r-project.org/> |
| Seurat | Stuart et al., 2019 | V4.0.5 | <https://satijalab.org/seurat/> |
| Hotspot | DeTomso & Yosef, 2021 | V0.9 | <https://yoseflab.github.io/Hotspot/> |
| OpenCV | Bradski, 2000 | V4.6.0.66 | <https://opencv.org/> |
| Monocle2 | Qiu et al., 2017 | V1.0.0 | <https://cole-trapnell-lab.github.io/monocle-release/docs/> |
| Spateo | Qiu et al., 2022 | V1.0.2 | <https://spateo-release.readthedocs.io/en/latest/> |
| clusterProfiler | Wu et al., 2021 | V4.6.2 | <https://github.com/YuLab-SMU/clusterProfiler> |
| SAW | BGI | V3.8.1 | <https://github.com/BGIResearch/SAW> |
| HISAT2 | Kim et al., 2019 | V2.1.0 | <http://daehwankimlab.github.io/hisat2/> |
| TopoVelo | Gu et al., 2025 | V0.0.1a1 | <https://github.com/welch-lab/TopoVelo> |
| seaborn | Waskom, M., 2021 | V0.13.2 | <https://github.com/mwaskom/seaborn> |
| Harmony | Korsunsky et al., 2019 | V0.1 | <https://github.com/immunogenomics/harmony> |
| CellChat | Suoqin Jin et al., 2021 | V1.5.0 | <https://github.com/sqjin/CellChat> |
| spacexr | Dylan M. Cable et al., 2021 | V1.0.0 | <https://github.com/dmcable/spacexr> |
| scvelo | Bergen et al., 2020 | V0.2.5 | <https://github.com/theislab/scvelo> |
| fgsea | Gennady Korotkevich et al., 2021 | V1.16.0 | <https://github.com/alserglab/fgsea> |
| GSVA | Hänzelmann S et al., 2013 | V2.4.1 | <https://github.com/rcastelo/GSVA> |

###Directory tree
```text
.
├── 00data_Processing/
│   ├── 1.1Creat_RNA.py
│   ├── 1.2cellsegmentation.py
│   ├── 1.2Create_seurat.r
│   ├── 1.3Creat_bin_data.r
│   ├── 1.3Create_bin_data.R
│   ├── 1.3tissuecut.py
│   ├── 1.4sample_area_spilt_and_digitization.py
│   ├── 1.4sample_spilt.py
│   ├── 1.5hotspot_ann_plot.py
│   ├── 1.5hotspot_data_pre.r
│   ├── 1.6human_ST_data_ann.r
│   ├── 1.6mouse_ST_data_ann.r
│   ├── 1.6mouse_ST_data_cellsubtype.r
│   ├── 1.6scRNA_ref_del.r
│   ├── 1.7gsva_ann.r
│   ├── 1.7hotspot_plaque_related_gene.csv
│   ├── 1.7result_data.RDS
│   ├── 1.8change.py
│   ├── 1.8pseudobulk_of_celltype_area.ipynb
│   └── 1.8pyscenic.sh
│
├── 01figure_Plotting/
│   ├── fig1/
│   │   ├── 1g.RDS
│   │   ├── 1H.RDS
│   │   ├── figure1F.r
│   │   ├── figure1G.py
│   │   ├── figure1H.r
│   │   ├── Figure1I.r
│   │   ├── figure1J.py
│   │   └── figure1K.r
│   ├── fig2/
│   │   ├── 2e.RDS
│   │   ├── 2g.RDS
│   │   ├── 2j_geneset.RDS
│   │   ├── figure2b_hotspot.py
│   │   ├── Figure2b.LPCA.R
│   │   ├── figure2c.r
│   │   ├── figure2d.r
│   │   ├── figure2e.r
│   │   ├── figure2f.py
│   │   ├── figure2g.r
│   │   ├── figure2h.r
│   │   ├── figure2j.r
│   │   └── plaque_gene.csv
│   ├── fig3/
│   │   ├── 3o.RDS
│   │   ├── figure3b.r
│   │   ├── figure3c.r
│   │   ├── figure3d.ipynb
│   │   ├── figure3ef.ipynb
│   │   ├── figure3g.ipynb
│   │   ├── figure3n.r
│   │   ├── figure3o.r
│   │   └── figure3p.ipynb
│   ├── fig4/
│   │   ├── 4a_p.RDS
│   │   ├── 4a.r
│   │   ├── 4a_tfs_score.RDS
│   │   ├── figure4d.ipynb
│   │   ├── scenic_select.RDS
│   │   └── scenic_type.RDS
│   ├── sup1/
│   │   ├── S1A.ipynb
│   │   ├── S1I.ipynb
│   │   ├── S1J.py
│   │   └── S1K.ipynb
│   ├── sup2/
│   │   └── S2.py
│   ├── sup3/
│   │   ├── S3A&S3B.ipynb
│   │   ├── S3C_heatmap.ipynb
│   │   └── S3D.py
│   ├── sup4/
│   │   ├── gs.csv
│   │   ├── let_ALL_save1.RDS
│   │   ├── pca_df.RDS
│   │   ├── sup4b.ipynb
│   │   ├── sup4c.ipynb
│   │   ├── sup4E.ipynb
│   │   ├── sup4F.ipynb
│   │   └── sup4G.ipynb
│   ├── sup5/
│   │   ├── sup5A.ipynb
│   │   ├── sup5C.ipynb
│   │   ├── sup5D,G,H.ipynb
│   │   ├── sup5E,F.ipynb
│   │   └── sup5K.ipynb
│   └── sup6/
│       ├── CELL_PROLIFERATION.csv
│       ├── sup6A.ipynb
│       ├── sup6B.ipynb
│       └── sup6C.ipynb
│
├── LICENSE
└── README.md
00data_Processing/ — Data preprocessing
This directory contains the main preprocessing pipeline. Scripts are numbered roughly in execution order (1.1 → 1.8).
Both Python and R are used.

## Workflow overview

The repository is organized into two stages: (i) spatial transcriptomics and single-cell data processing (`00data_Processing/`),
 and (ii) generation of main and supplementary figures (`01figure_Plotting/`). Scripts are numbered approximately in execution order.

---

## Step 1: Spatial transcriptomics and single-cell data processing

### 1.1 mRNA image generation and registration

The ssDNA nuclear image and the mRNA signal image are registered. mRNA images are generated with `00data_Processing/1.1Creat_RNA.py` (Python 3.8)
 from the Stereo-seq gene expression matrix (`*.txt` or `*.gem`; file format: <https://db.cngb.org/stomics/project/STT0000118/spatialGeneExpression>).
After generation, images are manually aligned in Photoshop or ImageJ using the reference lines provided by Stereo-seq.

### 1.2 Cell segmentation and Seurat object construction

Cell segmentation is performed on the calibrated ssDNA image with `00data_Processing/1.2cellsegmentation.py` (Python 3.8). Inputs
are the calibrated ssDNA image and the gene expression matrix (`*.txt` or `*.gem`). Before running, ensure that the signal image
and the ssDNA image are consistently oriented. Seurat objects are then constructed in R 4.2 from the segmentation-derived metadata
and the cell × gene matrix (`*_data.loom`, `*_metadata.csv`) using `00data_Processing/1.2Create_seurat.r` and `00data_Processing/1.3Creat_bin_data.r`.

### 1.3 Pseudo-cell (bin) construction

When cell segmentation is not required, pseudo-cells (bins) can be generated for rapid downstream analyses such as area statistics
and spatial enrichment. Use R 4.2 with `00data_Processing/1.3Create_bin_data.R` and `00data_Processing/1.3Create_seurat.r`, and
set the bin size as needed (e.g., bin50 = 50 spots ≈ 25 µm × 25 µm). Tissue segmentation is required before bin construction and
can be based on gene capture counts or, as in this study, on the calibrated ssDNA image or the cell-segmentation mask. Tissue
cutting is performed with `00data_Processing/1.3tissuecut.py` (Python 3.8); adjust parameters as needed. To remain consistent
with cell segmentation, the cell-segmentation mask is recommended. Note that `1.3Creat_bin_data.r` and `1.3Create_bin_data.R`
differ only in filename case and are two distinct files on Linux; their roles should be checked before execution.

### 1.4 Sample splitting, vessel region division, and digitization

Because two samples are mounted on one chip, sample splitting is required. A mask is drawn on the calibrated ssDNA image in
Photoshop or ImageJ and processed with `00data_Processing/1.4sample_spilt.py` (Python 3.8). For vessel region division, six
regions were defined and digitized into layers 1–100. Masks were drawn on the calibrated ssDNA image in Photoshop or ImageJ
and saved as `tif`, then processed with `00data_Processing/1.4sample_area_spilt_and_digitization.py` (Python 3.8).

### 1.5 Spatial co-expression module analysis and plaque-related gene identification

Hotspot (v0.9) was used for spatial co-expression module analysis with the top 5,000 highly variable genes. The Seurat object
is used to export the cell × gene expression matrix and spatial coordinates. Python 3.8 and R 4.2 were used. See
`00data_Processing/1.5hotspot_data_pre.r` and `00data_Processing/1.5hotspot_ann_plot.py`.

### 1.6 Cell type annotation

Cell types in human and mouse spatial transcriptomics data were annotated with RCTD. A processed single-cell reference is
required before RCTD. Single-cell reference processing is performed with `00data_Processing/1.6scRNA_ref_del.r` (R 4.2).
RCTD annotation of human and mouse spatial data is implemented in `00data_Processing/1.6human_ST_data_ann.r`,
`00data_Processing/1.6mouse_ST_data_ann.r`, and `00data_Processing/1.6mouse_ST_data_cellsubtype.r`.

### 1.7 GSVA enrichment of atherosclerotic progression

Based on the spatial co-expression modules, 44 gene sets related to atherosclerotic plaque progression were analyzed with GSVA.
Counts for each vessel region were first aggregated into pseudo-bulk profiles. R 4.2 was used; see `00data_Processing/1.7gsva_ann.r`.

### 1.8 Transcription factor analysis

Configure the environment according to the pySCENIC documentation (<https://github.com/aertslab/pySCENIC>). Analysis scripts are `00data_Processing/1.8pseudobulk_of_celltype_area.ipynb`, `00data_Processing/1.8change.py`, and `00data_Processing/1.8pyscenic.sh`.

---

## Step 2: Figure plotting scripts

### Figure 1

- `figure1F.r` — Gene-set scoring of single-cell data with Seurat `AddModuleScore()`. Input: mouse single-cell data processed by `00data_Processing/1.6scRNA_ref_del.r`.
- `figure1G.py` — Spatial cell mask image generation. Inputs: final annotated metadata (e.g., `*.csv`) and the segmentation mask
stored in `adata.uns["seg"]` by `00data_Processing/1.2cellsegmentation.py`. Python 3.8 and OpenCV 4.5.5 are recommended.
- `figure1H.r` — Cell-type compositional changes across stages. Cell proportions are computed per sample and plotted as line charts with ggplot.
- `Figure1I.r` — Visualization of vessel coefficient of variation (CV) and vessel thickness. Vessel length is measured in ImageJ.
Vessel area is derived from the bin50 pseudo-cell file (each bin50 = 625 µm²), and mean thickness = total area / vessel length.
CV is computed from digitized layers 1–100 using sample SD and mean; see `00data_Processing/1.4sample_area_spilt_and_digitization.py` for digitization.
- `figure1J.py` — Same as `figure1G.py`.
- `figure1K.r` — Gene-set scoring of human spatial transcriptomics data with Seurat `AddModuleScore()`. Input: the human Seurat object
after `00data_Processing/`.

### Figure 2

- `figure2b_hotspot.py` — See Section 1.5.
- `Figure2b.LPCA.R` — GSVA using published ChIP-seq data; see the code and Methods.
- `figure2c.r` — Gene-set scoring of mouse spatial data with `AddModuleScore()`. Gene set: `plaque_gene.csv`. Input: the
mouse spatial Seurat object after `00data_Processing/`.
- `figure2d.r` — Plaque progression gene-set scoring with `plaque_gene.csv`, plotted separately for “GC” and “LC”. Same input as above.
- `figure2e.r` — Changes in cell proportions along the atherosclerotic progression score; analysis follows Section 1.7.
- `figure2f.py` — Same principle as `figure1G.py`.
- `figure2g.r` — Mosaic plot of plaque severity across stages and regions in mouse samples. Data processed by `00data_Processing/`.
- `figure2h.r` — Density distribution of cell types from vessel lumen to outer wall across stages. Merge sample metadata after `00data_Processing/`.
- `figure2j.r` — Requires macrophage pseudo-bulk profiles. After `00data_Processing/`, extract macrophages from intima (layers 1–25)
and adventitia (layers 75–100) per sample and region, then build pseudo-bulk. Layer definition follows `00data_Processing/1.4sample_area_spilt_and_digitization.py`.

### Figure 3

- `figure3b.r` — Requires pseudo-bulk profiles per region and cell type (e.g., all SMC in GC-AA of sample 1). See `00data_Processing/1.8pseudobulk_of_celltype_area.ipynb`.
- `figure3c.r` — Uses the final Seurat object after `00data_Processing/`. Extract the normalized gene expression matrix and metadata.
The x-axis intima/adventitia corresponds to layer after inversion (`mer$layer <- abs(mer$layer-100)`).
- `figure3d.ipynb` — Pi16 expression in fibroblasts from GSE155513 and this study’s single-cell data. Processing script:
`00data_Processing/1.XX_Pi16_del.R`.
- `figure3ef.ipynb` — KDE analysis of PI16 and FN1 in human clinical samples. After `00data_Processing/`, extract cell
centroids (x, y) and PI16/FN1 expression, merge with R `cbind`, save as CSV, and run in Python.
- `figure3g.ipynb` — Human clinical data after `00data_Processing/`. Because the human dataset is large, layers are rounded
to integers and pseudo-bulk is built per sample at layers 0–100, then merged; see the code.
- `figure3n.r` — Mouse data after `00data_Processing/`. Select W6, W10, and W20 samples with plaque features, extract mSMC
in intima layers 1–25 and fibroblasts in adventitia layers 75–100, then merge (Seurat `merge()` or rebuild from counts + metadata).
- `figure3o.r` — Uses DEGs from `figure3n.r` to examine their distribution in ECM. ECM file: `3o.RDS`.
- `figure3p.ipynb` — Script provided by a co-author (documentation pending).

### Figure 4

- `4a.r` — Analysis follows `00data_Processing/1.8pseudobulk_of_celltype_area.ipynb`, `00data_Processing/1.8change.py`, and `00data_Processing/1.8pyscenic.sh`. Pseudo-bulk per cell type is required.
- `figure4d.ipynb` — GO enrichment of downstream targets of Sox4. The top 100 weighted targets of Sox4 were used; see the code.

### Supplementary Figure 1 (sup1)

- `S1A.ipynb` — PCA of mouse spatial samples. Two tissue sections were mounted per chip. The top 2,000 highly variable genes
from single-cell data were used. Pseudo-bulk was built for each spatial sample, PCA was run on the 2,000-gene matrix, and
Euclidean distances were computed between the two samples on the same chip and between samples within the same stage.
Within-chip and within-stage distances were compared.
- `S1I.ipynb` — Stacked bar plot of cell proportions in mouse spatial samples. Input: the Seurat object after
`00data_Processing/` or the final metadata.
- `S1J.py` — Spatial cell mask image generation. Inputs: final annotated metadata (e.g., `*.csv`) and the segmentation
mask in `adata.uns["seg"]` from `00data_Processing/1.2cellsegmentation.py`.
- `S1K.ipynb` — Line plot of cell-type proportions across regions in mouse spatial samples (e.g., proportion of SMC in
GC-AA of sample 1). Input: the Seurat object after `00data_Processing/` or the final metadata.

### Supplementary Figure 2 (sup2)

- `S2.py` — Spatial cell mask image generation comparing broad and subtype annotations. Inputs: final annotated metadata (e.g., `*.csv`)
and the segmentation mask in `adata.uns["seg"]` from `00data_Processing/1.2cellsegmentation.py`. Python 3.8 and OpenCV 4.5.5 are recommended.

### Supplementary Figure 3 (sup3)

- `S3A&S3B.ipynb` — Cross-species (mouse–human) single-cell integration with SAMap. Mouse–human ortholog mapping must be built first;
see the SAMap documentation (<https://github.com/atarashansky/SAMap>) and the code.
- `S3C_heatmap.ipynb` — Heatmap of cell-type marker expression in human spatial data. Seurat objects after `00data_Processing/`
are merged, and each cell type is downsampled to 2,000 cells for runtime; see the code.
- `S3D.py` — Spatial cell mask image generation. Inputs: final annotated metadata (e.g., `*.csv`) and the segmentation mask in
`adata.uns["seg"]` from `00data_Processing/1.2cellsegmentation.py`. Python 3.8 and OpenCV 4.5.5 are recommended.

### Supplementary Figure 4 (sup4)

- `sup4b.ipynb` — Spatial expression of plaque-related gene sets in human clinical samples, scored with Seurat `AddModuleScore()`.
Mouse genes must first be converted to human orthologs (see `gs.csv`). Input: the Seurat object after `00data_Processing/`.
- `sup4c.ipynb` — Density distribution of mSMC-1 and mSMC-2 in the LC region along blood-flow direction. The blood-flow column is defined in `00data_Processing/1.4sample_area_spilt_and_digitization.py`.
- `sup4E.ipynb` — Expression patterns of plaque-related gene sets in human clinical samples. Human clinical spatial transcriptomics samples
from Lai et al., Bleckwehl et al., and this study are collected, pseudo-bulk is built, and shared gene sets are extracted for analysis; see the code.
- `sup4F.ipynb` — Molecular staging of plaque severity in human clinical samples with GSVA. Samples from Lai et al., Bleckwehl et al., and this
study are collected, pseudo-bulk is built, and shared gene sets are extracted; see the code.
- `sup4G.ipynb` — Histogram of cell distribution from outer to inner vessel wall across plaque stages in human clinical samples. Each stage is
downsampled to the same cell number; see the code.

### Supplementary Figure 5 (sup5)

- `sup5A.ipynb` — Differential gene analysis of mouse macrophages across stages, displayed as an UpSet plot. All macrophages are extracted
and compared against stage 1 macrophages using Seurat `FindAllMarkers()` with `min.pct = 0.06` and `logfc.threshold = 0.5`; see the code.
- `sup5C.ipynb` — Differential gene analysis between intimal (layers 0–25) and adventitial (layers 75–100) macrophages at stage 1 and
stage 2 using Seurat `FindAllMarkers()`; see the code.
- `sup5D,G,H.ipynb` — D: spatial distribution of intimal and adventitial macrophages at stage 1 and stage 2. G, H: spatial expression
distribution of Mmp12 and Lgals3 along the inner and outer vessel wall; see the code.
- `sup5E,F.ipynb` — E: differential gene analysis of adventitial macrophages across stages. F: differential gene analysis between intimal
and adventitial macrophages within the same molecular stage.
- `sup5K.ipynb` — Expression distribution of Mmp9/Mmp12 in mouse and MMP9/MMP12 in human clinical samples along the inner and outer vessel
wall layers. Layer definition follows `00data_Processing/1.4sample_area_spilt_and_digitization.py`.

### Supplementary Figure 6 (sup6)

- `sup6A.ipynb` — Line plots of cell-type markers along plaque progression in human samples.
Pseudo-bulk construction and plaque scoring follow `sup4F.ipynb`.
- `sup6B.ipynb` — Signaling score of the CELL_PROLIFERATION pathway in mouse samples. Data after `00data_Processing/` are
aggregated into pseudo-bulk per sample and scored with Seurat `AddModuleScore()`.
- `sup6C.ipynb` — Spatial RNA velocity analysis on high-fat-diet 10-week samples. Spliced and unspliced matrices are
required and are retained during cell segmentation. Analysis uses TopoVelo; see the code.

---

## Notes on scripts

Several figure scripts share the same mask-generation logic (`figure1G.py`, `figure1J.py`, `figure2f.py`, `S1J.py`, `S2.py`, and `S3D.py`);
they differ mainly in input metadata and annotation level. The layer system (1–100) used throughout is defined in `00data_Processing/1.4sample_area_spilt_and_digitization.py`; scripts referring to “layer” should use that definition. Some intermediate `.RDS`
objects are included for reproducibility; if they are large, consider Git LFS or regenerate them from the preprocessing scripts.
