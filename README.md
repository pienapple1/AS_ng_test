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
This directory contains the main preprocessing pipeline. Scripts are numbered roughly in execution order (1.1 → 1.8). Both Python and R are used.

File	Description
1.1Creat_RNA.py	Create RNA-level data objects from raw spatial data.
1.2cellsegmentation.py	Cell segmentation for spatial data.
1.2Create_seurat.r	Create Seurat objects from segmented/processed data.
1.3tissuecut.py	Tissue cut / region extraction.
1.3Creat_bin_data.r, 1.3Create_bin_data.R	Create bin-level data objects. Note: filenames differ only by case; on Linux these are two separate files.
1.4sample_spilt.py, 1.4sample_area_spilt_and_digitization.py	Sample splitting, area splitting, and digitization of spatial regions.
1.5hotspot_data_pre.r	Prepare data for Hotspot analysis.
1.5hotspot_ann_plot.py	Hotspot annotation and plotting.
1.6human_ST_data_ann.r, 1.6mouse_ST_data_ann.r	Annotation of human and mouse spatial transcriptomics data.
1.6mouse_ST_data_cellsubtype.r	Mouse spatial data cell-subtype annotation.
1.6scRNA_ref_del.r	scRNA reference processing / deletion of unwanted cells or samples.
1.7gsva_ann.r	GSVA annotation and pathway scoring.
1.7hotspot_plaque_related_gene.csv	Gene list related to plaque/regions used in Hotspot analysis.
1.7result_data.RDS	Intermediate RDS object generated during preprocessing.
1.8change.py	Format conversion / data reshaping.
1.8pseudobulk_of_celltype_area.ipynb	Pseudobulk generation by cell type and area.
1.8pyscenic.sh	Shell script for running pySCENIC.
Suggested execution order: 1.1 → 1.2 → 1.3 → 1.4 → 1.5 → 1.6 → 1.7 → 1.8. Some steps depend on outputs from previous steps; adjust according to your local paths and data availability.

01figure_Plotting/ — Figure plotting
This directory contains scripts for main figures (fig1–fig4) and supplementary figures (sup1–sup6). Script names generally correspond to figure panels (e.g., figure2b_hotspot.py for Figure 2B; S1A.ipynb for Supplementary Figure 1A).

