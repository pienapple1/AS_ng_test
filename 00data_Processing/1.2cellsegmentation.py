###Author Jinpei Lin
###Times 2024.07
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import skimage
import sklearn
import spateo as st
import dynamo as dyn
import cv2
import anndata as ad
import pandas as pd
from scipy import stats
import os
import scanpy as sc
from collections import Counter
import anndata as ad
import loompy as lp
###set Set data path
path_gem="/00gem/"
path_DNA="/01ssDNA_adjust/"
dout="/02cellsegmentation/"
####read file
key="SS200000709TR_D1"
outdir = dout+key+"/"
if not os.path.exists(outdir):
    os.makedirs(outdir)
adata = st.io.read_bgi_agg(
    path_gem+key+'.gem2.gz', path_DNA+key+'_ssDNA.tif',prealigned=True
)
####Since Python needs to flip the image when reading it, if the information of the read image does not correspond to RNA, use .T to flip it.###That‘s importany
adata.layers["stain"]=adata.layers["stain"].T
st.cs.mask_nuclei_from_stain(adata)
st.pl.imshow(adata, 'stain_mask')####Check whether the imaging morphology is consistent with the RNA signal map
st.cs.find_peaks_from_mask(adata, 'stain', 7)
st.cs.watershed(adata, 'stain', 5, out_layer='watershed_labels')
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'watershed_labels', labels=True, alpha=0.5, ax=ax)
st.cs.stardist(adata,equalize=2.0, out_layer='stardist_labels')###This step requires a large amount of computing resources. If the computing resources are not enough to meet the computing needs, an error will be reported.
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'stardist_labels', labels=True, alpha=0.5, ax=ax)
st.cs.mask_nuclei_from_stain(adata)
st.pl.imshow(adata, 'stain_mask')
st.cs.find_peaks_from_mask(adata, 'stain', 7)
st.cs.watershed(adata, 'stain', 5, out_layer='watershed_labels')
    
st.cs.augment_labels(adata, 'watershed_labels', 'stardist_labels', out_layer='augmented_labels')

fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'augmented_labels', labels=True, alpha=0.5, ax=ax)
st.cs.mask_cells_from_stain(adata, out_layer='stain_cell_mask')
st.cs.watershed(
    adata, 'stain',
    mask_layer='stain_cell_mask',
    markers_layer='augmented_labels',
    out_layer='cell_labels',
)
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'cell_labels', labels=True, alpha=0.5, ax=ax)
st.cs.expand_labels(adata, 'augmented_labels', distance=10, max_area=600)#####After locating the nucleus, expand the nuclear area to capture the genetic information of the cytoplasm.
fig, ax = st.pl.imshow(adata, 'stain', save_show_or_return='return')
st.pl.imshow(adata, 'augmented_labels_expanded', labels=True, alpha=0.5, ax=ax)
cv2.imwrite(outdir+key+"_cell_shape_exp.tif", adata.layers['augmented_labels_expanded'].astype(np.uint8))####in this reaserach,we select layer["augmented_labels_expanded"]  As the result of cell segmentation
#####save cellsegmentation imformation h5ad
adata.write(outdir+key + "_cellsegmentation_shape.h5ad")

###creat cell*gene scanpy file
cell_adata = st.io.read_bgi(
    path_gem+key+'.gem2.gz',
    segmentation_adata=adata,
    labels_layer='augmented_labels_expanded',
)
cell_adata.uns["seg"] = adata.layers["augmented_labels_expanded"].copy()####Save array information of cell segmentation
###add coodrnate to obs
cell_adata.obs['x'] = cell_adata.obsm['spatial'][:,0]
cell_adata.obs['y'] = cell_adata.obsm['spatial'][:,1]
    
    
###papre to seurat file meta.data and gene*cell datframe
row_attrs={
    "Gene":np.array(cell_adata.var.index)
}
col_attrs = { 
    "CellID": np.array(cell_adata.obs.index),
    "nGene": np.array( np.sum(cell_adata.X.transpose()>0 , axis=0)).flatten(),
    "nUMI": np.array( np.sum(cell_adata.X.transpose() , axis=0)).flatten(),
}
lp.create(outdir+key+'_data.loom', cell_adata.X.transpose(), row_attrs, col_attrs)
cell_adata.obs.to_csv(outdir+key+'_metadata.csv')
    
####save file base on h5ad file
cell_adata.write(outdir+key + "_cell_gene.h5ad")