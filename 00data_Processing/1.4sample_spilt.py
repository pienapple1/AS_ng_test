import os
import matplotlib.pyplot as plt
import spateo as st
import dynamo as dyn
import scanpy as sc
import numpy as np
import plotly.express as px
import pandas as pd
import cv2


#####chip spilt
path1 = "/01h5ad_raw/"####seurat transfer h5ad
path2 = "/00spilt_sample/"
outdir = "/02csv_spilt_sample/"
for key in [i [:-5] for i in os.listdir(path1)]:
    adata = sc.read_h5ad(path1 +key+".h5ad")
    adata.obsm["spatial"]=adata.obs[["x","y"]].to_numpy()
    adata.uns["__type"] = "UMI"
    mask = cv2.imread(path2+key+".tif",cv2.IMREAD_GRAYSCALE)
    adata.obs['chip_id'] = 0
    for i in range(len(adata)):
        try:
            if mask[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
                adata.obs['chip_id'][i] = 1
        except:
            pass
    st.pl.space(adata, color=['chip_id'], pointsize=0.2)