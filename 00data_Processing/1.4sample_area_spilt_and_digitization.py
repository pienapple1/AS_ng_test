####Author Jinpei Lin
####Times 2024.07

import os
import matplotlib.pyplot as plt
import spateo as st
import dynamo as dyn
import scanpy as sc
import numpy as np
import plotly.express as px
import pandas as pd
import cv2
    
#####GC/LC spilt
path1 = "/DATA/User/liqian7/pienapple/04AS_updata/00data/07spilt_data/07h5ad_digitization/"
path2 = "/DATA/User/liqian7/pienapple/04AS_updata/00data/05mask/02digitization/"
outdir = "/DATA/User/liqian7/pienapple/04AS_updata/00data/07spilt_data/08csv_digitization/"
key="SS200000709TR_D1_1"
adata = sc.read_h5ad(path1 +key+".h5ad")
adata.obsm["spatial"]=adata.obs[["x","y"]].to_numpy()
adata.uns["__type"] = "UMI"
mask = cv2.imread(path2+key+".tif",cv2.IMREAD_GRAYSCALE)#######
adata.obs['roi'] = "GC"
for i in range(len(adata)):
    try:
        if mask[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
            adata.obs['roi'][i] = "LC"
    except:
        pass
st.pl.space(adata, color=['roi'], pointsize=0.2)
#####digitization
# Extract an area of interest
adata.obsm['spatial_bin20'] = adata.obsm['spatial']//20
cluster_label_image_lowres = st.dd.gen_cluster_image(adata, bin_size=1, spatial_key="spatial_bin20", cluster_key='roi', show=False)
cluster_label_list = np.unique(adata[adata.obs['roi'].isin(["LC"]), :].obs["cluster_img_label"])
contours, cluster_image_close, cluster_image_contour = st.dd.extract_cluster_contours(cluster_label_image_lowres, 
                                                                                      cluster_label_list, bin_size=1, k_size=15, show=False)
px.imshow(cluster_image_contour, width=500, height=500)

# User input to specify a gridding direction
pnt_XY = (634,263)###maxcolum maxlayer
pnt_Xy = (624,265)###maxcolum minlayer 
pnt_xY = (1008,311)##mincolum maxlayer
pnt_xy = (1031,323)##mincolum minlayer 
# Digitize the area of interest
st.dd.digitize(
    adata=adata,
    ctrs=contours,
    ctr_idx=0,
    pnt_xy=pnt_xy,
    pnt_xY=pnt_xY,
    pnt_Xy=pnt_Xy,
    pnt_XY=pnt_XY,
    spatial_key="spatial_bin20"
)
# Visualize digitized layers and columns
st.pl.space(
    adata,
    color=['digital_layer', 'digital_column'],
    ncols=2,
    pointsize=0.1,
    show_legend="upper left",
    figsize=(4, 3.5),
    color_key_cmap = "tab20",
    save_show_or_return='return',
)

####division spilt
%matplotlib inline
path1 = "/DATA/User/liqian7/pienapple/04AS_updata/00data/07spilt_data/04h5ad_division/"
path2 = "/DATA/User/liqian7/pienapple/04AS_updata/00data/05mask/01division/"
outdir = "/DATA/User/liqian7/pienapple/04AS_updata/00data/07spilt_data/05csv_division/"
for key in [i [:-5] for i in os.listdir(path1)]:
    adata = sc.read_h5ad(path1 +key+".h5ad")
    print(key)
    adata.obsm["spatial"]=adata.obs[["x","y"]].to_numpy()
    adata.uns["__type"] = "UMI"
    division1 = cv2.imread(path2+key+"_division1.tif",cv2.IMREAD_GRAYSCALE)
    division2 = cv2.imread(path2+key+"_division2.tif",cv2.IMREAD_GRAYSCALE)
    division3 = cv2.imread(path2+key+"_division3.tif",cv2.IMREAD_GRAYSCALE)
    adata.obs['division'] = "others"
    for i in range(len(adata)):
        try:
            if division1[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
                adata.obs['division'][i] = "AA"
            if division2[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
                adata.obs['division'][i] = "AR"
            if division3[int(adata.obsm['spatial'][i][0]), int(adata.obsm['spatial'][i][1])]==255:
                adata.obs['division'][i] = "DA"
        except:
            pass
    st.pl.space(adata, color=['division'], pointsize=0.2)