import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import skimage
import sklearn
import spateo as st
import cv2
import anndata as ad
import pandas as pd
from scipy import stats
import os


import scanpy as sc
from collections import Counter
import anndata as ad
import loompy as lp
import numpy as np

import scanpy as sc

#####set pathway
os.chdir("./01shape/")
rawdir = "./h5ad/"
outdir = "./01shape/"
csvdir = "./00csv/"

#####let us creat processed cycle code to finish cellsegmentation
for key in [i [:-4] for i in os.listdir(csvdir)]:
    print(key)
    key1 = key[:-2]
    print(key1)
    df = pd.read_csv(csvdir+key+".csv",index_col=0)
    data = sc.read_h5ad(rawdir+key1+"_cell_gene.h5ad")
    data.obs.index=data.obs.index.astype("int64")
    data.obs=pd.merge(data.obs,df,right_index=True,left_index=True,how="left")
    data.obs=data.obs.fillna("others")
    data.obs.index=data.obs.index.astype("object")
    data=data[data.obs['ann2'].isin(['DC',
    'mSMC',
    'Fib',
    'Macro',
    'Mast_cell',
    'EC_1',
    'EC_2',
    'SMC',
    'Monocyte',
    'NKT',
    'B_cell',
    'plasma cell'])]
    data.obs.index=data.obs.index.astype("int64")
    seg=data.uns["seg"]
    np.unique(data.obs["ann2"])
    height,width=seg.shape
    rgb_image=np.zeros((height,width,3),dtype=np.uint8)
    array_celltype=np.unique(data.obs['ann2'])
    celltype_dic={}
    for celltype in array_celltype:
        print(celltype)
        celltype_dic[celltype]=data.obs_names[data.obs['ann2']==celltype]
        i=-1
        for i in range(3):
            rgb_image[:,:,i][np.isin(seg,celltype_dic[celltype])]=rgb_dic[celltype][i]
    bgr_image=cv2.cvtColor(rgb_image,cv2.COLOR_RGB2BGR)
    cv2.imwrite(outdir+key+".tif",bgr_image)