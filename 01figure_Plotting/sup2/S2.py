###Celltype
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
os.chdir("./1B/")
rawdir = "./02h5ad/"
outdir = "./1B/"
csvdir = "./01csv/"
rgb_dic = {
    'EC':[35, 139, 69],
    'Fib':[139, 74, 75],
    'Macro':[223, 101, 176],
    'mSMC':[161, 207, 250],
    'SMC':[31, 120, 180],
    'T cell':[206, 18, 86]
}
for key in [i [:-4] for i in os.listdir(csvdira)]:
    print(key)
    key1 = key[:-2]
    print(key1)
    df = pd.read_csv(csvdir+key+".csv",index_col=0)
    data = sc.read_h5ad(rawdir+key1+"_cell_gene.h5ad")
    data.obs.index=data.obs.index.astype("int64")
    data.obs=pd.merge(data.obs,df,right_index=True,left_index=True,how="left")
    data.obs=data.obs.fillna("others")
    data.obs.index=data.obs.index.astype("object")
    data=data[data.obs['Celltype'].isin(['EC','Fib','Macro','mSMC','SMC','T cell'])]
    data.obs.index=data.obs.index.astype("int64")
    seg=data.uns["seg"]
    np.unique(data.obs["Celltype"])
    height,width=seg.shape
    rgb_image=np.zeros((height,width,3),dtype=np.uint8)
    array_celltype=np.unique(data.obs['Celltype'])
    celltype_dic={}
    for celltype in array_celltype:
        print(celltype)
        celltype_dic[celltype]=data.obs_names[data.obs['Celltype']==celltype]
        i=-1
        for i in range(3):
            rgb_image[:,:,i][np.isin(seg,celltype_dic[celltype])]=rgb_dic[celltype][i]
    bgr_image=cv2.cvtColor(rgb_image,cv2.COLOR_RGB2BGR)
    cv2.imwrite(outdir+key+".tif",bgr_image)



###Cellsubtype
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
os.chdir("./03ann2_shape/")
rawdir = "./h5ad_RAw/"
outdir = "./03ann2_shape/"
csvdir = "./00csv/"
rgb_dic = {
    'Spp1+Trem2 macro.':[255,112,171],
    'SMC_like macro.':[249,228,0],
    'Inflammatory macro.':[223,101,176],
    'DCs':[183,163,227],
    'Col15a1 FBs':[194,166,140],
    'Pi16 FBs':[139,74,75],
    'ECs':[35,139,69],
    'SMC_like ECs':[6,208,1],
    'Inflammatory ECs':[0,255,156],
    'mSMCs_1':[119,67,219],
    'mSMCs_2':[161,207,250],
    'Pericytes':[51,48,228],
    'Contractile SMCs':[31,120,180],
    'Pre_modulated_SMCs':[77,119,255],
    'Tcell':[206,18,86]
}

#####let us creat processed cycle code to finish cellsegmentation
#####let us creat processed cycle code to finish cellsegmentation
for key in [i [:-4] for i in os.listdir(csvdir)]:
    print(key)
    key1 = key[:-2]
    print(key1)
    df = pd.read_csv(csvdir+key+".csv",index_col=0)
    k = df["Week"].unique()
    s = df["sample"].unique()
    data = sc.read_h5ad(rawdir+key1+"_cell_gene.h5ad")
    data.obs.index=data.obs.index.astype("int64")
    data.obs=pd.merge(data.obs,df,right_index=True,left_index=True,how="left")
    data.obs=data.obs.fillna("others")
    data.obs.index=data.obs.index.astype("object")
    data=data[data.obs['ann2'].isin(['Spp1+Trem2 macro.','Inflammatory ECs','Tcell','Inflammatory macro.','DCs','mSMCs_1','Col15a1 FBs','SMC_like macro.','mSMCs_2','Pi16 FBs','Contractile SMCs','Pre_modulated_SMCs','ECs','Pericytes','SMC_like ECs'])]
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