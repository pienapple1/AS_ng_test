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
import numpy as np

import scanpy as sc

#####set pathway
os.chdir("/shape/")
rawdir = "/h5ad_RAw/"
outdir = "/shape/"
csvdir = "/meta_stage/"

#####let us creat processed cycle code to finish cellsegmentation
for key in [i [:-10] for i in os.listdir(csvdir)]:
    print(key)
    key1 = key[:-2]
    print(key1)
    df = pd.read_csv(csvdir+key+"_stage.csv",index_col=0)
    data = sc.read_h5ad(rawdir+key1+"_cell_gene.h5ad")
    data.obs.index=data.obs.index.astype("int64")
    data.obs=pd.merge(data.obs,df,right_index=True,left_index=True,how="left")
    data.obs=data.obs.fillna("others")
    data.obs.index=data.obs.index.astype("object")
    data1=data[data.obs['satge_p_all'].isin(['stage1','stage2','stage3','stage4'])]
    data1.obs.index=data1.obs.index.astype("int64")
    ###data1.write(outdir+key + "_cell_gene.h5ad")
    print("data1 is save")
    ###复制三份数组
    red = np.zeros_like(data1.uns["seg"])
    green = np.zeros_like(data1.uns["seg"])
    blue = np.zeros_like(data1.uns["seg"])
    ###提取 每个细胞类型的obs_name，即灰度值
    stage1_list=data1.obs_names[data1.obs['satge_p_all']=="stage1"]
    stage2_list=data1.obs_names[data1.obs['satge_p_all']=="stage2"]
    stage3_list=data1.obs_names[data1.obs['satge_p_all']=="stage3"]
    stage4_list=data1.obs_names[data1.obs['satge_p_all']=="stage4"]
    ###stage1_list 175,214,250
    red[np.isin(data1.uns["seg"], stage1_list)] = 167
    green[np.isin(data1.uns["seg"], stage1_list)] = 211
    blue[np.isin(data1.uns["seg"], stage1_list)] = 212
    print("stage1_list is over")
    
    ###stage2_list 31,120,180
    red[np.isin(data1.uns["seg"], stage2_list)] = 0
    green[np.isin(data1.uns["seg"], stage2_list)] = 155
    blue[np.isin(data1.uns["seg"], stage2_list)] = 158
    print("stage2_list is over")
    
    ###stage3_list 254,181,210
    red[np.isin(data1.uns["seg"], stage3_list)] = 228
    green[np.isin(data1.uns["seg"], stage3_list)] = 193
    blue[np.isin(data1.uns["seg"], stage3_list)] = 217
    print("stage3_list is over")
    
    ###stage4_list 206,18,86
    red[np.isin(data1.uns["seg"], stage4_list)] = 199
    green[np.isin(data1.uns["seg"], stage4_list)] = 93
    blue[np.isin(data1.uns["seg"], stage4_list)] = 171
    print("stage4_list is over")
    
    
    RGB = np.dstack((blue,green,red))
    
    cv2.imwrite(key+'_celltype.tif',RGB)
    print("tif is save")