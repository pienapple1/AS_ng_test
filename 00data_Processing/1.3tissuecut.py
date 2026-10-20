import cv2
import numpy as np
import math
import argparse
import os

path1 = "/mask/"
path2 = "/tissuecut_csv/"
xoffset=0
yoffset=0
bin=50
threshold=20
for key in [i [:-4] for i in os.listdir(path1)]:
    print(key)
    input = path1+key+".tif"
    output = path2+key+"/"
    if not os.path.exists(output):
        os.makedirs(output)
    bin_size=int(bin)
    os.chdir(output)
    img = cv2.imread(input,cv2.IMREAD_GRAYSCALE)
    raw_img = img.T
    cv2.imwrite('./1.grey.png', raw_img)
    offset=raw_img[int(xoffset):,int(yoffset):]
    cv2.imwrite('./2.offset.png', offset)
    ret, binary = cv2.threshold(offset, int(threshold), 255, cv2.THRESH_BINARY)
    cv2.imwrite('./3.binary.png', binary)
    kernel = np.ones((100, 100), np.uint8)
    dilated_img = cv2.dilate(binary, kernel)
    cv2.imwrite('./4.dilated.png', dilated_img)
    eroded_img = cv2.erode(dilated_img, kernel)
    cv2.imwrite('./5.eroded.png', eroded_img)
    contours, hierarchy = cv2.findContours(eroded_img, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    areas = []
    for i in contours:
        cur_area = cv2.contourArea(i)
        areas.append(cur_area)
        
    idx = np.argsort(areas)[::-1][:5]
    remain_contours = []
    for i in idx:
        remain_contours.append(contours[i])
    
    blank_img = np.uint8(np.full(eroded_img.shape, 0))
    rst = cv2.drawContours(blank_img.copy(),remain_contours,-1,(255,0,0),-1)
    cv2.imwrite('./6.rst.png', rst)
    
    rst_len_0=math.ceil(rst.shape[0]/bin_size)
    rst_len_1=math.ceil(rst.shape[1]/bin_size)
    pool=np.zeros([rst_len_0,rst_len_1])
    
    for i in np.arange(rst_len_0):
        i_min=i*bin_size
        i_max=(i+1)*bin_size
        
        for j in np.arange(rst_len_1):
            j_min=j*bin_size

            j_max=(j+1)*bin_size
            pool[i,j]=rst[i_min:i_max,j_min:j_max].mean()
            
    cv2.imwrite('./7.pooling.png', pool)
    ret, binary2 = cv2.threshold(pool, 255/2, 255, cv2.THRESH_BINARY)
    binary2=binary2.astype(int)
    cv2.imwrite('./8.binary2.png', binary2)
    
    bin_list=list()
    for i in np.arange(rst_len_0):
        for j in np.arange(rst_len_1):
            if binary2[i,j]==255:
                bin_list.append(np.array([str(i+1),str(j+1),str(i+1)+"_"+str(j+1)]))
                
    bin_array=np.array(bin_list)
    np.savetxt("tissue.txt",bin_array,fmt="%s")