# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 22:44:34 2024

@author: joanv
"""

import numpy as np
import pandas as pd

df = pd.read_csv('C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/INPUT/Tohoku/model/kinematic/100_samples/mean/Tohoku_mean_kinematic_model.csv')

nrows, ncols = 9,24
shape_geom = (nrows,ncols)
patch = 29


xsrc = np.arange((1/2)*patch,ncols*patch, patch)

yref = np.arange(-(nrows-1/2)*patch,0,patch)

dip = df['dip'].values[:nrows]
proj_dysrc = -patch*np.cos(dip*np.pi/180)
proj_ysrc = np.zeros_like(proj_dysrc)
for i in range(len(proj_ysrc)):
    proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 

ysrc = np.flip(proj_ysrc)

Xsrc,Ysrc =  np.meshgrid(xsrc,ysrc)



