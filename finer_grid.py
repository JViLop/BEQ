# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:07:30 2024

@author: joanv
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def arrays1(c):
    ncols = 24
    nrows = 9
    patch = 29
    dpatch = patch/c
    factor = 4
    control = 0
    
    offsetx =  (ncols//factor)*patch
    offsety =  (nrows//factor)*patch
    xstn = np.arange(-control*offsetx + dpatch, c*ncols*dpatch + control*offsetx,dpatch)
    ystn = np.arange(-control*offsety - c*(nrows-1/2)*dpatch, control*offsety,dpatch)
    return xstn, ystn

def arrays2(c):
    ncols = 24
    nrows = 9
    patch = 29
    dpatch = patch/c
    factor = nrows
    control = 1
    
    offsetx =  (ncols//factor)*patch
    offsety =  (nrows//factor)*patch
    xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
    ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
    ystn = np.flip(ystn)
    return xstn, ystn



ncols = 24
nrows = 9
patch = 29
dpatch = patch
factor = 4

offsetx =  (ncols//factor)*patch
offsety =  (nrows//factor)*patch
xstn0 = np.arange(patch/2 , ncols*patch,dpatch)
ystn0 = np.arange(- (nrows-1/2)*patch, 0,dpatch)
nx0,ny0 = len(xstn0),len(ystn0)
xstn1,ystn1 = arrays2(3
                      )

nx1,ny1 = len(xstn1),len(ystn1)
print(nx0,ny0)
print(nx1,ny1)
fig,ax = plt.subplots()
im =ax.pcolormesh(xstn0,ystn0,np.ones((ny0,nx0)),edgecolors='yellow',cmap='viridis')
im2 =ax.pcolormesh(xstn1,ystn1,np.ones((ny1,nx1)),edgecolors='white',alpha=0.2,cmap='viridis')
ax.set_aspect('equal', 'box')


step = int(len(ystn1)/len(ystn0))

zstn = np.flip(pd.read_csv('Tohoku_mean_kinematic_model.csv')['depth'].values)

unique_val = np.flip(np.unique(zstn))
depth = np.ones((len(ystn1),len(xstn1)))

# for i in range(0,len(ystn1),int(step)):
#     print(i,i+step)
#     depth[i:i+step,:] = unique_val[int(i/step)]*depth[i:i+step,:] 


ncols = 24
nrows = 9
patch = 29*1e3
factor = 4
control = 0
def set_stn(c,control=0):
    
    dpatch = patch/c
    offsetx =  (ncols//factor)*patch
    offsety =  (nrows//factor)*patch
    xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
    ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
    ystn = np.flip(ystn)
    return xstn, ystn

xstn1,ystn1 = set_stn(1)
xstn2,ystn2 = set_stn(2,control=1)