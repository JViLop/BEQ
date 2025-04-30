# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 14:27:26 2024

@author: vite0005
"""

import re
import obspy as obs
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.colors import TwoSlopeNorm


nmodels = 4
nsamples = 100
eq_name = 'Tohoku' 
nstations = 216
npoints = 256
shape = (nmodels,nstations,npoints,3)
data = np.zeros(shape)

'''
parent_dir = f'//wsl.localhost/Ubuntu/root/dynamic/MudPy/{eq_name}/{nsamples}_samples/'
for i in range(1,nmodels+1):
   imodel_dir = os.path.join(parent_dir,f'{i}/output/waveforms/{eq_name}')
   for m,comp in enumerate(['N','E','Z']):
       for k in range(nstations):
           file_dir = os.path.join(imodel_dir,'%.4d.LY%s.sac'%(k,comp))
           stream = obs.read(file_dir)
           trace = stream[0].data
           data[i-1,k,:,m] = stream[0].data 
           
           
# example final time
data_Ntf = data[:,:,-1,0].T
data_Etf = data[:,:,-1,1].T
data_Ztf= data[:,:,-1,2].T

cov_Ntf = np.cov(data_Ntf)        
           
'''
'''
for t in range(npoints):
    data_Nt = data[:,:,t,0]
    data_Nt = data_Nt.T
    
    data_Et = data[:,:,t,1]
    data_Et = data_Et.T
    
    data_Zt= data[:,:,t,2]
    data_Zt = data_Zt.T

    cov_Nt = np.cov(data_Et) 
    cov_Et = np.cov(data_Et) 
    cov_Zt = np.cov(data_Zt) 
           
'''      


path = 'statics.neu'
path = '//wsl.localhost/Ubuntu/root/dynamic/MudPy/Tohoku/100_samples/3/output/waveforms/Tohoku/_summary.Tohoku.txt'
nrows,ncols = (9,24)

patch = 29

x = np.arange(patch/2,ncols*patch,patch)
y = np.arange(patch/2,nrows*patch,patch)
y = -np.flip(y)

data = np.loadtxt(path)
strike = 194*np.pi/180
strike = strike - 90.0*np.pi/180.0
n = data[:,3]
e = data[:,4]
z = data[:,5]
dy = n*np.cos(strike) + e*np.sin(strike)
dx = -n*np.sin(strike) +e*np.cos(strike)
dz = z
d = {'x':dx,'y':dy,'z':dz}


fig,axes = plt.subplots(3,1,figsize=(5,8.5),dpi=400)
for i,key in enumerate(list(d.keys())):
    parameter = d[key]
    im = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(parameter.min(),-parameter.max())),vmax=(max(-parameter.min(),parameter.max()))))
    fig.colorbar(im,ax=axes[i],shrink=0.32)
    axes[i].set_aspect('equal','box')
    axes[i].set_title(f'{key} displacement')
plt.subplots_adjust(hspace=-0.4)
fig.suptitle(f'{eq_name} FK final displacement',y=0.85,fontweight='bold')
#plt.tight_layout()
plt.savefig(f'{eq_name}_fk_final_displacement.png' )
    
        

















           