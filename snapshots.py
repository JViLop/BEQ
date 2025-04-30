# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 09:22:37 2024

@author: vite0005
"""

import re
import obspy as obs
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
import sys

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']

nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patches = [29,17,18,10,15]
geoms = [(nrows[i],ncols[i]) for i in range(len(names))]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,90,90,107,360-99]
ramps = [0,3,0,0,9]
factors =[3,3,3,2,2]
strikes = [194,-13.58,4,293,27.05]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,nramps,factors,strikes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i]
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i]
        model[name]['nramp'] = nramps[i]
        model[name]['factor'] = factors[i]
        model[name]['factor'] = factors[i]
        model[name]['strike'] = strikes[i]
    return model



models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,ramps,factors,strikes)



# snapshots

ith = {'x':0,'y':1,'z':2}
     

nepochs = 6

for comp in ['x','y','z']:
    fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20),dpi=500)
    
    for n,eq_name in enumerate(names):
        patch = models[eq_name]['patch']
        nrows,ncols = models[eq_name]['geom']
        strike  = models[eq_name]['strike']
        
        x = np.arange(patch/2,ncols*patch,patch)
        y = np.arange(patch/2,nrows*patch,patch)
        y = -np.flip(y,axis=0)
        
        
        strike = strike*np.pi/180
        strike = strike - 90.0*np.pi/180.0
        nmodels = 10
        nsamples = 100
        
        nstations = nrows*ncols
        npoints = 256
        shape = (nmodels,nstations,npoints,3)
        data = np.zeros(shape)
        
        parent_dir = '/home/josevilo/Dynamic/MudPy/examples/%s/%s_samples'%(eq_name,nsamples)
        for i in range(1,nmodels+1):
           imodel_dir = parent_dir+f'/{i}/output/waveforms/{eq_name}'
           for l,orientation in enumerate(['N','E','Z']):
               for k0 in range(nstations):
                   file_dir =imodel_dir+'/%.4d.LY%s.sac'%(k0,orientation)
                   stream = obs.read(file_dir)
             
                   trace = stream[0].data
        
                   data[i-1,k0,:,l] = stream[0].data 
                   
    
        N = data[:,:,:,0]  
        E = data[:,:,:,1]  
        
        data[:,:,:,0] =  -N*np.sin(strike) + E*np.cos(strike)      
        data[:,:,:,1] =  N*np.cos(strike) + E*np.sin(strike) 
        
        DATA = np.mean(data,axis=0)
        for k,epoch in enumerate(list(np.logspace(12,np.log10(256),num=nepochs,base=10))):
                 
            m = ith[comp]
            if epoch==256:
                epoch-=1
            G0 = DATA[:,int(epoch),m] 
            
            cax1 = axes[n][k].pcolormesh(x, y, G0.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,vmin = min(min(G0),-max(G0)),vmax = max(-min(G0),max(G0)),cmap='bwr')
            if n==0:
                epoch+=1
                axes[n][k].set_title(f'{eq_name}' + ' $t = %s s$'%(round(epoch,0)))
            else:
                
                axes[n][0].set_title(f'{eq_name}')
            
            axes[n][k].set_aspect('equal','box')
            fig.colorbar(cax1,shrink=0.2,label='$d_{%s} (m)$'%(comp))
            fig.tight_layout()
    fig.savefig(f'snapshot_all_{comp}.png')
