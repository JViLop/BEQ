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
import h5py

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
scales = [0.5,0.8e-1,1e-1,0.75e-1,0.75e-1]
shrinks = [0.7,0.2,0.2,0.15,0.7]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,nramps,factors,strikes,scales,shrinks):
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
        model[name]['scale'] = scales[i]
        model[name]['shrink'] = shrinks[i]
    return model



models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,ramps,factors,strikes,scales,shrinks)



# snapshots

ith = {'x':0,'y':1,'z':2}
     

nepochs = 6

nmodels = 1
nsamples = 100
npoints = 256
names = ['Iquique','Pedernales','Gorkha']
for n,eq_name in enumerate(names):
    patch = models[eq_name]['patch']
    nrows,ncols = models[eq_name]['geom']
    strike  = models[eq_name]['strike']
    scale =  models[eq_name]['scale']
    shrink = models[eq_name]['shrink']
    x = np.arange(patch/2,ncols*patch,patch)
    y = np.arange(patch/2,nrows*patch,patch)
    y = -np.flip(y,axis=0)
    
    X,Y = np.meshgrid(x,y)
    
    strike = strike*np.pi/180
    strike = strike - 90.0*np.pi/180.0
    nstations = nrows*ncols
    shape = (nmodels,nstations,npoints,3)
    data = np.zeros(shape)
    
    parent_dir = os.path.join(os.getcwd(),f'Dynamic_Simulations/{eq_name}/{nsamples}_samples')
    for i in range(1,nmodels+1):
       imodel_dir = parent_dir+f'/{i}/output/waveforms/{eq_name}'
       for l,orientation in enumerate(['N','E','Z']):
           for t in range(0,nstations):
               file_dir = imodel_dir+'/%.4d.LY%s.sac'%(t,orientation)
               stream = obs.read(file_dir)
               trace = stream[0].data
    
               data[i-1,t,:,l] = stream[0].data 
               

    N = data[:,:,:,0]  
    E = data[:,:,:,1]  
    Z = data[:,:,:,2]  
    
    #data[:,:,:,0] =  -N*np.sin(strike) + E*np.cos(strike)      
    #data[:,:,:,1] =  N*np.cos(strike) + E*np.sin(strike) 
    current_dir = os.getcwd()
    folder_dir = os.path.join(current_dir,f'Dynamic_Simulations/{eq_name}/{nsamples}_samples/H5_{nmodels}')
    os.makedirs(folder_dir, exist_ok =True)
    
    file_dir = os.path.join(folder_dir,f'{eq_name}_N_nmodels_{nmodels}.h5')
    hf = h5py.File(file_dir, 'w')
    hf.create_dataset('N', data=N)
    hf.close()
    
    file_dir = os.path.join(folder_dir,f'{eq_name}_E_nmodels_{nmodels}.h5')
    hf = h5py.File(file_dir, 'w')
    hf.create_dataset('E', data=E)
    hf.close()
    
    file_dir = os.path.join(folder_dir,f'{eq_name}_Z_nmodels_{nmodels}.h5')
    hf = h5py.File(file_dir, 'w')
    hf.create_dataset('Z', data=Z)
    hf.close()
    
    
