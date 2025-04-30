# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 12:14:57 2024

@author: joanv
"""

import re
import numpy as np
import pandas as pd
import h5py as h5

def array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()
f = open('model_test', "r")
content = f.readlines()

model = content[12:]
keys = re.findall('\w+',content[11],flags=re.DOTALL)
fmt = re.compile('\d+.\d+')
layers = [list(map(float,re.findall(fmt,m))) for m in model]
layers = np.array(layers).T
model_dict = dict(zip(keys,layers))

    
def get_mu(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VS'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        return mu
    
df = pd.read_csv('Tohoku_mean_kinematic_model.csv')

nrows,ncols = 9,24
shape = (nrows,ncols)
patchsize = 29



Uperp = array_formatter(df['U_perp'].values,shape)
Uparallel =  array_formatter(df['U_parallel'].values,shape)
SLIP = array_formatter(df['Slip'].values,shape)
RAKE = array_formatter(np.arctan2(Uparallel,Uperp)*(180/np.pi),shape)
DIP = array_formatter(df['dip'].values,shape)
DEPTH = array_formatter(df['depth'].values,shape)
DURATION = array_formatter(df['Tr'].values,shape)

M0_total = 0
for i in range(nrows*ncols):
    mu = get_mu(model_dict,DEPTH[i])
    A = (patchsize*1e3)**2
    M0 = A*mu*SLIP[i]
    print(f'depth = {round(DEPTH[i],3)}, mu = {round(mu*1e-9,3)}, M0 = {round(M0,3)}')
    M0_total += M0
    

Mw = (2/3)*(np.log10(M0_total) - 9.1)
    
    