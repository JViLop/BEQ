# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 22:51:31 2024

@author: joanv
"""



import pandas as pd
import numpy as np
import os

from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import h5py as h5
import io
import os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

working_dir = os.getcwd()

def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)
    
def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()
    
def set_stn(c,factor,patch,control=0):
      
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn
  
def proj_ysrc_coords(patch,dip):
      proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
      proj_ysrc = np.zeros_like(proj_dysrc)
      for i in range(len(proj_ysrc)):
          proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
      ysrc = np.flip(proj_ysrc)
        
      return ysrc


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,0,0,0,0]
z_offsets = [2.8,0,0,0,0]
nramps = [0,0,0,0,0]
def build_model(names,geoms,patches,arrow_sizes,nparams,rakes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['z_offset'] = z_offsets[i]
        model[name]['nramp'] = nramps[i]
    return model


def get_mu(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VP'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        return mu

models = build_model(names,geoms,patches,arrow_sizes,nparams,rakes)

km = 1000

name = 'Tohoku'
model_type = 'kinematic'
nsamples = 100



nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
nparam = models[name]['nparam']

shape = (nrows,ncols)
patch = models[name]['patch']
nramp = models[name]['nramp']




    
df = pd.read_csv(f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')




Uperp = two_array_formatter(df['U_perp'].values,shape)
Uparallel =  two_array_formatter(df['U_parallel'].values,shape)


U = np.sqrt(Uperp**2+Uparallel**2)

fig,ax = plt.subplots()
th = 10
ths = [th]
xsrc = np.arange(patch/2 , ncols*patch,patch)
ysrc = np.arange(-(nrows-1/2)*patch,0,patch)
Xsrc ,Ysrc = np.meshgrid(xsrc,ysrc)
im= ax.pcolormesh(xsrc,ysrc,U, edgecolors='k',cmap='YlOrRd')
cs =  ax.contour(Xsrc,Ysrc,U,levels=np.array(ths)*np.max(U)/100,linewidths=[1,1.5,2.0],colors='black',linestyles='dashed')
ax.set_aspect('equal','box')

ax.set_xlabel('Along-strike distance (km)')
fig.colorbar(im,ax=ax,shrink=0.5,label='U(m)')
ax.set_ylabel('Along-dip distance (km)')
ax.set_title('Tohoku',fontweight='bold')
fig.savefig(f'contour_{th}.png',dpi=700)