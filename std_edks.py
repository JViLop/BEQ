# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 18:56:39 2024

@author: joanv
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt



names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,0,0,0,0]


def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
    return model

models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes)

name = 'Tohoku'
edks_file = f'edks_{name}_displacement_nsamples_100.h5'
okada_file = f'okada_{name}_displacement_nsamples_100.h5'
key = 'displacement'

def read_h5file(file_name,key):
       # self.Multiple_Okada_displacement()

       f = h5py.File(file_name,'r')
       dset = np.array(f[key])
       return dset
    
def cov(file_name,key):
    dset  = read_h5file(file_name,key)
    covariance = np.cov(dset.transpose())
    
    nparameters = dset.shape[1]
    
    cov1 = covariance[:nparameters//3,:nparameters//3]
    cov2 = covariance[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    cov3 = covariance[2*nparameters//3:,2*nparameters//3:]

    # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
    # standard deviation (= square root of variance)

    std1 = np.sqrt(cov1.diagonal())
    std2 = np.sqrt(cov2.diagonal())
    std3 = np.sqrt(cov3.diagonal())
    
    return std1,std2,std3



def corr(file_name,key):
    dset  = read_h5file()
    correlation = np.corrcoef(dset.transpose())
    
    nparameters = dset.shape[1]
    corr1 = correlation[:nparameters//3,:nparameters//3]
    corr2 = correlation[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    corr3 = correlation[2*nparameters//3:,2*nparameters//3:]

    return corr1, corr2, corr3

def set_stn(c,factor,patch,nrows,ncols,control=0):
      
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn
  
def plot_cov(name,file_name,key,method = 'EDKS'):
    std1,std2,std3 = cov(file_name,key)


    std = {'Along-strike':std1,
            'Trench-normal':std2,
            'Vertical':std3}
    
    patch = models[name]['patch']
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(1,4,patch,nrows,ncols,control=1)
    nrows = len(y)
    ncols = len(x)
    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=900)
    for i,parameter_id in enumerate(std.keys()):
        parameter = std[parameter_id]

        title = f"{name} Uncertainty in {method} Surface Displacement "

        im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow',linewidth=0.25)


        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i].margins(0)        
        fig.colorbar(im0, ax=axes[i],shrink=0.65,label='Uncertainty (m)')
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=12)
        

plot_cov(name,edks_file,key)
plot_cov(name,okada_file,key,method='Okada')
# def cov(self):
#     dset  = self.read_h5file()
#     self.covariance = np.cov(dset.transpose())
    
#     nparameters = dset.shape[1]
    
#     self.cov1 = self.covariance[:nparameters//3,:nparameters//3]
#     self.cov2 = self.covariance[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
#     self.cov3 = self.covariance[2*nparameters//3:,2*nparameters//3:]

#     # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
#     # standard deviation (= square root of variance)

#     self.std1 = np.sqrt(self.cov1.diagonal())
#     self.std2 = np.sqrt(self.cov2.diagonal())
#     self.std3 = np.sqrt(self.cov3.diagonal())