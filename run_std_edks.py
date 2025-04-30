# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 18:56:39 2024

@author: joanv
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

from matplotlib.colors import TwoSlopeNorm


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [0,0.5,0.5,0.5,0.5]
wspaces = [0.2,0.3,0.4,0.3,0.4]
sizes = [(14,16),(10,11),(10,12),(10,11),(14,16)]
shrinks = [0.5,0.5,0.5,0.5,0.2]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['hspace'] = hspaces[i]
        model[name]['wspace'] = wspaces[i]
        model[name]['size'] = sizes[i]
        model[name]['shrink'] = shrinks[i]
    return model

models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks)

nsamples = 100
name = str(sys.argv[1])
working_dir = os.getcwd()
edks_file = os.path.join(working_dir,f'{name}/{nsamples}_samples/EDKS/EDKS_{name}_displacement_nsamples_{nsamples}.h5')
okada_file = os.path.join(working_dir,f'{name}/{nsamples}_samples/Okada/Okada_{name}_displacement_nsamples_{nsamples}.h5')
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

def average(file_name,key):
    dset  = read_h5file(file_name,key)
    mean = np.mean(dset.transpose(),axis=1)
    
    nparameters = dset.shape[1]
    
    x = mean[:nparameters//3]
    y = mean[nparameters//3:2*nparameters//3]
    z = mean[2*nparameters//3:]

    # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
    # standard deviation (= square root of variance)

    
    return x,y,z

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
    hspace = models[name]['hspace']

    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(2,4,patch,nrows,ncols,control=1)
    nrows = len(y)
    ncols = len(x)
    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=900)
    for i,parameter_id in enumerate(std.keys()):
        parameter = std[parameter_id]

        title = f"{parameter_id} Displacement "

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
    
    plt.subplots_adjust(hspace=hspace)
    fig.suptitle(f'{name } uncertainty in {method} Surface Displacement')
    os.makedirs(os.path.join(working_dir,f'{name}/{nsamples}_samples/{method}/cov'),exist_ok=True)
    figure_dir = os.path.join(os.path.join(working_dir,f'{name}/{nsamples}_samples/{method}/cov'),f'{method}_{name}_uncertainty.png')	
    fig.savefig(figure_dir,dpi=700)
        
	
plot_cov(name,edks_file,key)
plot_cov(name,okada_file,key,method='Okada')


def plot_mean_cov(name,file_name,key,method = 'EDKS'):
    std1,std2,std3 = cov(file_name,key)
    d1,d2,d3 = average(file_name,key)
    
    d = {'Along-strike':d1,
            'Trench-normal':d2,
            'Vertical':d3}
    

    std = {'Along-strike':std1,
            'Trench-normal':std2,
            'Vertical':std3}
    
    patch = models[name]['patch']
    hspace = models[name]['hspace']
    wspace = models[name]['wspace']
    size = models[name]['size']
    shrink = models[name]['shrink']
    
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(2,4,patch,nrows,ncols,control=1)
    nrows = len(y)
    ncols = len(x)
    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, axes = plt.subplots(3,2,figsize=size,dpi=700)
    symbol_d = ['$x$ (m)','$y$ (m)','$z$ (m)']
    symbol_sigma = ['$\sigma_{x}$ (m)','$\sigma_{y}$ (m)','$\sigma_{z}$ (m)']
    for i,parameter_id in enumerate(std.keys()):
        sigma = std[parameter_id]
        val = d[parameter_id]
        title = f"{parameter_id}"

        im0 = axes[i][0].pcolormesh(x,y,val.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(val.min(),-val.max())),vmax=(max(-val.min(),val.max()))))
        im1 = axes[i][1].pcolormesh(x,y,sigma.reshape(nrows,ncols),edgecolors='k', cmap='rainbow',linewidth=0.25)

        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i][0].margins(0)
        axes[i][1].margins(0)
        fig.colorbar(im0, ax=axes[i][0],shrink=shrink,label= symbol_d[i])
        fig.colorbar(im1, ax=axes[i][1],shrink=shrink,label= symbol_sigma[i])

        axes[i][0].set_ylabel('Down-dip distance (km)',fontsize=10)
        axes[i][0].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i][0].set_aspect('equal', 'box')
        axes[i][0].set_title(title,fontweight='bold')
        axes[i][0].tick_params(labelsize=10)
        
        axes[i][1].set_ylabel('Down-dip distance (km)',fontsize=10)
        axes[i][1].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i][1].set_aspect('equal', 'box')
        axes[i][1].set_title(title,fontweight='bold')
        axes[i][1].tick_params(labelsize=10)
    
    plt.subplots_adjust(wspace=wspace,hspace=hspace)
    fig.suptitle(f'{name}',fontsize=13,fontweight='bold')
    os.makedirs(os.path.join(working_dir,f'{name}/{nsamples}_samples/{method}'),exist_ok=True)
    figure_dir = os.path.join(os.path.join(working_dir,f'{name}/{nsamples}_samples/{method}'),f'{method}_{name}_mean_and_uncertainty_.png')	
    fig.savefig(figure_dir,dpi=700)
        
	
plot_mean_cov(name,edks_file,key)
plot_mean_cov(name,okada_file,key,method='Okada')
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
