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
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
#hspaces = [0,0.5,0.5,0.5,0.5]
#wspaces = [0.2,0.3,0.4,0.3,0.25]
sizes = [(12,9),(10,10),(10,10),(10,10),(12,9)]
#shrinks = [0.5,0.7,0.70,0.70,0.75]

hspaces = [-0.1,0,-0.1,-0.1,0.05]
wspaces = [0.30,0.5,0.4,0.5,0.3]
sizes = [(14,6),(10,8),(10,6.0),(10,8),(14,7)]
shrinks = [0.5,0.6,0.6,0.6,0.75]
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


name = str(sys.argv[1])
nsamples = str(sys.argv[2])
working_dir = os.getcwd()

mean_file = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')
okada_file = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/Stress change/{nsamples}_samples/data_ensemble/{name}_Stress change_nsamples_{nsamples}.h5')
key = 'Stress change'


def get_U_hypo(file_name,shape):
       # self.Multiple_Okada_displacement()
       df = pd.read_csv(file_name)
       df = df.drop(df.columns[0],axis=1)
       UPERP = np.flip(df['U_perp'].values.reshape(shape[0],shape[1],order='F'),axis=0)
       UPARALLEL = np.flip(df['U_parallel'].values.reshape(shape[0],shape[1],order='F'),axis=0)
       U = np.sqrt(UPERP**2  + UPARALLEL**2)
       hyp_as = df['Hypo_as'].values[0]
       hyp_dip = -df['Hypo_dd'].values[0]
       return U, hyp_as, hyp_dip



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
            'Along-dip':std2,
            'Normal':std3}
    
    patch = models[name]['patch']
    hspace = models[name]['hspace']

    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(1,1,patch,nrows,ncols,control=0)
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
        
	
#plot_cov(name,edks_file,key)
#plot_cov(name,okada_file,key,method='Okada')


def plot_mean_cov(name,file_name,file_name_mean,key,method = 'EDKS'):
    std1,std2,std3 = cov(file_name,key)
    d1,d2,d3 = average(file_name,key)
    
    d = {'Along-strike':d2,
            'Along-dip':d3,
            'Normal':d1}
    

    std = {'Along-strike':std2,
            'Along-dip':std3,
            'Normal':std1}
    
    rel_std = {'Along-strike':100*std2/d2,
            'Along-dip':100*std3/d3,
            'Normal':100*std1/d1}
    
    patch = models[name]['patch']
    hspace = models[name]['hspace']
    wspace = models[name]['wspace']
    size = models[name]['size']
    shrink = models[name]['shrink']
    
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    

    U,hypo_as,hypo_dd = get_U_hypo(file_name_mean,(nrows,ncols))
   


    x,y = set_stn(1,1,patch,nrows,ncols,control=0)
    nrows = len(y)
    ncols = len(x)
    X, Y = np.meshgrid(x,y)

    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, axes = plt.subplots(3,3,figsize=size,dpi=700,sharex=True)
    symbol_d = ['$\Delta\sigma_{\perp}$ (MPa)','$\Delta\sigma_{\parallel}$ (MPa)','$\Delta\sigma_{n}$ (MPa)']
    symbol_sigma = ['$s$ (MPa)','$s$ (MPa)','$s$ (MPa)']
    for i,parameter_id in enumerate(std.keys()):
        sigma = std[parameter_id]
        val = d[parameter_id]
        rel = rel_std[parameter_id]
        title = f"{parameter_id}"

        im0 = axes[i][0].pcolormesh(x,y,val.reshape(nrows,ncols),cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(val.min(),-val.max())),vmax=(max(-val.min(),val.max()))))
        im1 = axes[i][1].pcolormesh(x,y,sigma.reshape(nrows,ncols), cmap='rainbow')
        im2 = axes[i][2].pcolormesh(x,y,rel.reshape(nrows,ncols), cmap='Blues',vmin=0,vmax=100)

        #im0 = axes[i][0].pcolormesh(x,y,val.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(val.min(),-val.max())),vmax=(max(-val.min(),val.max()))))
        #im1 = axes[i][1].pcolormesh(x,y,sigma.reshape(nrows,ncols),edgecolors='k', cmap='rainbow',linewidth=0.25)
        cs0 = axes[i][0].contour(X,Y,U.reshape(nrows,ncols),levels=7,linewidths=0.4,colors='black',alpha=0.75,linestyles='dashed')
        cs1 = axes[i][1].contour(X,Y,U.reshape(nrows,ncols),levels=7,linewidths=0.4,colors='black',alpha=0.75,linestyles='dashed')
        cs2 = axes[i][2].contour(X,Y,U.reshape(nrows,ncols),levels=7,linewidths=0.4,colors='black',alpha=0.75,linestyles='dashed')
        axes[i][0].scatter(hypo_as,hypo_dd,marker='*',s=125,color='yellow',edgecolor='black',linewidth=0.5)
        axes[i][1].scatter(hypo_as,hypo_dd,marker='*',s=125,color='yellow',edgecolor='black',linewidth=0.5)
        axes[i][2].scatter(hypo_as,hypo_dd,marker='*',s=125,color='yellow',edgecolor='black',linewidth=0.5)
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i][0].margins(0)
        axes[i][1].margins(0)
        axes[i][2].margins(0)
        fig.colorbar(im0, ax=axes[i][0],shrink=shrink,label= symbol_d[i])
        fig.colorbar(im1, ax=axes[i][1],shrink=shrink,label= symbol_sigma[i])
        fig.colorbar(im2, ax=axes[i][2],shrink=shrink,label= 'rel. error (%)',extend='max')


        axes[1][0].set_ylabel('Down-dip distance (km)',fontsize=9)
        axes[2][i].set_xlabel('Along-strike distance (km)',fontsize=9)
        axes[i][0].set_aspect('equal', 'box')
        axes[i][0].set_title(title,fontweight='bold')
        axes[i][0].tick_params(labelsize=9)
        axes[0][i].text(-0.1,1.15,'abc'[i],fontweight='bold',fontsize=20,transform=axes[0][i].transAxes)
        
        # axes[i][1].set_ylabel('Trench-normal distance (km)',fontsize=10)
        # axes[i][1].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i][1].set_aspect('equal', 'box')
        axes[i][1].set_title(title,fontweight='bold')
        axes[i][1].tick_params(labelsize=9)
        
        # axes[i][2].set_ylabel('Trench-normal distance (km)',fontsize=10)
        # axes[i][2].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i][2].set_aspect('equal', 'box')
        axes[i][2].set_title(title,fontweight='bold')
        axes[i][2].tick_params(labelsize=9)
        
    
    plt.subplots_adjust(wspace=wspace,hspace=hspace)
    #fig.suptitle(f'{name}',fontsize=15,fontweight='bold')
    figure_dir = os.path.join(os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/Stress change/{nsamples}_samples/figures'),f'{method}_{name}_mean_and_uncertainty_n.pdf')	
    fig.savefig(figure_dir)
        
	
#plot_mean_cov(name,edks_file,key)
plot_mean_cov(name,okada_file,mean_file,key,method='Okada')

