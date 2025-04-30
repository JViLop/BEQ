# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 17:09:08 2024

@author: joanv
"""
import re
import obspy as obs
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from matplotlib.colors import TwoSlopeNorm
import h5py


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [0.8,0.8,0.8,0.8,0.8]
wspaces = [0.2,0.3,0.4,0.3,0.2]
sizes = [(14,16),(10,11),(10,12),(10,11),(14,16)]
shrinks = [0.7,0.7,0.8,0.7,0.7]
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



Gorkha_file = 'Gorkha_Stress change_nsamples_100.h5'
Tohoku_file =  'Tohoku_Stress change_nsamples_100.h5'
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
  

        
	
def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()

def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)




def plot_coeff():
    key = 'Stress change'
    fig, axes = plt.subplots(5,1,figsize=(5,10),dpi=700)
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)

        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = one_array_formatter(df['U_perp'].values,shape)
        Uparallel =  one_array_formatter(df['U_parallel'].values,shape)
        
        U = np.sqrt(Uperp**2 + Uparallel**2)
        Tr = one_array_formatter(df['Tr'].values,shape)
        file_name =  f'{name}_Stress change_nsamples_100.h5'
        std1,std2,std3 = cov(file_name,key)
        d1,d2,d3 = average(file_name,key)
        
        d = {'Along-strike':d2,
                'Along-dip':d3,
                'Normal':d1}
        shear = np.sqrt((d['Along-strike'])**2 + (d['Along-dip'])**2)
        normal = abs(d['Normal'])
        mu = shear/normal

        
        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
        X,Y = np.meshgrid(x,y)
        nrows = len(y)
        ncols = len(x)
        # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"

        im0 = axes[i].pcolormesh(x,y,mu.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='rainbow',vmin = 0,vmax=10)
        axes[i].contour(X,Y,U.reshape(nrows,ncols),levels=6,colors='k',linewidths=0.5,linestyles='dashed')
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i].margins(0)
        fig.colorbar(im0, ax=axes[i],shrink=shrink,label= '$\mu$',extend='max')

        axes[i].set_ylabel('Down-dip distance (km)',fontsize=8)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=8)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold',fontsize=9)
        axes[i].tick_params(labelsize=8)
        #fig.align_ylabels(fig.get_axes())
        plt.subplots_adjust(hspace=hspace)
    figure_dir = 'friction_coeff.png'
    fig.savefig(figure_dir,dpi=700)
    
plot_coeff() 

def plot_static():
    key = 'Stress change'
    fig, axes = plt.subplots(5,1,figsize=(5,12),dpi=700)
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)

        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = one_array_formatter(df['U_perp'].values,shape)
        Uparallel =  one_array_formatter(df['U_parallel'].values,shape)
        
        U = np.sqrt(Uperp**2 + Uparallel**2)
        Tr = one_array_formatter(df['Tr'].values,shape)
        file_name =  f'{name}_Stress change_nsamples_100.h5'
        d1,d2,d3 = average(file_name,key)
        patch = models[name]['patch']
        d = {'Along-strike':d2,
                'Along-dip':d3,
                'Normal':d1}	
        
        stress_drop = (1e-6)*(30e9*U)/(patch*1e3)
        
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
        X,Y = np.meshgrid(x,y)
        nrows = len(y)
        ncols = len(x)
        im0 = axes[i].pcolormesh(x,y,stress_drop.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='rainbow')
        axes[i].contour(X,Y,U.reshape(nrows,ncols),levels=6,colors='k',linewidths=0.5,linestyles='dashed')
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        
        axes[i].margins(0)
        fig.colorbar(im0, ax=axes[i],shrink=shrink,label='$\Delta\sigma$ (MPa)')
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"

        axes[i].set_ylabel('Down-dip distance (km)',fontsize=10)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        plt.subplots_adjust(hspace=hspace)
        axes[i].tick_params(labelsize=10)
    figure_dir = 'static_stress.png'
    fig.savefig(figure_dir,dpi=700)
plot_static()



def plot_dynamic(vel = 1.5):
    key = 'Stress change'
    fig, axes = plt.subplots(5,1,figsize=(5,12),dpi=700)
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)

        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = one_array_formatter(df['U_perp'].values,shape)
        Uparallel =  one_array_formatter(df['U_parallel'].values,shape)
        
        U = np.sqrt(Uperp**2 + Uparallel**2)
        Tr = one_array_formatter(df['Tr'].values,shape)
        file_name =  f'{name}_Stress change_nsamples_100.h5'
        d1,d2,d3 = average(file_name,key)
        patch = models[name]['patch']
        d = {'Along-strike':d2,
                'Along-dip':d3,
                'Normal':d1}	
        
        dyn_stress_drop = (2*1e-6)*(30e9*U)/(2*Tr*vel*1e3)
        
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
        X,Y = np.meshgrid(x,y)
        nrows = len(y)
        ncols = len(x)
        im0 = axes[i].pcolormesh(x,y,dyn_stress_drop.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='rainbow')
        axes[i].contour(X,Y,U.reshape(nrows,ncols),levels=6,colors='k',linewidths=0.5,linestyles='dashed')
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        
        axes[i].margins(0)
        fig.colorbar(im0, ax=axes[i],shrink=shrink,label='$\Delta\sigma_d$ (MPa)')
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"

        axes[i].set_ylabel('Down-dip distance (km)',fontsize=10)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        plt.subplots_adjust(hspace=hspace)
        axes[i].tick_params(labelsize=10)
    figure_dir = 'dynamic_stress.png'
    fig.savefig(figure_dir,dpi=700)
plot_dynamic()


def plot_relation():
    key = 'Stress change'
    fig, ax = plt.subplots(dpi=700)
    markers = ['s','^','v','*','o']
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)

        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = one_array_formatter(df['U_perp'].values,shape)
        Uparallel =  one_array_formatter(df['U_parallel'].values,shape)
        
        U = np.sqrt(Uperp**2 + Uparallel**2)

        file_name =  f'{name}_Stress change_nsamples_100.h5'
        d1,d2,d3 = average(file_name,key)
        patch = models[name]['patch']

        M0 = (np.sum(U)*30e9*(patch*1e3)**2)
        S = np.ones_like(M0)*ncols*nrows*(patch*1e3)**2
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"

        ax.scatter(M0*(2/3),S,marker=markers[i],s=65,ec='black',label=f'{title}')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.legend()
        title = 'Area-Moment relation'
        ax.set_xlabel('Scalar Moment $(M_0)$',fontsize=10)
        ax.set_ylabel('Area $(S)$',fontsize=10)
        ax.set_title(title,fontweight='bold')
        ax.tick_params(labelsize=10)
    figure_dir = 'relation_stress.png'
    fig.savefig(figure_dir,dpi=700)
plot_relation()



def plot_along_strike_U():
    key = 'Stress change'
    fig, axes = plt.subplots(5,1,figsize=(5,12),dpi=700)
    markers = ['s','^','v','*','o']
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)

        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = two_array_formatter(df['U_perp'].values,shape)
        Uparallel =  two_array_formatter(df['U_parallel'].values,shape)
        
        U = np.sqrt(Uperp**2 + Uparallel**2)
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
       

        file_name =  f'{name}_Stress change_nsamples_100.h5'
        d1,d2,d3 = average(file_name,key)
        patch = models[name]['patch']

        M0 = (np.sum(U)*30e9*(patch*1e3)**2)
        S = np.ones_like(M0)*ncols*nrows*(patch*1e3)**2
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"
        for j in range(nrows):
            U_strike = U[j,:]
            axes[i].plot(x,U_strike,marker = '.',ms=1.25,lw=0.5,color='blue')
         
        axes[i].set_xlabel('Distance along strike $(km)$',fontsize=10)
        axes[i].set_ylabel('Slip $(m)$',fontsize=10)
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=10)
        plt.subplots_adjust(hspace=hspace)
    figure_dir = 'slip_strike.png'
    fig.savefig(figure_dir,dpi=700)
plot_along_strike_U()

def plot_along_strike_rate():
    key = 'Stress change'
    fig, axes = plt.subplots(5,1,figsize=(5,12),dpi=700)
    markers = ['s','^','v','*','o']
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)

        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = two_array_formatter(df['U_perp'].values,shape)
        Uparallel =  two_array_formatter(df['U_parallel'].values,shape)
        Tr = two_array_formatter(df['Tr'].values,shape)
        U = np.sqrt(Uperp**2 + Uparallel**2)
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
        
        slip_rate = U/Tr
        

        file_name =  f'{name}_Stress change_nsamples_100.h5'
        d1,d2,d3 = average(file_name,key)
        patch = models[name]['patch']

        M0 = (np.sum(U)*30e9*(patch*1e3)**2)
        S = np.ones_like(M0)*ncols*nrows*(patch*1e3)**2
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"
        for j in range(nrows):
            rate_strike = slip_rate[j,:]
            axes[i].plot(x,rate_strike,marker = '.',ms=1.25,lw=0.5,color='blue')
        axes[i].set_xlabel('Distance along strike $(km)$',fontsize=10)
        axes[i].set_ylabel('Slip rate $(m/s)$',fontsize=10)
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=10)
        plt.subplots_adjust(hspace=hspace)
    figure_dir = 'slip_rate_strike.png'
    fig.savefig(figure_dir,dpi=700)
plot_along_strike_rate()




def plot_diff_stress(pf,mu,rho):
    g = 9.8
    key = 'Stress change'
    fig, axes = plt.subplots(5,2,figsize=(10,10),dpi=700)
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)
        print(df)
        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = one_array_formatter(df['U_perp'].values,shape)
        Uparallel =  one_array_formatter(df['U_parallel'].values,shape)
        depth  = one_array_formatter(df['depth'].values,shape)
        U = np.sqrt(Uperp**2 + Uparallel**2)
        Tr = one_array_formatter(df['Tr'].values,shape)
        file_name =  f'{name}_Stress change_nsamples_100.h5'
        std1,std2,std3 = cov(file_name,key)
        d1,d2,d3 = average(file_name,key)
        
        rhs = 1e-6*(((1/(np.sqrt(1+mu**2) - mu)**2)-1)*((depth*1e3)*rho*g - pf))
        d = {'Along-strike':d2,
                'Along-dip':d3,
                'Normal':d1}

        diff = np.abs(d['Along-dip']) - np.abs(d['Normal'])
      
        print(f'{name}',rhs)
    
        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
        X,Y = np.meshgrid(x,y)
        nrows = len(y)
        ncols = len(x)
        # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"

        im0 = axes[i][0].pcolormesh(x,y,diff.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='rainbow')
        axes[i][0].contour(X,Y,U.reshape(nrows,ncols),levels=6,colors='k',linewidths=0.5,linestyles='dashed')
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i][0].margins(0)
        fig.colorbar(im0, ax=axes[i][0],shrink=shrink,label= '$\sigma_{1} - \sigma_{3}$ (MPa)')

        axes[i][0].set_ylabel('Down-dip distance (km)',fontsize=8)
        axes[i][0].set_xlabel('Along-strike distance (km)',fontsize=8)
        axes[i][0].set_aspect('equal', 'box')
        axes[i][0].set_title(title,fontweight='bold',fontsize=9)
        axes[i][0].tick_params(labelsize=8)
        
        im0 = axes[i][1].pcolormesh(x,y,rhs.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='rainbow')
        axes[i][1].contour(X,Y,U.reshape(nrows,ncols),levels=6,colors='k',linewidths=0.5,linestyles='dashed')
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i][1].margins(0)
        fig.colorbar(im0, ax=axes[i][1],shrink=shrink,label= 'Sibson Criteria (MPa)')

        axes[i][1].set_ylabel('Down-dip distance (km)',fontsize=8)
        axes[i][1].set_xlabel('Along-strike distance (km)',fontsize=8)
        axes[i][1].set_aspect('equal', 'box')
        axes[i][1].set_title(title,fontweight='bold',fontsize=9)
        axes[i][1].tick_params(labelsize=8)
        #fig.align_ylabels(fig.get_axes())
        plt.subplots_adjust(hspace=hspace)
    figure_dir = 'slip_criteria.png'
    fig.savefig(figure_dir,dpi=700)
    
plot_diff_stress(10*1e6,0.6,3000) 



def coulomb_stress(pf,mu,rho):
    g = 9.8
    key = 'Stress change'
    fig, axes = plt.subplots(5,1,figsize=(10,10),dpi=700)
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        mean_dir = f'INPUT/{name}/model/kinematic/100_samples/mean'
        mean_file = os.path.join(mean_dir,f'{name}_mean_kinematic_model.csv')
        
        df = pd.read_csv(mean_file)
        print(df)
        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        shrink = models[name]['shrink']  
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        shape = (nrows,ncols)      
        Uperp = one_array_formatter(df['U_perp'].values,shape)
        Uparallel =  one_array_formatter(df['U_parallel'].values,shape)
        depth  = one_array_formatter(df['depth'].values,shape)
        U = np.sqrt(Uperp**2 + Uparallel**2)
        Tr = one_array_formatter(df['Tr'].values,shape)
        file_name =  f'{name}_Stress change_nsamples_100.h5'
        std1,std2,std3 = cov(file_name,key)
        d1,d2,d3 = average(file_name,key)
        
        d = {'Along-strike':d2,
                'Along-dip':d3,
                'Normal':d1}

        coulomb =  d['Along-dip']+ mu*(d['Normal'])
      

    
        patch = models[name]['patch']
        hspace = models[name]['hspace']
        wspace = models[name]['wspace']
        shrink = models[name]['shrink']
        
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        
        x,y = set_stn(1,4,patch,nrows,ncols,control=0)
        X,Y = np.meshgrid(x,y)
        nrows = len(y)
        ncols = len(x)
        # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        mag = {'Tohoku':9.1,'Illapel':8.3,'Iquique':8.1,'Pedernales':7.8,'Gorkha':7.8}
        title = "$M_{w}%s$"%(mag[name])+f" {name}"

        im0 = axes[i].pcolormesh(x,y,coulomb.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,cmap='bwr',vmin=min(np.min(coulomb),-np.max(coulomb)),vmax=max(-np.min(coulomb),np.max(coulomb)))
        axes[i].contour(X,Y,U.reshape(nrows,ncols),levels=6,colors='k',linewidths=0.5,linestyles='dashed')
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i].margins(0)
        fig.colorbar(im0, ax=axes[i],shrink=shrink,label= '$\Delta CFS$ (MPa)')

        axes[i].set_ylabel('Down-dip distance (km)',fontsize=8)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=8)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold',fontsize=9)
        axes[i].tick_params(labelsize=8)
        
        
        plt.subplots_adjust(hspace=hspace)
    figure_dir = 'coulomb.png'
    fig.savefig(figure_dir,dpi=700)
    
coulomb_stress(10*1e6,0.6,3000) 


  
'''
def plot_coeff(name,file_name,key,method = 'Okada'):
    std1,std2,std3 = cov(file_name,key)
    d1,d2,d3 = average(file_name,key)
    
    d = {'Along-strike':d2,
            'Along-dip':d3,
            'Normal':d1}
    shear = np.sqrt((d['Along-strike'])**2 + (d['Along-dip'])**2)
    normal = d['Normal']
    mu = shear/abs(normal)

    
    patch = models[name]['patch']
    size = models[name]['size']
    shrink = models[name]['shrink']
    title = 'Coefficient of Friction'
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(1,4,patch,nrows,ncols,control=0)
    nrows = len(y)
    ncols = len(x)
    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, ax = plt.subplots(figsize=size,dpi=700)
    im0 = ax.pcolormesh(x,y,mu.reshape(nrows,ncols),edgecolors='k',linewidth=0.25,vmin=0,vmax =10,cmap='rainbow')

    # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
    #X0, Y0 = np.meshgrid(x0, y0)
    #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

    ax.margins(0)
    fig.colorbar(im0, ax=ax,shrink=shrink,label='$\mu$')

    ax.set_ylabel('Down-dip distance (km)',fontsize=10)
    ax.set_xlabel('Along-strike distance (km)',fontsize=10)
    ax.set_aspect('equal', 'box')
    ax.set_title(title,fontweight='bold')
    ax.tick_params(labelsize=10)
    

    fig.suptitle(f'{name}',fontsize=15,fontweight='bold',y = 0.82)
    figure_dir = f'{name}_{title}.png'
    fig.savefig(figure_dir,dpi=700)


        
plot_coeff('Tohoku',Tohoku_file,'Stress change')
plot_coeff('Gorkha',Gorkha_file,'Stress change')
'''