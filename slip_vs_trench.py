# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:31:14 2024

@author: joanv
"""

# Changes of slip 

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import h5py
import matplotlib.pyplot as plt
import os

import matplotlib as mpl

names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [0,0.5,0.5,0.5,0.4]
wspaces = [0.30,0.5,0.5,0.5,0.3]
sizes = [(14,6),(10,8),(10,6.0),(10,8),(14,7)]
shrinks = [0.5,0.7,0.7,0.75,0.75]
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

def array_formatter(array,geom):
        return np.flip(array.reshape(geom,order='F'),axis=0).flatten()
samples = 100
file_folder ='INPUT'
fig,axes = plt.subplots(5,1,figsize=(10,10))
for i,name in enumerate(names):
    patch = models[name]['patch']
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1] 
    event_folder = os.path.join(f'INPUT/{name}/model/kinematic',f'{samples}_samples/mean')
    file_dir = os.path.join(event_folder,f'{name}_mean_kinematic_model.csv')
    df = pd.read_csv(file_dir)
    df = df.drop(df.columns[0],axis=1)
    
    Slip = array_formatter(df['Slip'].values,(nrows,ncols))
    std_Slip = array_formatter(df['std_U'].values,(nrows,ncols))
    xsrc = np.arange(patch/2 , ncols*patch,patch)
    ysrc = np.arange(-(nrows-1/2)*patch,0,patch)
    Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)

    n0cols,n0rows = len(xsrc),len(ysrc)
    
    max_slip_patch = np.argmin(abs(Slip-max(Slip)))

    
    col0_max_slip,row0_max_slip = max_slip_patch%n0cols, max_slip_patch//n0cols
    offset = 75
    
    # displacement is defined on customized geometry overlapping prescribed surface 
    ncol_target = col0_max_slip
    nrow_target = row0_max_slip 
    max_slip_patch_in_stn = nrow_target*ncols + ncol_target 
    if name=='Gorkha':
        ysrc = ysrc  - offset

    m = {'m':Slip.reshape(nrows,ncols)[:,ncol_target] ,'err_m':std_Slip.reshape(nrows,ncols)[:,ncol_target]}
    
    y0  = ysrc
    dysrcleft = np.abs(y0[0] - y0[1])/2
    dysrcright = np.abs(y0[-1] -y0[-2])/2 
    ysrcleft = np.array([y0[0]-dysrcleft,y0[0]])
    ysrcright = np.array([y0[-1],y0[-1] + dysrcright])
    mleft = np.array([m['m'][0],m['m'][0]])
    mright = np.array([m['m'][-1],m['m'][-1]])
    

        # Create the step function
        
    axes[i].plot(ysrcleft, mleft,linewidth=0.8,color='black')
    axes[i].plot(ysrcright, mright,linewidth=0.8,color='black')
    axes[i].fill_between(ysrcleft, mleft-m['err_m'][0],mleft+m['err_m'][0],step='pre',facecolor='grey',alpha=0.75)
    axes[i].fill_between(ysrcright, mright-m['err_m'][-1],mright+m['err_m'][-1],step='post',facecolor='grey',alpha=0.75)
    axes[i].plot(ysrc,m['m'],drawstyle='steps-mid',marker = '.',linewidth=1.25,ms=2.5, color='black')
    axes[i].fill_between(ysrc,m['m'] + m['err_m'] ,m['m'] - m['err_m'],step='mid',facecolor='grey',alpha=0.75)
    axes[i].text(0.05,0.8,f'{name}',transform=axes[i].transAxes,fontsize=12.5,fontweight='bold')
    if name !='Gorkha':
        axes[i].axvline(x = 0,ls='--',linewidth=0.6)
    axes[i].axhline(y = 0,ls='-',linewidth=0.4,color='black')
    axes[i].set_ylabel('Slip (m) ',fontsize=9)
    if name =='Gorkha':

        axes[i].text(-offset, 4, "Trench",
        ha="center", va="center", rotation=0, size=8,
        bbox=dict(boxstyle="rarrow,pad=0.3",
                  fc="white", ec="black", lw=0.6))
        # ax3.annotate((0,0.5),(2,0.5),ha="center", va="center", rotation=0, size=4,
        # bbox=dict(boxstyle="rarrow,pad=0.3",fc="blue", ec="black", lw=2))
    else:
        axes[i].annotate('Trench',(0,0.5),(2,0.5),fontsize=9,rotation = 90)
    axes[i].set_xlabel('Distance from trench (km)',fontsize=9)
    axes[i].tick_params(labelsize=10)
    fig.align_ylabels(fig.get_axes())
    fig.tight_layout()
    
fig.savefig('change_slip_vs_trench.png',dpi=600)



fig,axes = plt.subplots(5,1,figsize=(10,10))
for i,name in enumerate(names):
    patch = models[name]['patch']
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1] 
    event_folder = os.path.join(f'INPUT/{name}/model/kinematic',f'{samples}_samples/mean')
    file_dir = os.path.join(event_folder,f'{name}_mean_kinematic_model.csv')
    df = pd.read_csv(file_dir)
    df = df.drop(df.columns[0],axis=1)
    
    Slip = array_formatter(df['Slip'].values,(nrows,ncols))
    std_Slip = array_formatter(df['std_U'].values,(nrows,ncols))
    xsrc = np.arange(patch/2 , ncols*patch,patch)
    ysrc = np.arange(-(nrows-1/2)*patch,0,patch)
    Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)

    n0cols,n0rows = len(xsrc),len(ysrc)
    
    max_slip_patch = np.argmin(abs(Slip-max(Slip)))

    
    col0_max_slip,row0_max_slip = max_slip_patch%n0cols, max_slip_patch//n0cols
    offset = 75
    
    # displacement is defined on customized geometry overlapping prescribed surface 
    ncol_target = col0_max_slip
    nrow_target = row0_max_slip 
    max_slip_patch_in_stn = nrow_target*ncols + ncol_target 
    if name=='Gorkha':
        ysrc = ysrc  - offset

    m = {'m':Slip.reshape(nrows,ncols)[:,ncol_target] ,'err_m':std_Slip.reshape(nrows,ncols)[:,ncol_target]}
    rel_m = m['err_m']/np.max(m['m'])
    y0  = ysrc
    dysrcleft = np.abs(y0[0] - y0[1])/2
    dysrcright = np.abs(y0[-1] -y0[-2])/2 
    ysrcleft = np.array([y0[0]-dysrcleft,y0[0]])
    ysrcright = np.array([y0[-1],y0[-1] + dysrcright])
    rel_mleft = np.array([rel_m[0],rel_m[0]])
    rel_mright = np.array([rel_m[-1],rel_m[-1]])
    

        # Create the step function
        
    axes[i].plot(ysrcleft, rel_mleft,linewidth=1.25,color='black')
    axes[i].plot(ysrcright, rel_mright,linewidth=1.25,color='black')
    axes[i].plot(ysrc,rel_m,drawstyle='steps-mid',marker = '.',linewidth=1.25,ms=2.5, color='black')
    #axes[i].fill_between(ysrc,m['m'] + m['err_m'] ,m['m'] - m['err_m'],step='mid',facecolor='grey',alpha=0.75)
    axes[i].text(0.05,0.8,f'{name}',transform=axes[i].transAxes,fontsize=12.5,fontweight='bold')
    if name !='Gorkha':
        axes[i].axvline(x = 0,ls='--',linewidth=0.6)
    axes[i].axhline(y = 0,ls='-',linewidth=0.4,color='black')
    axes[i].axvline(x = ysrc[nrow_target],ls='-',linewidth=0.75,color='blue')
    axes[i].annotate('$U_{max}$',(ysrc[nrow_target],0.75*np.mean(rel_m)),(ysrc[nrow_target],0.75*np.mean(rel_m)),fontsize=9,rotation = 0)
    axes[i].set_ylabel('$\sigma_{U}/U_{max}$',fontsize=12)
    if name =='Gorkha':

        axes[i].text(-offset, 0.75*np.mean(rel_m), "Trench",
        ha="center", va="center", rotation=0, size=8,
        bbox=dict(boxstyle="rarrow,pad=0.3",
                  fc="white", ec="black", lw=0.6))
        # ax3.annotate((0,0.5),(2,0.5),ha="center", va="center", rotation=0, size=4,
        # bbox=dict(boxstyle="rarrow,pad=0.3",fc="blue", ec="black", lw=2))
    else:
        axes[i].annotate('Trench',(0,0.75*np.mean(rel_m)),(2,0.75*np.mean(rel_m)),fontsize=9,rotation = 90)
    axes[i].set_xlabel('Distance from trench (km)',fontsize=12)
    axes[i].tick_params(labelsize=12)
    fig.align_ylabels(fig.get_axes())
    fig.tight_layout()
    
fig.savefig('change_rel_slip_vs_trench.png',dpi=600)