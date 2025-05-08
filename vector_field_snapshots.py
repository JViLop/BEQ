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

names =  ['Iquique','Pedernales','Gorkha']
nx = len(names)
width = 4 
fig, axes = plt.subplots(nx,nepochs,figsize=(int(nepochs*width),int(width*nx)),dpi=700)
nsamples = 100
nmodels = 1
npoints = 256

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
    print(DATA.shape)
    for k,epoch in enumerate(list(np.logspace(3,np.log2(256),num=nepochs,base=2))):
             
        if epoch==256:
            epoch-=1
        
        Ux = DATA[:,int(epoch),0] 
        Uy = DATA[:,int(epoch),1] 
        Uz = DATA[:,int(epoch),2] 
        U = np.sqrt(Ux**2 + Uy**2)
        print(max(Ux),max(Uy))

        if nx ==1:
            ax=axes[k]

        else:
            ax=axes[n][k]

        cax1 = ax.pcolormesh(x, y, Uz.reshape(nrows,ncols),linewidth=0.25,vmin = min(min(Uz),-max(Uz)),vmax = max(-min(Uz),max(Uz)),cmap='bwr')

        q = ax.quiver(X,Y,Ux,Uy,scale=scale,scale_units ='x', units='width',width=0.0035,headwidth=2.5,headlength=4.0)
        ax.set_clip_on(False)
        
        ax.quiverkey(q, 0.87,1.09 , U=np.ceil(U.max()),
                      label=' {} m'.format(int(np.ceil(U.max()))), labelpos='E',angle=0)
        factor = 1/q.scale    # scaling factor from UV unit to XY units
        
        offsetXY = np.column_stack((Ux,Uy))
        
        # minus sign in Y coordinate corrects for y-axis inversion
        
        offsetXY = np.array([[xy[0],xy[1]] for xy in offsetXY])
        #coord = q.XY + offsetXY*factor
        # width and height had to be formatted in row-like order (initially in column-like order)
#        ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
#                        width=dictionary['std_U_perp'][i]*factor,
#                        height=dictionary['std_U_parallel'][i]*factor,
#                        angle=0,alpha=0.5,fill=False,edgecolor='k')
#                for i in range(q.N)]
#        for e in ells:
#            ax.add_artist(e)

  


        ax.text(0.05, 1.25, '$t = %s s$'%(int(round(epoch,0))), transform=ax.transAxes,fontsize = 9.5,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4'))
        
        ax.set_aspect('equal','box')
        if k==0:
            ax.set_title(f'{eq_name}', fontweight='bold',fontsize=14)
            ax.set_ylabel('Trench-normal distance (km)',fontsize=7)
            ax.set_xlabel('Along-strike distance (km)',fontsize=7)
            
        fig.colorbar(cax1,shrink=shrink,label='$d_{z} (m)$',pad=0.03,orientation='horizontal')
        fig.tight_layout()
fig.savefig('snapshot_xyz.png')


