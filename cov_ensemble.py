# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 14:27:26 2024

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
eq_name = str(sys.argv[1])


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
   for m,comp in enumerate(['N','E','Z']):
       for k in range(nstations):
           file_dir =imodel_dir+'/%.4d.LY%s.sac'%(k,comp)
           stream = obs.read(file_dir)
     
           trace = stream[0].data

           data[i-1,k,:,m] = stream[0].data 
           



N = data[:,:,:,0]  
E = data[:,:,:,1]  

data[:,:,:,0] =  -N*np.sin(strike) + E*np.cos(strike)      
data[:,:,:,1] =  N*np.cos(strike) + E*np.sin(strike)  
# example final time
xtf = data[:,:,-1,0].T
ytf = data[:,:,-1,1].T
ztf= data[:,:,-1,2].T

cov_xtf = np.cov(xtf)        
cov_ytf = np.cov(ytf) 
cov_ztf = np.cov(ztf)      
      
uncertainty = {'x':np.diag(cov_xtf),'y':np.diag(cov_ytf),'z':np.diag(cov_ztf)}


dx = data[0,:,-1,0]
dy = data[0,:,-1,1]
z = data[0,:,-1,2]

dz = z
d = {'x':dx,'y':dy,'z':dz}
fig,axes = plt.subplots(3,1,figsize=(6,10),dpi=400)
for i,key in enumerate(list(d.keys())):
    parameter = d[key]
    #im = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),cmap='bwr')

    im = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(parameter.min(),-parameter.max())),vmax=(max(-parameter.min(),parameter.max()))))
    fig.colorbar(im,ax=axes[i],shrink=0.32)
    axes[i].set_aspect('equal','box')
    axes[i].set_title(f'{key} displacement')
plt.subplots_adjust(hspace=0.2)
fig.suptitle(f'{eq_name} final displacement ($n$={nmodels})',y=0.93,fontweight='bold')
#plt.tight_layout()
plt.savefig(f'{eq_name}_final_displacement_n_{nmodels}.png' )




fig,axes = plt.subplots(3,1,figsize=(6,10),dpi=400)
for i,key in enumerate(list(uncertainty.keys())):
    parameter = uncertainty[key]
    im = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),cmap='rainbow')
    fig.colorbar(im,ax=axes[i],shrink=0.32)
    axes[i].set_aspect('equal','box')
    axes[i].set_title(f'Uncertainty {key} displacement')
plt.subplots_adjust(hspace=0.2)
fig.suptitle(f'{eq_name} Uncertainty final displacement',y=0.93,fontweight='bold')
#plt.tight_layout()
plt.savefig(f'{eq_name}_uncertainty_fk_final_displacement.png' )
    






std_xt = np.zeros((nrows,ncols,npoints))
std_yt = np.zeros((nrows,ncols,npoints))
std_zt = np.zeros((nrows,ncols,npoints))

for t in range(npoints):
    data_xt = data[:,:,t,0]
    data_xt = data_xt.T
    
    data_yt = data[:,:,t,1]
    data_yt = data_yt.T
    
    data_zt= data[:,:,t,2]
    data_zt = data_zt.T

    cov_xt = np.cov(data_xt) 
    cov_yt = np.cov(data_yt) 
    cov_zt = np.cov(data_zt)
    
    std_xt[:,:,t] = np.diag(cov_xt).reshape(nrows,ncols)
    std_yt[:,:,t] = np.diag(cov_yt).reshape(nrows,ncols)
    std_zt[:,:,t] = np.diag(cov_zt).reshape(nrows,ncols)
    
    
    
std = {'x':std_xt,'y':std_yt,'z':std_zt}
import numpy as np
from matplotlib import pyplot as plt, animation
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True


G0 = std_xt[:,:,:]
nframes = G0.shape[2]

    
'''     
for comp in list(std.keys()):         
    fig, ax = plt.subplots(dpi=300)
    parameter = std[comp]
    G0 = parameter[:,:,:]
    
    axtext = fig.add_axes([0.0,0.95,0.1,0.05])
    axtext.axis("off")

    cax = ax.pcolormesh(x, y, G0[:, :, 0],vmin = np.min(G0),vmax =  np.max(G0),cmap='rainbow')
    ax.set_title(f'Uncertainty {comp} displacement')
    ax.set_aspect('equal','box')
    time = axtext.text(0.5,0.25, '$t=$'+str(0)+' s', ha="left", va="top")
    fig.colorbar(cax,shrink=0.4,label='$\sigma(d_{%s}) (m)$'%(comp))
    
    def animate(i):
        cax.set_array(G0[:,:,i].flatten())
        # arr = G0[:,:,i]
        # vmax     = np.max(arr)
        # vmin     = np.min(arr)

        # cf = ax.pcolormesh(x,y,arr, vmax=vmax, vmin=vmin)
        # fig.colorbar(cf, ax=ax)
        time.set_text('$t=$'+str(i)+' s')
    plt.tight_layout()
    anim = animation.FuncAnimation(fig, animate,frames =nframes,interval =50)
    anim.save(filename=f"{comp}_uncertainty_cb.gif", writer="pillow")
    
'''

        
     
fig, axes = plt.subplots(3,1,figsize=(6,10),dpi=300)


axtext = fig.add_axes([0.0,0.95,0.1,0.05])
axtext.axis("off")
comp = 'x'
parameter = std[comp]
G01 = parameter[:,:,:]
cax1 = axes[0].pcolormesh(x, y, G01[:, :, 0],edgecolors='k',linewidth=0.25,vmin = np.min(G01),vmax =  np.max(G01),cmap='rainbow')
axes[0].set_title(f'Uncertainty {comp} displacement')
axes[0].set_aspect('equal','box')
time = axtext.text(0.5,-0.25, '$t=$'+str(0)+' s', ha="left", va="top")
fig.colorbar(cax1,shrink=0.3,label='$\sigma(d_{%s}) (m)$'%(comp))

comp = 'y'
parameter = std[comp]
G02 = parameter[:,:,:]
cax2 = axes[1].pcolormesh(x, y, G02[:, :, 0],edgecolors='k',linewidth=0.25,vmin = np.min(G02),vmax =  np.max(G02),cmap='rainbow')
axes[1].set_title(f'Uncertainty {comp} displacement')
axes[1].set_aspect('equal','box')
fig.colorbar(cax2,shrink=0.3,label='$\sigma(d_{%s}) (m)$'%(comp))

comp = 'z'
parameter = std[comp]
G03 = parameter[:,:,:]
cax3 = axes[2].pcolormesh(x, y, G03[:, :, 0],edgecolors='k',linewidth=0.25,vmin = np.min(G03),vmax =  np.max(G03),cmap='rainbow')
axes[2].set_title(f'Uncertainty {comp} displacement')
axes[2].set_aspect('equal','box')
fig.colorbar(cax3,shrink=0.3,label='$\sigma(d_{%s}) (m)$'%(comp))

def animate(i):
    cax1.set_array(G01[:,:,i].flatten())
    cax2.set_array(G02[:,:,i].flatten())
    cax3.set_array(G03[:,:,i].flatten())
    # arr = G0[:,:,i]
    # vmax     = np.max(arr)
    # vmin     = np.min(arr)

    # cf = ax.pcolormesh(x,y,arr, vmax=vmax, vmin=vmin)
    # fig.colorbar(cf, ax=ax)
    time.set_text('$t=$'+str(i)+' s')
plt.subplots_adjust(hspace=-0.2)
anim = animation.FuncAnimation(fig, animate,frames =nframes,interval =40)
anim.save(filename=f"{eq_name}_uncertainty_{nmodels}.gif", writer="pillow")
    


     
fig, axes = plt.subplots(3,1,figsize=(6,10),dpi=300)


axtext = fig.add_axes([0.0,0.95,0.1,0.05])
axtext.axis("off")
comp = 'x'
parameter = std[comp]
G01 = parameter[:,:,:]
cax1 = axes[0].pcolormesh(x, y, G01[:, :, 0],edgecolors='k',linewidth=0.25,vmin = np.min(G01[:,:,-1]),vmax =  np.max(G01[:,:,-1]),cmap='rainbow')
axes[0].set_title(f'Uncertainty {comp} displacement')
axes[0].set_aspect('equal','box')
time = axtext.text(0.5,-0.25, '$t=$'+str(0)+' s', ha="left", va="top")
fig.colorbar(cax1,shrink=0.3,label='$\sigma(d_{%s}) (m)$'%(comp))

comp = 'y'
parameter = std[comp]
G02 = parameter[:,:,:]
cax2 = axes[1].pcolormesh(x, y, G02[:, :, 0],edgecolors='k',linewidth=0.25,vmin = np.min(G02[:,:,-1]),vmax =  np.max(G02[:,:,-1]),cmap='rainbow')
axes[1].set_title(f'Uncertainty {comp} displacement')
axes[1].set_aspect('equal','box')
fig.colorbar(cax2,shrink=0.3,label='$\sigma(d_{%s}) (m)$'%(comp))

comp = 'z'
parameter = std[comp]
G03 = parameter[:,:,:]
cax3 = axes[2].pcolormesh(x, y, G03[:, :, 0],edgecolors='k',linewidth=0.25,vmin = np.min(G03[:,:,-1]),vmax =  np.max(G03[:,:,-1]),cmap='rainbow')
axes[2].set_title(f'Uncertainty {comp} displacement')
axes[2].set_aspect('equal','box')
fig.colorbar(cax3,shrink=0.3,label='$\sigma(d_{%s}) (m)$'%(comp))

def animate(i):
    cax1.set_array(G01[:,:,i].flatten())
    cax2.set_array(G02[:,:,i].flatten())
    cax3.set_array(G03[:,:,i].flatten())
    # arr = G0[:,:,i]
    # vmax     = np.max(arr)
    # vmin     = np.min(arr)

    # cf = ax.pcolormesh(x,y,arr, vmax=vmax, vmin=vmin)
    # fig.colorbar(cf, ax=ax)
    time.set_text('$t=$'+str(i)+' s')
plt.subplots_adjust(hspace=-0.2)
anim = animation.FuncAnimation(fig, animate,frames =nframes,interval =40)
anim.save(filename=f"{eq_name}_uncertainty_{nmodels}_cb_final_time.gif", writer="pillow")

'''


path = 'statics.neu'
path = '//wsl.localhost/Ubuntu/root/dynamic/MudPy/Tohoku/100_samples/3/output/waveforms/Tohoku/_summary.Tohoku.txt'
nrows,ncols = (9,24)

patch = 29

x = np.arange(patch/2,ncols*patch,patch)
y = np.arange(patch/2,nrows*patch,patch)
y = -np.flip(y)

data = np.loadtxt(path)
strike = 194*np.pi/180
strike = strike - 90.0*np.pi/180.0
n = data[:,3]
e = data[:,4]
z = data[:,5]
dy = n*np.cos(strike) + e*np.sin(strike)
dx = -n*np.sin(strike) +e*np.cos(strike)
dz = z
d = {'x':dx,'y':dy,'z':dz}


fig,axes = plt.subplots(3,1,figsize=(5,8.5),dpi=400)
for i,key in enumerate(list(d.keys())):
    parameter = d[key]
    im = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(parameter.min(),-parameter.max())),vmax=(max(-parameter.min(),parameter.max()))))
    fig.colorbar(im,ax=axes[i],shrink=0.32)
    axes[i].set_aspect('equal','box')
    axes[i].set_title(f'{key} displacement')
plt.subplots_adjust(hspace=-0.4)
fig.suptitle(f'{eq_name} FK final displacement',y=0.85,fontweight='bold')
#plt.tight_layout()
plt.savefig(f'{eq_name}_fk_final_displacement.png' )
    
        

'''















           