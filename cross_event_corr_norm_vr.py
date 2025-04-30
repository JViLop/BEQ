# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 18:30:06 2025

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:16:39 2025

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 15:20:16 2024

@author: joanv
"""


import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from matplotlib.colors import TwoSlopeNorm
import re


def model_args(name,ncols,nrows,patchs,nparameters,nramps,magnitudes):
    config = dict()
    for i,n in enumerate(name):
        config[n] = {}
        config[n]['ncols'] = ncols[i]
        config[n]['nrows'] = nrows[i]
        config[n]['patch'] = patchs[i]
        config[n]['nparams'] = nparameters[i] 
        config[n]['nramp'] = nramps[i]
        config[n]['Mw'] = magnitudes[i]
        
            
    return config

    
def get_corr(x,y,npatches,geom):
    z = np.concatenate((x,y))
    corr = np.corrcoef(z)
    diag_corr = np.diagonal(corr[:Np,Np:2*Np])
    diag_corr = np.flip(diag_corr.reshape(geom[0],geom[1],order='F'),axis=0).flatten()
    return diag_corr


def get_vs(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VP'][i]    # in km/s
        return vs

def get_mean(x,npatches,geom):
    mean = np.mean(x,axis=1)
    mean_arr = np.flip(mean.reshape(geom[0],geom[1],order='F'),axis=0).flatten()
    return mean_arr

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']
magnitudes = [9.0,8.1,8.3,7.8,7.8]
nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patchs = [29,17,18,10,15]
nparams = [866,533,682,650,331]


nramps = [0,3,0,0,9]


xlimits = [(0,8),(0,3.5),(0,5),(0,4.5),(0,5)]
ylimits = [(0,3.5),(0,1),(0,4),(0,5.5),(0,3.5)]

models = model_args(names,ncols,nrows,patchs,nparams,nramps,magnitudes)
samples = 100
plt.ioff()


def scale_var(x):
    return (x-np.min(x))/(np.max(x) - np.min(x))

import matplotlib as mpl


lims = {'Tohoku':(0.0,12.5),
              'Iquique':(0.25,4.0),
              'Illapel':(-0.5,5.0),
              'Pedernales':(-5.0,-5.0),
              'Gorkha':(1.25,-0.05)
              }
#or th in [0,5,10,15,20,25,30,35,40]:

#for th in [30]:
th = 20


'''
fig, axes = plt.subplots(2,2, figsize=(14,7.5))

n_lines = 5
cmap = mpl.colormaps['winter']



for j,name in enumerate(list(models.keys())):
    # Take colors at regular intervals spanning the colormap.
    colors = ['black','orange','red','blue','brown'] 
    colors= list(reversed(colors))
    
    
    earth_model_dir = os.path.join(f'1d_models/model_{name}')
    f = open(earth_model_dir, "r")
    content = f.readlines()
    df = pd.read_csv(f'INPUT/{name}/{samples}_samples/{name}_mean_kinematic_model.csv')

    #nlayers =int(sys.argv[4])

    ncols, nrows = models[name]['ncols'],models[name]['nrows']

    shape = (nrows,ncols)

    earth_model = content[12:]
    keys = re.findall('\w+',content[11],flags=re.DOTALL)
    fmt = re.compile('\d+.\d+')
    layers = [list(map(float,re.findall(fmt,m))) for m in earth_model]
    #print(layers)
    layers = np.array(layers).T
    layers0 = np.zeros(layers.shape)
    layers0[0,:] = layers[0,:]
    layers0[1,:] = layers[2,:]
    layers0[2,:] = layers[1,:]
    layers0[3,:] = layers[3,:]
    layers0[4,:] = layers[5,:]
    layers0[5,:] = layers[4,:]
    model_dict = dict(zip(keys,layers0))
    #print(layers0)
    
    model_df = pd.DataFrame(model_dict)
    
    print(model_df)
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0).flatten() 
    VS =  np.ones_like(DEPTH)

    for q in range(len(VS)):
        vs = get_vs(dict(model_df),DEPTH[q])
        VS[q] = vs
    
# for name in ['Tohoku']:
    mean_data = f'INPUT/{name}/{samples}_samples/{name}_mean_kinematic_model.csv'
    df = pd.read_csv(mean_data)
    
    
    data = f'INPUT/{name}/{samples}_samples/{name}_kinematic_n_{samples}.dat'
    ncols, nrows = models[name]['ncols'],models[name]['nrows']
    
    Slip = np.flip(df['Slip'].values.reshape((nrows,ncols),order='F'),axis=0).flatten()
    Np = ncols*nrows
    patch = models[name]['patch']
    nramp = models[name]['nramp']
    mag = models[name]['Mw']
    data = np.fromfile(data,'double').reshape((models[name]['nparams'],samples))
    
    #colors = np.linspace(0,1,Np)
    #cmap = plt.colormaps.get_cmap('jet')
    
    xsrc = np.arange((1/2)*patch,ncols*patch,patch)
    ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
    Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
    xsrc_flat = Xsrc.flatten()
    ysrc_flat = Ysrc.flatten()
    
    
    Tr = data[2*Np + nramp:3*Np+nramp,:]
    
    Vr = data[3*Np+nramp:4*Np+nramp,:]
    Vr_norm = Vr.T/VS
    Vr_norm = Vr_norm.T
    
    U1 = data[:Np,:]
    
    
    U2 = data[Np:2*Np,:]
    
    U = np.sqrt(U1**2 + U2**2)
    slip_velocity = np.sqrt((data[:Np,:])**2 + (data[Np:2*Np,:])**2)/data[2*Np:3*Np,:]
    #slip_velocity[np.isnan(slip_velocity)] = 0.0
    U_mean = np.mean(U,axis=1)
    th_ratio = th/100
    ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]
    Vr_mean = np.mean(Vr,axis=1)
    Vr_norm_mean = np.mean(Vr_norm,axis=1) 
    slip_velocity_mean = np.mean(slip_velocity,axis=1)
    #lws=[0.1,0.4,0.7,1.0]
    lns = ['dashed']*n_lines
    #colors = ['red','orange','green','blue']
    alphas = [0.3,0.4,0.6,0.8,1.0]
    alphas = list(reversed(alphas))
    #alphas = [1.0,0.8,0.6,0.4]

    units = [('(km/s)','(m)')]
    labels= [('$V_{r}/V_{s}$','$U$'),('$T_{r}$','$U$'),('$V_{r}/V_{s}$','$T_{r}$'),('$V_{r}/V_{s}$','$U/T_{r}$')]
    units = [('','(m)',),('(s)','(m)'),('','(s)'),('','(m/s)')]

    for k, pair in enumerate([(Vr_norm,U),(Tr,U),(Vr_norm,Tr),(Vr_norm,slip_velocity)]):
    
        
        for i in ids_greater:
            
            # axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
            data = {'x':pair[0][i,:],'y':pair[1][i,:]}
            data = pd.DataFrame(data)
        
            sns.kdeplot(data,x='x',y='y',ax =axes[k//2][k%2],levels = [1-0.68],linewidths=0.3,color=colors[j],alpha=alphas[j])
        x_mean = np.mean(pair[0],axis=1)
        y_mean = np.mean(pair[1],axis=1)
        #corrcoef = np.corrcoef(x_mean[ids_greater],y_mean[ids_greater])
        #ax.text(0.0, 0.90, r'$th =%s $'%(round(corrcoef[0][1],3)),verticalalignment='bottom', transform = ax.transAxes,horizontalalignment='right',color='green', fontsize=6.5,fontweight='bold')
        if k//2==0 and k%2==0: 
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=4.5,label ='$M_w$'+f'{mag} {name}')
            axes[k//2][k%2].legend(labelspacing=0.1, fontsize=7)
        else:
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=4.5)

        #ax.scatter(np.mean(x_mean[ids_greater]),np.mean(y_mean[ids_greater]),s=8,marker = 's',color=colors[j],label=f'{name}')
        #ax.scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=3,label =r'$Corr_{%s} =%s $'%(name,round(corrcoef[0][1],3)))

        #axes[k][0].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m] ,alpha = alphas[m],s=3,label =r'$Corr =%s $'%(round(corrcoef[0][1],3)))
        #axes[k//2][k%2].set_title(f'{labels[k][1]} vs {labels[k][0]} for ' + f'$U >{th_ratio}$'+'$U_{max}$',fontsize=8.75,fontweight='bold')
        axes[k//2][k%2].set_xlabel(f'{labels[k][0]} {units[k][0]}')
        axes[k//2][k%2].set_ylabel(f'{labels[k][1]} {units[k][1]}')
        axes[k//2][k%2].set_box_aspect(1)
        axes[k//2][k%2].text(0.05,1.04,'abcd'[k],transform=axes[k//2][k%2].transAxes,fontsize=18,fontweight='bold') 
        plt.subplots_adjust(hspace = 0.25,wspace=-0.5)

        
folder = os.path.join(os.getcwd(),f'Cross_event/{samples}')
os.makedirs(folder,exist_ok=True)
file_name = os.path.join(folder,f'cross_event_corr_samples_{samples}_th_{th}_vrnorm.pdf')
fig.savefig(file_name)
plt.close(fig)
'''

fig, axes = plt.subplots(2,4, figsize=(10,6),layout = 'constrained')

n_lines = 5
cmap = mpl.colormaps['winter']

name = 'Tohoku'

for j,name in enumerate([name]):
    # Take colors at regular intervals spanning the colormap.
    colors = ['black','orange','red','blue','brown'] 
    colors= list(reversed(colors))
    
    
    earth_model_dir = os.path.join(f'1d_models/model_{name}')
    f = open(earth_model_dir, "r")
    content = f.readlines()
    df = pd.read_csv(f'INPUT/{name}/{samples}_samples/{name}_mean_kinematic_model.csv')

    #nlayers =int(sys.argv[4])

    ncols, nrows = models[name]['ncols'],models[name]['nrows']

    shape = (nrows,ncols)

    earth_model = content[12:]
    keys = re.findall('\w+',content[11],flags=re.DOTALL)
    fmt = re.compile('\d+.\d+')
    layers = [list(map(float,re.findall(fmt,m))) for m in earth_model]
    #print(layers)
    layers = np.array(layers).T
    layers0 = np.zeros(layers.shape)
    layers0[0,:] = layers[0,:]
    layers0[1,:] = layers[2,:]
    layers0[2,:] = layers[1,:]
    layers0[3,:] = layers[3,:]
    layers0[4,:] = layers[5,:]
    layers0[5,:] = layers[4,:]
    model_dict = dict(zip(keys,layers0))
    #print(layers0)
    
    model_df = pd.DataFrame(model_dict)
    
    print(model_df)
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0).flatten() 
    VS =  np.ones_like(DEPTH)

    for q in range(len(VS)):
        vs = get_vs(dict(model_df),DEPTH[q])
        VS[q] = vs
    
# for name in ['Tohoku']:
    mean_data = f'INPUT/{name}/{samples}_samples/{name}_mean_kinematic_model.csv'
    df = pd.read_csv(mean_data)
    
    
    data = f'INPUT/{name}/{samples}_samples/{name}_kinematic_n_{samples}.dat'
    ncols, nrows = models[name]['ncols'],models[name]['nrows']
    
    Slip = np.flip(df['Slip'].values.reshape((nrows,ncols),order='F'),axis=0).flatten()
    Np = ncols*nrows
    patch = models[name]['patch']
    nramp = models[name]['nramp']
    mag = models[name]['Mw']
    data = np.fromfile(data,'double').reshape((models[name]['nparams'],samples))
    
    #colors = np.linspace(0,1,Np)
    #cmap = plt.colormaps.get_cmap('jet')
    
    xsrc = np.arange((1/2)*patch,ncols*patch,patch)
    ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
    Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
    xsrc_flat = Xsrc.flatten()
    ysrc_flat = Ysrc.flatten()
    
    
    Tr = data[2*Np + nramp:3*Np+nramp,:]
    
    Vr = data[3*Np+nramp:4*Np+nramp,:]
    Vr_norm = Vr.T/VS
    Vr_norm = Vr_norm.T
    
    U1 = data[:Np,:]
    
    
    U2 = data[Np:2*Np,:]
    
    U = np.sqrt(U1**2 + U2**2)
    slip_velocity = np.sqrt((data[:Np,:])**2 + (data[Np:2*Np,:])**2)/data[2*Np+nramp:3*Np+nramp,:]
    #slip_velocity[np.isnan(slip_velocity)] = 0.0
    U_mean = np.mean(U,axis=1)
    th_ratio = th/100
    ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]
    
    Vr_mean = np.mean(Vr,axis=1)
    Vr_norm_mean = np.mean(Vr_norm,axis=1) 
    Vr_norm_mean_ordered = np.flip(Vr_norm_mean.reshape((nrows,ncols),order='F'),axis=0).flatten()
    slip_velocity_mean = np.mean(slip_velocity,axis=1)
    ids_greater_umean = np.where(Slip>th_ratio*np.max(Slip))[0]
    ids_greater_vnorm = np.where(Vr_norm_mean_ordered>1)[0]
    ids_needed = np.intersect1d(ids_greater_umean,ids_greater_vnorm)
    #lws=[0.1,0.4,0.7,1.0]
    lns = ['dashed']*n_lines
    #colors = ['red','orange','green','blue']
    alphas = [0.3,0.4,0.6,0.8,1.0]
    alphas = list(reversed(alphas))
    #alphas = [1.0,0.8,0.6,0.4]

    units = [('(km/s)','(m)')]
    labels= [('$V_{r}/V_{s}$','$U$'),('$T_{r}$','$U$'),('$V_{r}/V_{s}$','$T_{r}$'),('$V_{r}/V_{s}$','$U/T_{r}$')]
    units = [('','(m)',),('(s)','(m)'),('','(s)'),('','(m/s)')]

    for k, pair in enumerate([(Vr_norm,U),(Tr,U),(Vr_norm,Tr),(Vr_norm,slip_velocity)]):
    
        
        for i in ids_greater:
            
            # axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
            data = {'x':pair[0][i,:],'y':pair[1][i,:]}
            data = pd.DataFrame(data)
        
            sns.kdeplot(data,x='x',y='y',ax =axes[k//2][k%2],levels = [1-0.68],linewidths=0.3,color=colors[j],alpha=alphas[j])
        x_mean = np.mean(pair[0],axis=1)
        y_mean = np.mean(pair[1],axis=1)
        #corrcoef = np.corrcoef(x_mean[ids_greater],y_mean[ids_greater])
        #ax.text(0.0, 0.90, r'$th =%s $'%(round(corrcoef[0][1],3)),verticalalignment='bottom', transform = ax.transAxes,horizontalalignment='right',color='green', fontsize=6.5,fontweight='bold')
        if k//2==0 and k%2==0: 
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=4.5,label ='$M_w$'+f'{mag} {name}')
            axes[k//2][k%2].legend(labelspacing=0.1, fontsize=7)
        else:
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=4.5)

        #ax.scatter(np.mean(x_mean[ids_greater]),np.mean(y_mean[ids_greater]),s=8,marker = 's',color=colors[j],label=f'{name}')
        #ax.scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=3,label =r'$Corr_{%s} =%s $'%(name,round(corrcoef[0][1],3)))

        #axes[k][0].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m] ,alpha = alphas[m],s=3,label =r'$Corr =%s $'%(round(corrcoef[0][1],3)))
        #axes[k//2][k%2].set_title(f'{labels[k][1]} vs {labels[k][0]} for ' + f'$U >{th_ratio}$'+'$U_{max}$',fontsize=8.75,fontweight='bold')
        axes[k//2][k%2].set_xlabel(f'{labels[k][0]} {units[k][0]}')
        axes[k//2][k%2].set_ylabel(f'{labels[k][1]} {units[k][1]}')
        axes[k//2][k%2].set_box_aspect(1)
        
        corr = get_corr(pair[0],pair[1],Np, (nrows,ncols))
        axes[k//2][k%2].text(0.05,1.04,'abcd'[k],transform=axes[k//2][k%2].transAxes,fontsize=18,fontweight='bold') 
        m = axes[k//2][k%2+2].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
        axes[k//2][k%2+2].scatter(xsrc_flat[ids_greater_umean],ysrc_flat[ids_greater_umean],marker='.',color='k',s = 12)
        #im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='viridis')
        #cs =  axes[k//2][k%2+2].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=[0.2,0.5,0.8,1.25],colors='black',linestyles='dashed')
        cs =  axes[k//2][k%2+2].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels= 4,colors='black',linestyles='dashed',linewidths=0.02)

        #cs =  axes[k][1].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=0.4,colors=colors,alpha=alphas,linestyles='dashed')

        #axes[k][1].scatter(xsrc[id_th_col],ysrc[id_th_row],marker='.',color='k',s = 10,label = f'$U > {th_ratio}$'+'$U_{max}$')
        #axes[1].set_title('$U/T_r$',fontsize=9,fontweight='bold')

        axes[k//2][k%2+2].set_xlabel('Along-strike (km)',fontsize=8)
        axes[k//2][k%2+2].set_ylabel('Along-dip (km)',fontsize=8) 
        axes[k//2][k%2+2].set_aspect('equal', 'box')

folder = os.path.join(os.getcwd(),f'Cross_event/{samples}')
os.makedirs(folder,exist_ok=True)
file_name = os.path.join(folder,f'cross_event_corr_samples_{samples}_th_{th}_vrnorm_spatial_{name}.pdf')
fig.savefig(file_name)
plt.close(fig)
        
# folder = os.path.join(os.getcwd(),f'Cross_event/{samples}')
# os.makedirs(folder,exist_ok=True)
# file_name = os.path.join(folder,f'cross_event_corr_samples_{samples}_th_{th}_vrnorm.pdf')
# fig.savefig(file_name)
# plt.close(fig)
    

