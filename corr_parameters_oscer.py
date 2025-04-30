# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 13:08:18 2025

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:00:34 2025

@author: joanv
"""

import matplotlib.pyplot as plt 
import numpy as np
import matplotlib as mpl
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Illapel/model/kinematic/step_052-001.h5'

import h5py

import os 
from matplotlib.colors import TwoSlopeNorm


def model_args(name,ncols,nrows,patchs,nparameters,ramps,nsamples,shrinks,lims):
    config = dict()
    for i,n in enumerate(name):
        config[n] = {}
        config[n]['ncols'] = ncols[i]
        config[n]['nrows'] = nrows[i]
        config[n]['patch'] = patchs[i]
        config[n]['nparams'] = nparameters[i] 
        config[n]['ramps'] = ramps[i]
        config[n]['nsamples'] = nsamples[i]
        config[n]['shrink'] = shrinks[i]
        config[n]['lims'] = lims[i]
            
    return config

    
def get_corr(x,y,npatches,geom):
    z = np.concatenate((x,y))
    corr = np.corrcoef(z)
    diag_corr = np.diagonal(corr[:Np,Np:2*Np])
    diag_corr = np.flip(diag_corr.reshape(geom[0],geom[1],order='F'),axis=0).flatten()
    return diag_corr


def get_mean(x,npatches,geom):
    mean = np.mean(x,axis=1)
    mean_arr = np.flip(mean.reshape(geom[0],geom[1],order='F'),axis=0).flatten()
    return mean_arr

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']

nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patchs = [29,17,18,10,15]
nparams = [866,533,682,650,331]

ramps = [0,3,0,0,9]

 
# names = ['Iquique']

# nrows = [11]
# ncols = [12]
# patchs = [17]
# nparams = [533]


# names = ['Gorkha']

# nrows = [9]
# ncols = [18]
# patchs = [10]
# nparams = [650]


# names = ['Tohoku']

# nrows = [9]
# ncols = [24]
# patchs = [29]
# nparams = [866]

# Pedernales parameters below

# names = ['Pedernales']
# nrows = [8]
# ncols = [10]
# patchs = [15]
# nparams = [331]
# ramps = [9]

nsamples = 100
nsamples = np.ones(len(names),dtype=int)*nsamples
shrinks =  [0.35,0.20,0.25,0.3,0.125]
lims = [(0.0,12.5),(0.25,4.0),(-0.5,5.0),(1,2),(1.25,-0.05)]
models = model_args(names,ncols,nrows,patchs,nparams,ramps,nsamples,shrinks,lims)

plt.ioff()


def scale_var(x):
    return (x-np.min(x))/(np.max(x) - np.min(x))



id_event = int(sys.argv[1])
name = names[id_event]


fig = plt.figure(figsize=(14,7.5),layout="constrained")
gs0 = fig.add_gridspec(1, 2)

gs00 = gs0[0].subgridspec(2, 2,wspace=0.05)
gs01 = gs0[1].subgridspec(4, 1)

ths = [10,20,30]
n_lines = len(ths)
cmap = mpl.colormaps['winter']

# Take colors at regular intervals spanning the colormap.
colors = cmap(np.linspace(0, 1, n_lines))
colors = np.flip(colors,axis=0)

axes1 = []
axes2 = []
for k in range(4):
    
    ax1 = fig.add_subplot(gs00[k//2, k%2])
    
    ax2 = fig.add_subplot(gs01[k])
    axes1.append(ax1)
    axes2.append(ax2)

for m,th in enumerate(ths):
# for name in ['Tohoku']:
    mean_data = f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv'
    df = pd.read_csv(mean_data)
    
    samples = models[name]['nsamples']
    
    data = f'INPUT/{name}/model/kinematic/{samples}_samples/bin_data/{name}_kinematic_n_{samples}.dat'
    ncols, nrows = models[name]['ncols'],models[name]['nrows']
    
    Slip = np.flip(df['Slip'].values.reshape((nrows,ncols),order='F'),axis=0).flatten()
    Np = ncols*nrows
    patch = models[name]['patch']
    data = np.fromfile(data,'double').reshape((models[name]['nparams'],samples))
    
    #colors = np.linspace(0,1,Np)
    #cmap = plt.colormaps.get_cmap('jet')
    
    xsrc = np.arange((1/2)*patch,ncols*patch,patch)
    ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
    Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
    xsrc_flat = Xsrc.flatten()
    ysrc_flat = Ysrc.flatten()
    scaler = MinMaxScaler()
    
    nramp = models[name]['ramps']
    
    Tr = data[2*Np+nramp:3*Np+nramp,:]
    
    Vr = data[3*Np+nramp:4*Np+nramp,:]
    
    U1 = data[:Np,:]
    
  

    U2 = data[Np:2*Np,:]
    
    U = np.sqrt(U1**2 + U2**2)
    slip_velocity = U/Tr
    #slip_velocity[np.isnan(slip_velocity)] = 0.0
    U_mean = np.mean(U,axis=1)
    th_ratio = th/100
    ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]
    Vr_mean = np.mean(Vr,axis=1)
    slip_velocity_mean = np.mean(slip_velocity,axis=1)
    #lws=[0.1,0.4,0.7,1.0]
    lws=[0.5,1,1.7]
    lns = ['dashed']*n_lines
    #colors = ['red','orange','green','blue']
    alphas = [0.4,0.6,0.8,1.0]
    labels= [('$V_{r}$','$U$'),('$T_{r}$','$U$'),('$V_{r}$','$T_{r}$'),('$V_{r}$','$U/T_{r}$')]
    units = [('(km/s)','(m)',),('(s)','(m)'),('(km/s)','(s)'),('(km/s)','(m/s)')]
    print(name)
    for k, pair in enumerate([(Vr,U),(Tr,U),(Vr,Tr),(Vr,slip_velocity)]):
        ax1 = axes1[k]
        ax1.set_box_aspect(1)
        ax2 = axes2[k]

            
        for i in ids_greater:
            
            # axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
            data = {'x':pair[0][i,:],'y':pair[1][i,:]}
            data = pd.DataFrame(data)
        
            sns.kdeplot(data,x='x',y='y',ax =ax1,levels = [1-0.68],linewidths=0.5,color=colors[m])
        #corrcoef = np.corrcoef(pair[0][i,:],pair[1][i,:])
        x_mean = np.mean(pair[0],axis=1)
        y_mean = np.mean(pair[1],axis=1)
        corrcoef = np.corrcoef(x_mean[ids_greater],y_mean[ids_greater])

        #axes[k][0].text(0.97, 0.90, r'$Corr =%s $'%(round(corrcoef[0][1],3)),
        #    verticalalignment='bottom', transform = axes[k][0].transAxes,horizontalalignment='right',
        #    color='green', fontsize=6.5,fontweight='bold')
        
        ax1.scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m],alpha = alphas[m],s=3)
        print(labels[k])
        print('Corr =%s $'%(round(corrcoef[0][1],3)))
        #axes[k][0].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m] ,alpha = alphas[m],s=3,label =r'$Corr =%s $'%(round(corrcoef[0][1],3)))
        ax1.set_xlabel(f'{labels[k][0]} {units[k][0]}')
        ax1.set_ylabel(f'{labels[k][1]} {units[k][1]}')
        #axes[k][0].set_aspect('equal','box')
        if k==0 or k==1:
            
            ax1.axhline(y=th_ratio*np.max(U_mean),linestyle=lns[m],linewidth = lws[m],color='k')
            ax1.text(models[name]['lim'][k],(th_ratio + 0.015)*np.max(U_mean),f'${th_ratio}$'+'$U_{max}$',fontsize=5)
        ax1.legend(loc= 'upper right',fontsize=5)
        #axes[k][0].set_xlim([-0.1,1.1])
        #axes[k][0].set_ylim([-0.1,1.1])
        
        
    
    
    
        corr = get_corr(pair[0],pair[1],Np, (nrows,ncols))
        id_th_col = ids_greater//nrows 
        id_th_row = nrows - ids_greater%nrows -1
        im = ax2.pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
        #im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='viridis')
        #cs =  axes[k//2][k%2+2].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=[0.2,0.5,0.8,1.25],colors='black',linestyles='dashed')
        cs =  ax2.contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=[0.25,0.7,1.25],colors='black',linestyles='dashed')
        props = dict(boxstyle='round', facecolor='white', alpha=0.5)

        str_tuple = labels[k]
        ax2.text(0.03,0.15,str_tuple[1] + '-' + str_tuple[0], transform=ax2.transAxes,fontsize=8,verticalalignment='top', bbox=props)


        #cs =  axes[k][1].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=0.4,colors=colors,alpha=alphas,linestyles='dashed')
      
        #axes[k][1].scatter(xsrc[id_th_col],ysrc[id_th_row],marker='.',color='k',s = 10,label = f'$U > {th_ratio}$'+'$U_{max}$')
        #axes[1].set_title('$U/T_r$',fontsize=9,fontweight='bold')
        if k==0 and m==0:
            ax2.set_xlabel('Along-strike (km)',fontsize=7)
            ax2.set_ylabel('Along-dip (km)',fontsize=7)                         
        if name=='Illapel' and k%2==1 and k//2==1:
            
            ax1.set_ylim([0,2.25])
            
        ax2.set_aspect('equal', 'box')
        if k==3 and m==0:
            fig.colorbar(im,ax=ax2,shrink=models[name]['shrink'],pad=0.3,label= 'Correlation',orientation='horizontal')

                
              
        
folder = os.path.join(os.getcwd(),'global_corr/{samples}')
os.makedirs(folder,exist_ok=True)
file_name = os.path.join(folder,f'{name}_samples_{samples}_test.pdf')
fig.savefig(file_name)
plt.close(fig)



