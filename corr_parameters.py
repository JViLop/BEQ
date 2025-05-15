# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:00:34 2025

@author: joanv
"""

import matplotlib.pyplot as plt 
import numpy as np
import matplotlib as mpl
import pandas as pd
import seaborn as sns
from sklearn.preprocessing import MinMaxScaler
import h5py
import os 
from matplotlib.colors import TwoSlopeNorm


def model_args(name,ncols,nrows,patchs,nparameters,ramps,spacings):
    config = dict()
    for i,n in enumerate(name):
        config[n] = {}
        config[n]['ncols'] = ncols[i]
        config[n]['nrows'] = nrows[i]
        config[n]['patch'] = patchs[i]
        config[n]['nparams'] = nparameters[i] 
        config[n]['ramps'] = ramps[i]
        config[n]['spacings'] = spacings[i]
            
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


def scale_var(x):
    return (x-np.min(x))/(np.max(x) - np.min(x))


names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']
nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patchs = [29,17,18,10,15]
nparams = [866,533,682,650,331]
ramps = [0,3,0,0,9]
xlimits = [(0,8),(0,3.5),(0,5),(0,4.5),(0,5)]
ylimits = [(0,3.5),(0,1),(0,4),(0,5.5),(0,3.5)]
spacings = [0,-0.25,-0.125,-0.05,-0.2]

models = model_args(names,ncols,nrows,patchs,nparams,ramps,spacings)
samples = 100





lims = {'Tohoku':(0.0,12.5),
              'Iquique':(0.25,4.0),
              'Illapel':(-0.5,5.0),
              'Pedernales':(1,2),
              'Gorkha':(1.25,-0.05)}

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']


for j,name in enumerate(names):
    #fig = plt.figure(figsize=(14,8),layout="constrained")
    fig = plt.figure(figsize=(14,9))
    spacing = models[name]['spacings']
    gs0 = fig.add_gridspec(1, 2,wspace=spacing)

    gs00 = gs0[0].subgridspec(2, 2,wspace=0.27,hspace= 0.2)
    gs01 = gs0[1].subgridspec(4, 1,hspace=0.3)

    # thresholds
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

        
        
        data_dir = os.path.join(os.getcwd(),f'INPUT/{name}/model/kinematic/{samples}_samples/bin_data/{name}_kinematic_n_{samples}.dat')
        ncols, nrows = models[name]['ncols'],models[name]['nrows']
        
        Np = ncols*nrows
        patch = models[name]['patch']
        data = np.fromfile(data_dir,'double').reshape((models[name]['nparams'],samples))
        
  
        xsrc = np.arange((1/2)*patch,ncols*patch,patch)
        ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
        Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
        xsrc_flat = Xsrc.flatten()
        ysrc_flat = Ysrc.flatten()
        scaler = MinMaxScaler()
        lws=[0.5,1,1.7]
        lns = ['dashed']*n_lines
        
        # transparency of each 
        alphas = [0.4,0.6,0.8,1.0]
        labels= [('$V_{r}$','$U$'),('$T_{r}$','$U$'),('$V_{r}$','$T_{r}$'),('$V_{r}$','$U/T_{r}$')]
        units = [('(km/s)','(m)',),('(s)','(m)'),('(km/s)','(s)'),('(km/s)','(m/s)')]
        nramp = models[name]['ramps']        

        
        # actual data 
        
        Tr = data[2*Np+nramp:3*Np+nramp,:]
        Vr = data[3*Np+nramp:4*Np+nramp,:]
        U1 = data[:Np,:]
        U2 = data[Np:2*Np,:]
        U = np.sqrt(U1**2 + U2**2) 
        slip_velocity = U/Tr

        U_mean = np.mean(U,axis=1)
        th_ratio = th/100
        
        # retrieve patch ids satisfying threshold condition
        ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]

        Slip = np.flip(U_mean.reshape((nrows,ncols),order='F'),axis=0).flatten()
    

        for k, pair in enumerate([(Vr,U),(Tr,U),(Vr,Tr),(Vr,slip_velocity)]):
            ax1 = axes1[k]
            ax1.set_box_aspect(1)
            ax2 = axes2[k]

                
            for i in ids_greater:
                
                data = {'x':pair[0][i,:],'y':pair[1][i,:]}
                data = pd.DataFrame(data)
            
                sns.kdeplot(data,x='x',y='y',ax =ax1,levels = [1-0.68],linewidths=0.5,color=colors[m])
            
            x_mean = np.mean(pair[0],axis=1)
            y_mean = np.mean(pair[1],axis=1)
            corrcoef = np.corrcoef(x_mean[ids_greater],y_mean[ids_greater])


            ax1.scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m],alpha = alphas[m],s=3,label =r'$\rho =%s $'%(round(corrcoef[0][1],3)))

            ax1.set_xlabel(f'{labels[k][0]} {units[k][0]}')
            ax1.set_ylabel(f'{labels[k][1]} {units[k][1]}')
            if k==0 or k==1:
                
                ax1.axhline(y=th_ratio*np.max(U_mean),linestyle=lns[m],linewidth = lws[m],color='k')
                ax1.text(lims[name][k],(th_ratio + 0.015)*np.max(U_mean),f'${th_ratio}$'+'$U_{max}$',fontsize=5)
            ax1.legend(loc= 'upper right',fontsize=5)
        

        
            corr = get_corr(pair[0],pair[1],Np, (nrows,ncols))
            # id_th_col = ids_greater//nrows 
            # id_th_row = nrows - ids_greater%nrows -1
            im = ax2.pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))

            cs =  ax2.contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=[0.25,0.7,1.25],colors='black',linestyles='dashed')
            props = dict(boxstyle='round', facecolor='white', alpha=0.5)

            str_tuple = labels[k]
            ax2.text(0.025,0.15,str_tuple[1] + '-' + str_tuple[0], transform=ax2.transAxes,fontsize=8,verticalalignment='top', bbox=props)


            if k==0 and m==0:
                ax2.set_xlabel('Along-strike (km)',fontsize=7)
                ax2.set_ylabel('Along-dip (km)',fontsize=7)                         
            elif k!=0:
                ax2.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
                
            if name=='Illapel' and k%2==1 and k//2==1:
                
                ax1.set_ylim([0,2.25])
                
            ax2.set_aspect('equal', 'box')
            if k==3 and m==0:
                box = ax2.get_position()
                cb_ax = fig.add_axes([box.x0 ,0.85*box.y0,box.width,box.height/13])
                fig.colorbar(im,ax = ax2, orientation='horizontal',cax=cb_ax,ticks = [-1,-0.5,0, 0.5, 1],label='Correlation')    

                    
                  
    print(name)

            
       
    folder = os.path.join(os.getcwd(),'individual_corr')
    os.makedirs(folder,exist_ok=True)
    file_name = os.path.join(folder,f'{name}_samples_{samples}.pdf')
    fig.savefig(file_name,bbox_inches='tight')
    plt.close(fig)



