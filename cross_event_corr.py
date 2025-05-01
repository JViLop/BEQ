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
import matplotlib as mpl
import os 




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


def get_mean(x,npatches,geom):
    mean = np.mean(x,axis=1)
    mean_arr = np.flip(mean.reshape(geom[0],geom[1],order='F'),axis=0).flatten()
    return mean_arr

def scale_var(x):
    return (x-np.min(x))/(np.max(x) - np.min(x))




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




lims = {'Tohoku':(0.0,12.5),
              'Iquique':(0.25,4.0),
              'Illapel':(-0.5,5.0),
              'Pedernales':(-5.0,-5.0),
              'Gorkha':(1.25,-0.05)
              }
    
# threshold employed 
th = 30


fig, axes = plt.subplots(2,2, figsize=(14,7.5))

n_lines = 5
cmap = mpl.colormaps['winter']

for j,name in enumerate(list(models.keys())):
    # Take colors at regular intervals spanning the colormap.
    colors = ['black','orange','red','blue','brown'] 
    colors= list(reversed(colors))
    
    
    data_dir = os.path.join(os.getcwd(),f'INPUT/{name}/model/kinematic/{samples}_samples/bin_data/{name}_kinematic_n_{samples}.dat')
    ncols, nrows = models[name]['ncols'],models[name]['nrows']
    
    Np = ncols*nrows
    patch = models[name]['patch']
    nramp = models[name]['nramp']
    mag = models[name]['Mw']
    data = np.fromfile(data_dir,'double').reshape((models[name]['nparams'],samples))
    
    
    xsrc = np.arange((1/2)*patch,ncols*patch,patch)
    ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
    Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
    xsrc_flat = Xsrc.flatten()
    ysrc_flat = Ysrc.flatten()
    
    
    Tr = data[2*Np + nramp:3*Np+nramp,:]
    Vr = data[3*Np+nramp:4*Np+nramp,:] 
    U1 = data[:Np,:]
    U2 = data[Np:2*Np,:]  
    U = np.sqrt(U1**2 + U2**2)
    slip_velocity = np.sqrt((data[:Np,:])**2 + (data[Np:2*Np,:])**2)/data[2*Np:3*Np,:]
    U_mean = np.mean(U,axis=1)
    th_ratio = th/100
    
    ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]


    lns = ['dashed']*n_lines
    
    alphas = [0.3,0.4,0.6,0.8,1.0]
    alphas = list(reversed(alphas))
    

    labels= [('$V_{r}$','$U$'),('$T_{r}$','$U$'),('$V_{r}$','$T_{r}$'),('$V_{r}$','$(U/T_{r})$')]
    units = [('(km/s)','(m)',),('(s)','(m)'),('(km/s)','(s)'),('(km/s)','(m/s)')]

    for k, pair in enumerate([(Vr,U),(Tr,U),(Vr,Tr),(Vr,slip_velocity)]):
    
        
        for i in ids_greater:
            
            data = {'x':pair[0][i,:],'y':pair[1][i,:]}
            data = pd.DataFrame(data)
        
            sns.kdeplot(data,x='x',y='y',ax =axes[k//2][k%2],levels = [1-0.68],linewidths=0.3,color=colors[j],alpha=alphas[j])
        x_mean = np.mean(pair[0],axis=1)
        y_mean = np.mean(pair[1],axis=1)

        if k//2==0 and k%2==0: 
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=4.5,label ='$M_w$'+f'{mag} {name}')
            axes[k//2][k%2].legend(labelspacing=0.1, fontsize=7)
        else:
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[j],alpha = alphas[j],s=4.5)

        axes[k//2][k%2].set_xlabel(f'{labels[k][0]} {units[k][0]}')
        axes[k//2][k%2].set_ylabel(f'{labels[k][1]} {units[k][1]}')
        axes[k//2][k%2].set_box_aspect(1)

plt.subplots_adjust(hspace = 0.25,wspace=-0.55)     
folder = os.path.join(os.getcwd(),f'cross_event_corr')
os.makedirs(folder,exist_ok=True)
file_name = os.path.join(folder,f'cross_event_corr_samples_{samples}_th_{th}.pdf')
fig.savefig(file_name,bbox_inches='tight')
plt.close(fig)






