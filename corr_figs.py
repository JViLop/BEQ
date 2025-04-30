# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 15:20:16 2024

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 11:14:25 2024

@author: joanv
"""

import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler

data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Illapel/model/kinematic/step_052-001.h5'

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from matplotlib.colors import TwoSlopeNorm
### Illapel ###

# f = h5py.File(data,'r')

# Tr = np.array(f['ParameterSets']['risetime']).T
# Vr = np.array(f['ParameterSets']['rupturevelocity']).T
# U = np.sqrt(np.array(f['ParameterSets']['dipslip']).T**2 + np.array(f['ParameterSets']['strikeslip']).T**2)

# T_Vr = np.concatenate((Tr,Vr))
# corr_T_Vr = np.corrcoef(T_Vr)
# slip_velocity = U/Tr


# Vr_slip_velocity = np.concatenate((Vr,slip_velocity))
# corr_T_Vr = np.corrcoef(Vr_slip_velocity)
# corr = np.diagonal(corr_T_Vr[:170,170:340])
# fig,ax =plt.subplots()
# parameter = np.flip(corr.reshape(10,17,order='F'),axis=0)
# im0 = ax.imshow(parameter,cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
# fig.colorbar(im0, ax=ax,shrink=0.65,label='Correlation')

# mean_slip_velocity = np.mean(slip_velocity,axis=1)
# mean_Vr  = np.mean(Vr,axis=1)
#plt.scatter(mean_Vr ,mean_slip_velocity)
# plt.scatter(slip_velocity,Vr)

# plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)

# plt.xlabel('Vr')
# plt.ylabel('Slip Velocity')
# plt.title('Illapel')
# plt.show()
# plt.close()

# plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Tr')
# plt.title('Illapel')
# plt.show()
# plt.close()
# plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('U')
# plt.title('Illapel')
# plt.show()
# plt.close()
# plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Tr')
# plt.ylabel('U')
# plt.title('Illapel')
# plt.show()

### Iquique ###

# data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Iquique/model/kinematic/step_056.h5'

# f = h5py.File(data,'r')
# Np= 132
# nramp = 3 
# d = np.array(f['Sample Set']).T
# Tr = d[2*Np+nramp:3*Np+nramp,:]
# Vr = d[3*Np+nramp:4*Np+nramp,:]
# U_dip = d[:Np,:]
# U_stk = d[Np:2*Np,:]
# U = np.sqrt(U_dip**2 + U_stk**2)
# slip_velocity = U/Tr


# plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Slip Velocity')
# plt.title('Iquique')
# plt.show()
# plt.close()

# plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Tr')
# plt.title('Iquique')
# plt.show()
# plt.close()
# plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('U')
# plt.title('Iquique')
# plt.show()
# plt.close()
# plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Tr')
# plt.ylabel('U')
# plt.title('Iquique')
# plt.show()

# ### Gorkha ###
# data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Gorkha/model/kinematic/step_012.h5'
# f = h5py.File(data,'r')

# Np= 162
# nramp = 0 
# d = np.array(f['Sample Set']).T
# Tr = d[2*Np+nramp:3*Np+nramp,:]
# Vr = d[3*Np+nramp:4*Np+nramp,:]
# U_dip = d[:Np+nramp,:]
# U_stk = d[Np+nramp:2*Np+nramp,:]
# U = np.sqrt(U_dip**2 + U_stk**2)


# slip_velocity = U/Tr

# plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Slip Velocity')
# plt.title('Gorkha')
# plt.show()
# plt.close()

# plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Tr')
# plt.title('Gorkha')
# plt.show()
# plt.close()
# plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('U')
# plt.title('Gorkha')
# plt.show()
# plt.close()
# plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Tr')
# plt.ylabel('U')
# plt.title('Gorkha')
# plt.show()


# ### Tohoku ###


def model_args(name,ncols,nrows,patchs,nparameters,ramps):
    config = dict()
    for i,n in enumerate(name):
        config[n] = {}
        config[n]['ncols'] = ncols[i]
        config[n]['nrows'] = nrows[i]
        config[n]['patch'] = patchs[i]
        config[n]['nparams'] = nparameters[i] 
        config[n]['ramps'] = ramps[i]
            
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

xlimits = [(0,8),(0,3.5),(0,5),(0,4.5),(0,5)]
ylimits = [(0,3.5),(0,1),(0,4),(0,5.5),(0,3.5)]

models = model_args(names,ncols,nrows,patchs,nparams,ramps)
samples = 100
plt.ioff()


def scale_var(x):
    return (x-np.min(x))/(np.max(x) - np.min(x))
'''
 Correlation between U and Slip Rate 
'''
'''
for th in [0,5,10,15,20,25,30,35,40]:
    for j,name in enumerate(list(models.keys())):
    # for name in ['Tohoku']:
        mean_data = f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv'
        df = pd.read_csv(mean_data)
        
        
        data = f'INPUT/{name}/model/kinematic/{samples}_samples/bin_data/{name}_kinematic_n_{samples}.dat'
        ncols, nrows = models[name]['ncols'],models[name]['nrows']
        
        Slip = np.flip(df['Slip'].values.reshape((nrows,ncols),order='F'),axis=0).flatten()
        Np = ncols*nrows
        patch = models[name]['patch']
        data = np.fromfile(data,'double').reshape((models[name]['nparams'],samples))
        
        colors = np.linspace(0,1,Np)
        cmap = plt.colormaps.get_cmap('jet')
        
        xsrc = np.arange((1/2)*patch,ncols*patch,patch)
        ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
        Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
        xsrc_flat = Xsrc.flatten()
        ysrc_flat = Ysrc.flatten()
    
        Tr = data[2*Np:3*Np,:]
        Vr = data[3*Np:4*Np,:]
        U1 = data[:Np,:]
        U2 = data[Np:2*Np,:]
        U = np.sqrt(U1**2 + U2**2)
        slip_velocity = U/Tr
        slip_velocity[np.isnan(slip_velocity)] = 0.0
        U_mean = np.mean(U,axis=1)
        th_ratio = th/100
        ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]
        Vr_mean = np.mean(Vr,axis=1)
        slip_velocity_mean = np.mean(slip_velocity,axis=1)
    
        fig, axes = plt.subplots(1,2,figsize = (10,4),dpi = 800)
        for i in ids_greater:
            
            # axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
            data = {'vr':Vr[i,:],'utr':slip_velocity[i,:]}
            data = pd.DataFrame(data)
        
            sns.kdeplot(data,x='vr',y='utr',ax =axes[0],levels = [1-0.97, 1-0.68],linewidths=0.5)
        corrcoef = np.corrcoef(Vr[i,:],slip_velocity[i,:])
        axes[0].text(0.97, 0.90, r'$Corr =%s $'%(round(corrcoef[0][1],3)),
            verticalalignment='bottom', transform = axes[0].transAxes,horizontalalignment='right',
            color='green', fontsize=6.5,fontweight='bold')
        axes[0].scatter(Vr_mean[ids_greater],slip_velocity_mean[ids_greater], color='k' ,s=4,label = f'$U > {th_ratio}$'+'$U_{max}$')
        axes[0].set_title('$ U/T_r$ vs $V_r$ ',fontsize=9,fontweight='bold')
        axes[0].set_xlabel('$V_r$ (km/s)')
        axes[0].set_ylabel('$U/T_r$ (m/s)')
        axes[0].set_aspect('equal','box')
        axes[0].legend(loc= 'upper left',fontsize=6.5)
        axes[0].set_xlim(xlimits[j])
        axes[0].set_ylim(ylimits[j])
        corr = get_corr(slip_velocity,Vr,Np, (nrows,ncols))
        id_th_col = ids_greater//nrows 
        id_th_row = nrows - ids_greater%nrows -1
    
        im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
        #im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='viridis')
        
        cs =  axes[1].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
      
        axes[1].scatter(xsrc[id_th_col],ysrc[id_th_row],marker='.',color='k',s = 10,label = f'$U > {th_ratio}$'+'$U_{max}$')
        axes[1].set_title('Correlation per patch between $V_r$ and $U/T_r$',fontsize=9,fontweight='bold')
        #axes[1].set_title('$U/T_r$',fontsize=9,fontweight='bold')

        axes[1].set_xlabel('Along-strike (km)',fontsize=8)
        axes[1].set_ylabel('Along-dip (km)',fontsize=8)                         
        axes[1].set_aspect('equal', 'box')
        axes[1].legend(fontsize=6.5)
        fig.suptitle(f'{name} '+'Kinematic Model for $N_{samples}$'+f'$={samples}$',y=0.95,fontweight = 'bold')
        fig.colorbar(im,ax=axes[1],shrink=0.6,label= 'Correlation',orientation='horizontal')
        fig.tight_layout() 
        folder = os.path.join(os.getcwd(),'global_corr')
        os.makedirs(folder,exist_ok=True)
        file_name = os.path.join(folder,f'{name}_correlation_Vr_U_over_Tr_th_{th}_samples_{samples}_.png')
        fig.savefig(file_name,dpi=400)
        plt.close(fig)

'''
'''
Correlation between U and Duration

'''
'''
for th in [0,5,10,15,20,25,30,35,40]:
    for j,name in enumerate(list(models.keys())):
    # for name in ['Tohoku']:
        mean_data = f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv'
        df = pd.read_csv(mean_data)
        
        
        data = f'INPUT/{name}/model/kinematic/{samples}_samples/bin_data/{name}_kinematic_n_{samples}.dat'
        ncols, nrows = models[name]['ncols'],models[name]['nrows']
        
        Slip = np.flip(df['Slip'].values.reshape((nrows,ncols),order='F'),axis=0).flatten()
        Np = ncols*nrows
        patch = models[name]['patch']
        data = np.fromfile(data,'double').reshape((models[name]['nparams'],samples))
        
        colors = np.linspace(0,1,Np)
        cmap = plt.colormaps.get_cmap('jet')
        
        xsrc = np.arange((1/2)*patch,ncols*patch,patch)
        ysrc = np.arange(-(nrows-1/2)*patch,0, patch)
        Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
        xsrc_flat = Xsrc.flatten()
        ysrc_flat = Ysrc.flatten()
    
        Tr = data[2*Np:3*Np,:]
        Vr = data[3*Np:4*Np,:]
        U1 = data[:Np,:]
        U2 = data[Np:2*Np,:]
        U = np.sqrt(U1**2 + U2**2)
        slip_velocity = U/Tr
        slip_velocity[np.isnan(slip_velocity)] = 0.0
        U_mean = np.mean(U,axis=1)
        th_ratio = th/100
        ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]
        Tr_mean = np.mean(Tr,axis=1)
        slip_velocity_mean = np.mean(slip_velocity,axis=1)
    
        fig, axes = plt.subplots(1,2,figsize = (10,4),dpi = 800)
        for i in ids_greater:
            
            # axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
            data = {'tr':Tr[i,:],'u':U[i,:]}
            data = pd.DataFrame(data)
        
            sns.kdeplot(data,x='tr',y='u',ax =axes[0],levels = [1-0.97, 1-0.68],linewidths=0.5)
        corrcoef = np.corrcoef(Tr[i,:],U[i,:])
        axes[0].text(0.97, 0.90, r'$Corr =%s $'%(round(corrcoef[0][1],3)),
            verticalalignment='bottom', transform = axes[0].transAxes,horizontalalignment='right',
            color='green', fontsize=6.5,fontweight='bold')
        axes[0].scatter(Tr_mean[ids_greater],U_mean[ids_greater], color='k' ,s=4,label = f'$U > {th_ratio}$'+'$U_{max}$')
        axes[0].set_title('$ U$ vs $T_r$ ',fontsize=9,fontweight='bold')
        axes[0].set_xlabel('$T_r$ (s)')
        axes[0].set_ylabel('$U$ (m/s)')
        axes[0].set_aspect('equal','box')
        axes[0].legend(loc= 'upper left',fontsize=6.5)
        #axes[0].set_xlim(xlimits[j])
        #axes[0].set_ylim(ylimits[j])
        corr = get_corr(U,Tr,Np, (nrows,ncols))
        id_th_col = ids_greater//nrows 
        id_th_row = nrows - ids_greater%nrows -1
    
        im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
        #im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='viridis')
        
        cs =  axes[1].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
      
        axes[1].scatter(xsrc[id_th_col],ysrc[id_th_row],marker='.',color='k',s = 10,label = f'$U > {th_ratio}$'+'$U_{max}$')
        axes[1].set_title('Correlation per patch between $T_r$ and $U$',fontsize=9,fontweight='bold')
        #axes[1].set_title('$U/T_r$',fontsize=9,fontweight='bold')

        axes[1].set_xlabel('Along-strike (km)',fontsize=8)
        axes[1].set_ylabel('Along-dip (km)',fontsize=8)                         
        axes[1].set_aspect('equal', 'box')
        axes[1].legend(fontsize=6.5)
        fig.suptitle(f'{name} '+'Kinematic Model for $N_{samples}$'+f'$={samples}$',y=0.95,fontweight = 'bold')
        fig.colorbar(im,ax=axes[1],shrink=0.6,label= 'Correlation',orientation='horizontal')
        fig.tight_layout() 
        folder = os.path.join(os.getcwd(),'global_corr/U_Tr')
        os.makedirs(folder,exist_ok=True)
        file_name = os.path.join(folder,f'{name}_correlation_Tr_U_Tr_th_{th}_samples_{samples}_.png')
        fig.savefig(file_name,dpi=400)
        plt.close(fig)

'''
'''
Correlation between U and Rupture Speed
'''
import matplotlib as mpl


lims = {'Tohoku':(0.0,12.5),
              'Iquique':(0.25,4.0),
              'Illapel':(-0.5,5.0),
              'Pedernales':(1,2),
              'Gorkha':(1.25,-0.05)
              }
#or th in [0,5,10,15,20,25,30,35,40]:

#for th in [30]:
#for j,name in enumerate(list(models.keys())):
    
    
for j,name in enumerate(['Illapel']):
    fig, axes = plt.subplots(2,4,figsize = (14,7),dpi = 400)
    ths = [10,20,30]
    n_lines = len(ths)
    cmap = mpl.colormaps['winter']

    # Take colors at regular intervals spanning the colormap.
    colors = cmap(np.linspace(0, 1, n_lines))
    colors = np.flip(colors,axis=0)
    for m,th in enumerate(ths):
    # for name in ['Tohoku']:
        mean_data = f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv'
        df = pd.read_csv(mean_data)
        
        
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
        labels= [('$V_{r}$','$U$'),('$T_{r}$','$U$'),('$V_{r}$','$T_{r}$'),('$V_{r}$','$(U/T_{r})$')]
        units = [('(km/s)','(m)',),('(s)','(m)'),('(km/s)','(s)'),('(km/s)','(m/s)')]
        for k,pair in enumerate([(Vr,U),(Tr,U),(Vr,Tr),(Vr,slip_velocity)]):
            
                
            for i in ids_greater:
                
                # axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
                data = {'x':pair[0][i,:],'y':pair[1][i,:]}
                data = pd.DataFrame(data)
            
                sns.kdeplot(data,x='x',y='y',ax =axes[k//2][k%2],levels = [1-0.68],linewidths=0.5,color=colors[m])
            #corrcoef = np.corrcoef(pair[0][i,:],pair[1][i,:])
            x_mean = np.mean(pair[0],axis=1)
            y_mean = np.mean(pair[1],axis=1)
            corrcoef = np.corrcoef(x_mean[ids_greater],y_mean[ids_greater])

            #axes[k][0].text(0.97, 0.90, r'$Corr =%s $'%(round(corrcoef[0][1],3)),
            #    verticalalignment='bottom', transform = axes[k][0].transAxes,horizontalalignment='right',
            #    color='green', fontsize=6.5,fontweight='bold')
            
            axes[k//2][k%2].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m],alpha = alphas[m],s=3,label =r'$Corr =%s $'%(round(corrcoef[0][1],3)))

            #axes[k][0].scatter(x_mean[ids_greater],y_mean[ids_greater], color=colors[m] ,alpha = alphas[m],s=3,label =r'$Corr =%s $'%(round(corrcoef[0][1],3)))
            axes[k//2][k%2].set_title(f'{labels[k][1]} vs {labels[k][0]} ',fontsize=8.75,fontweight='bold')
            axes[k//2][k%2].set_xlabel(f'{labels[k][0]} {units[k][0]}')
            axes[k//2][k%2].set_ylabel(f'{labels[k][1]} {units[k][1]}')
            #axes[k][0].set_aspect('equal','box')
            if k==0 or k==1:
                
                axes[k//2][k%2].axhline(y=th_ratio*np.max(U_mean),linestyle=lns[m],linewidth = lws[m],color='k')
                axes[k//2][k%2].text(lims[name][k],(th_ratio + 0.015)*np.max(U_mean),f'${th_ratio}$'+'$U_{max}$',fontsize=5)
            axes[k//2][k%2].legend(loc= 'upper right',fontsize=5)
            #axes[k][0].set_xlim([-0.1,1.1])
            #axes[k][0].set_ylim([-0.1,1.1])
            
            
        
        
        
            corr = get_corr(pair[0],pair[1],Np, (nrows,ncols))
            id_th_col = ids_greater//nrows 
            id_th_row = nrows - ids_greater%nrows -1
            im = axes[k//2][k%2+2].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
            #im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='viridis')
            #cs =  axes[k//2][k%2+2].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=[0.2,0.5,0.8,1.25],colors='black',linestyles='dashed')
            cs =  axes[k//2][k%2+2].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=[0.25,0.7,1.25],colors='black',linestyles='dashed')

            #cs =  axes[k][1].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels=np.array(ths)*np.max(U_mean)/100,linewidths=0.4,colors=colors,alpha=alphas,linestyles='dashed')
          
            #axes[k][1].scatter(xsrc[id_th_col],ysrc[id_th_row],marker='.',color='k',s = 10,label = f'$U > {th_ratio}$'+'$U_{max}$')
            axes[k//2][k%2+2].set_title(f'Correlation per patch between {labels[k][1]} and {labels[k][0]}',fontsize=8.75,fontweight='bold')
            #axes[1].set_title('$U/T_r$',fontsize=9,fontweight='bold')
    
            axes[k//2][k%2+2].set_xlabel('Along-strike (km)',fontsize=8)
            axes[k//2][k%2+2].set_ylabel('Along-dip (km)',fontsize=8)                         
            if name=='Illapel':
                
                axes[1][1].set_ylim([0,2.25])
                
            axes[k//2][k%2+2].set_aspect('equal', 'box')
            if m==0:
                if name!= 'Iquique':
                    fig.colorbar(im,ax=axes[k//2][k%2+2],shrink=0.5,pad = 0.2,label= 'Correlation',orientation='horizontal')
                else:
                    fig.colorbar(im,ax=axes[k//2][k%2+2],shrink=0.5,pad=0.2,label= 'Correlation',orientation='horizontal')

                    
            
            
    fig.suptitle(f'{name} '+'Kinematic Model for $N_{samples}$'+f'$={samples}$',y=0.98,fontweight = 'bold')
    print(name)
    if name!= 'Iquique':
        
        fig.tight_layout() 
    else:
        plt.subplots_adjust(wspace = 0.4,hspace=0.4)
            
            
    folder = os.path.join(os.getcwd(),'global_corr/definite_figs')
    os.makedirs(folder,exist_ok=True)
    file_name = os.path.join(folder,f'adjusted_correlation_{name}_{th}_samples_{samples}_.png')
    fig.savefig(file_name,dpi=400)
    plt.close(fig)
