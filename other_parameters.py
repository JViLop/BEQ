# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 18:14:15 2024

@author: joanv
"""

  
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


def model_args(name,ncols,nrows,patchs,nparameters):
    config = dict()
    for i,n in enumerate(name):
        config[n] = {}
        config[n]['ncols'] = ncols[i]
        config[n]['nrows'] = nrows[i]
        config[n]['patch'] = patchs[i]
        config[n]['nparams'] = nparameters[i] 
            
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

names = ['Tohoku','Iquique','Illapel','Gorkha']
nrows = [9,11,10,9]
ncols = [24,12,17,18]
patchs = [29,17,18,10]
nparams = [866,533,682,650]

xlimits = [(0,8),(0,3.5),(0,5),(0,4.5)]
ylimits = [(0,3.5),(0,1),(0,4),(0,5.5)]
th = 10
models = model_args(names,ncols,nrows,patchs,nparams)
samples = 100


for j,name in enumerate(list(models.keys())):
    
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
    U_mean = np.mean(U,axis=1)
    th_ratio = th/100
    ids_greater = np.where(U_mean>th_ratio*np.max(U_mean))[0]
    Vr_mean = np.mean(Vr,axis=1)
    slip_velocity_mean = np.mean(slip_velocity,axis=1)
    
    fig, axes = plt.subplots(1,2,figsize = (10,4),dpi = 800)
    for i in ids_greater:
        
        axes[0].scatter(Vr[i,:],slip_velocity[i,:],s=2.5,alpha= 0.2,color=cmap(colors[i]))
   
    axes[0].scatter(Vr_mean[ids_greater],slip_velocity_mean[ids_greater], color='k' ,s=5,label = f'$U > {th_ratio}$'+'$U_{max}$')
    axes[0].set_title('$ U/T_r$ vs $V_r$ ',fontsize=9,fontweight='bold')
    axes[0].set_xlabel('$V_r$ (km/s)')
    axes[0].set_ylabel('$U/T_r$ (m/s)')
    axes[0].set_aspect('equal','box')
    axes[0].set_xlim(xlimits[j])
    axes[0].set_ylim(ylimits[j])
    axes[0].legend(fontsize=6.5)
    corr = get_corr(slip_velocity,Vr,Np, (nrows,ncols))
    
    id_th_col = ids_greater//nrows 
    id_th_row = nrows - ids_greater%nrows -1

    im = axes[1].pcolormesh(xsrc,ysrc,corr.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
    cs =  axes[1].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
  
    axes[1].scatter(xsrc[id_th_col],ysrc[id_th_row],marker='.',color='k',s = 10,label = f'$U > {th_ratio}$'+'$U_{max}$')
    axes[1].set_title('Correlation per patch between $V_r$ and $U/T_r$',fontsize=9,fontweight='bold')
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
    



    '''
    
    
    
    # d_data = 'correlation_Tohoku/Tohoku_displacement_nsamples_100.h5'
    # s_data = 'correlation_Tohoku/Tohoku_Stress_change_nsamples_100.h5'
    
    # D = h5py.File(d_data ,'r')
    # S = h5py.File(s_data ,'r')
    
    # nS = S['Stress change'].shape[1]//3
    
    # S_dip = np.array(S['Stress change']).T[2*nS:,:]
    
    # nD = D['displacement'].shape[1]//3
    
    # Dx = np.array(D['displacement']).T[:nD,:]
    # Dy = np.array(D['displacement']).T[nD:2*nD,:]
    # Dz = np.array(D['displacement']).T[2*nD:,:]
    
    # D = np.sqrt(Dx**2 + Dy**2 + Dz**2)
    # D_above = D.reshape(nrows + 4,ncols+ 12,100)[2:-2,6:-6,:]
    # D_above = D_above.reshape(Np,100)
    
    c1 = get_corr(Vr,Tr,Np, (nrows,ncols))
    c2 = get_corr(Vr,U,Np, (nrows,ncols))
    c3 = get_corr(slip_velocity,Vr,Np, (nrows,ncols))
    c4 = get_corr(Tr,U,Np, (nrows,ncols))
    # c5 = get_corr(Vr,S_dip)
    # c6 = get_corr(Vr,D_above)
    
    
    
    
    names = [['Vr','Tr'],['Vr','U'],['Vr','U/Tr'],['Tr','U'],['Vr','$\Delta\sigma_{dip}$'],['Vr','$|d|$']]
    # fig,ax =plt.subplots(3,2,figsize =(12,7.5),dpi = 900)
    # for i,c in enumerate([c1,c2,c3,c4,c5,c6]):
    #     im = ax[i//2][i%2].pcolormesh(xsrc,ysrc,c.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
    #     cs =  ax[i//2][i%2].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 6,linewidths=0.5,colors='black')
    #     x  = name[i][0]
    #     y  = name[i][1]
    #     ax[i//2][i%2].set_title(f'Correlation between {x} and {y}',fontsize=10,fontweight='bold')
    #     ax[i//2][i%2].set_xlabel('Along-strike (km)',fontsize=8)
    #     ax[i//2][i%2].set_ylabel('Along-dip (km)',fontsize=8)                         
    #     ax[i//2][i%2].set_aspect('equal', 'box')
        
    #     fig.colorbar(im,ax=ax[i//2][i%2],shrink=0.7,pad=0.05,label= 'Correlation',orientation='vertical')
    
    # fig.savefig('correlation_Tohoku/correlation.png')
    names = [['Vr','Tr'],['Vr','U/Tr']]
    fig,axes =plt.subplots(2,1,figsize =(5,4.5),dpi = 900)
    for i,c in enumerate([c1,c3]):
        im = axes[i].pcolormesh(xsrc,ysrc,c.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
        cs =  axes[i].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
        x  = names[i][0]
        y  = names[i][1]
        axes[i].set_title(f'{name} Correlation between {x} and {y}',fontsize=9.5,fontweight='bold')
        axes[i].set_xlabel('Along-strike (km)',fontsize=8)
        axes[i].set_ylabel('Along-dip (km)',fontsize=8)                         
        axes[i].set_aspect('equal', 'box')
        
        fig.colorbar(im,ax=axes[i],shrink=0.7,pad=0.05,label= 'Correlation',orientation='vertical')
    
    
    # fig.savefig(f'correlation_{name}/correlation.png')
    
    
    
    fig,axes =plt.subplots(3,1,figsize =(5,6),dpi = 900)
    Tr_mean,Vr_mean,U_Tr_mean = get_mean(Tr,Np, (nrows,ncols)),get_mean(Vr,Np, (nrows,ncols)),get_mean(slip_velocity,Np, (nrows,ncols))
    labels = ['Tr','Vr','U/Tr']
    units = ['s','km/s','m/s']
    for i,c in enumerate([Tr_mean,Vr_mean,U_Tr_mean]):
        im = axes[i].pcolormesh(xsrc,ysrc,c.reshape(nrows,ncols),edgecolors='none',cmap='viridis',alpha=0.85)
        # cs =  axes[i].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
        axes[i].set_title(f'{name} Mean {labels[i]}',fontsize=9.5,fontweight='bold')
        axes[i].set_xlabel('Along-strike (km)',fontsize=8)
        axes[i].set_ylabel('Along-dip (km)',fontsize=8)                         
        axes[i].set_aspect('equal', 'box')
        fig.tight_layout()
        fig.colorbar(im,ax=axes[i],shrink=0.7,pad=0.05,label= f'{labels[i]}' + f' ({units[i]})',orientation='vertical')
    
    
    
    #### Correlation of max. slip rate and rupture velocity
    ids_max_slip_velocity = np.argmax(slip_velocity,axis=0)
    
    vals, counts = np.unique(ids_max_slip_velocity, return_counts=True)
    
    #find mode
    mode = np.argwhere(counts == np.max(counts))
    
    id_max_U_Tr = vals[mode][0][0] 
    
    def corr_given_patch(x,y,id_patch = 1,patches =Np,geom = (nrows,ncols)):
        z = np.concatenate((x,y))
        corr = np.corrcoef(z)
        corr_patch = corr[id_patch,Np:2*Np]
        flipped_corr = np.flip(corr_patch.reshape(geom[0],geom[1],order='F'),axis=0).flatten()
        return flipped_corr
    
    max_U_Tr_corr = corr_given_patch(slip_velocity,Vr,id_max_U_Tr)
    n_max_col = id_max_U_Tr//nrows
    n_max_row = nrows - id_max_U_Tr%nrows - 1
    c = max_U_Tr_corr
    fig,ax =plt.subplots(figsize =(5,6),dpi = 900)
    
    im = ax.pcolormesh(xsrc,ysrc,c.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
    # cs =  axes[i].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
    ax.scatter(xsrc[n_max_col],ysrc[n_max_row],marker = '.',c='k')
    ax.set_title(f'{name} Correlation between max. U/Tr and Vr',fontsize=9.5,fontweight='bold')
    ax.set_xlabel('Along-strike (km)',fontsize=8)
    ax.set_ylabel('Along-dip (km)',fontsize=8)                         
    ax.set_aspect('equal', 'box')
    fig.colorbar(im,ax=ax,shrink=0.25,pad=0.05,label= 'Correlation',orientation='vertical')
    
    #### Correlation of min. rupture velocity and slip rate
    ids_min_Vr = np.argmin(Vr,axis=0)
    
    vals, counts = np.unique(ids_min_Vr, return_counts=True)
    
    #find mode
    mode = np.argwhere(counts == np.max(counts))
    
    id_min_Vr = vals[mode][0][0] 
    
    
    min_Vr_corr = corr_given_patch(Vr, slip_velocity,id_min_Vr)
    n_min_col = id_min_Vr//nrows
    n_min_row = nrows - id_min_Vr%nrows - 1
    c = min_Vr_corr
    fig,ax =plt.subplots(figsize =(5,6),dpi = 900)
    
    im = ax.pcolormesh(xsrc,ysrc,c.reshape(nrows,ncols),edgecolors='none',cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
    # cs =  axes[i].contour(Xsrc,Ysrc,Slip.reshape(nrows,ncols),levels = 7,linewidths=0.4,colors='black')
    ax.scatter(xsrc[n_min_col],ysrc[n_min_row],marker = '.',c='k')
    ax.set_title(f'{name} Correlation between min. Vr and slip velocity',fontsize=9.5,fontweight='bold')
    ax.set_xlabel('Along-strike (km)',fontsize=8)
    ax.set_ylabel('Along-dip (km)',fontsize=8)                         
    ax.set_aspect('equal', 'box')
    fig.colorbar(im,ax=ax,shrink=0.25,pad=0.05,label= 'Correlation',orientation='vertical')
    
    
    


# plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Slip Velocity')
# plt.title('Tohoku')
# plt.show()
# plt.close()

# plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('Tr')
# plt.title('Tohoku')
# plt.show()

# plt.close()
# plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Vr')
# plt.ylabel('U')
# plt.title('Tohoku')
# plt.show()
# plt.close()

# plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
# plt.xlabel('Tr')
# plt.ylabel('U')
# plt.title('Tohoku')
# plt.show()

### Tohoku ###
'''