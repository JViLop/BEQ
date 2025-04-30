# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 10:26:54 2024

@author: joanv
"""
import re
import obspy as obs
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
import sys
from scipy.interpolate import RegularGridInterpolator
import matplotlib as mpl
import h5py

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']

nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patches = [29,17,18,10,15]
magnitudes = [9.0,8.1,8.3,7.8,7.8]
geoms = [(nrows[i],ncols[i]) for i in range(len(names))]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,90,90,107,360-99]
ramps = [0,3,0,0,9]
factors =[3,3,3,2,2]
strikes = [194,-13.58,4,293,27.05]
scales = [1.0,0.8e-1,1e-1,1.25e-1,0.70e-1]
shrinks = [0.2,0.4,0.4,0.3,0.4]
hspaces = [-0.5,0.3,0.3,0.3,0.1]
def model_dict(names,geoms,patches,magnitudes,arrow_sizes,nparams,rakes,nramps,factors,strikes,scales,shrinks,hspaces):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i]
        model[name]['patch'] = patches[i]
        model[name]['magnitude'] = magnitudes[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i]
        model[name]['nramp'] = nramps[i]
        model[name]['factor'] = factors[i]
        model[name]['factor'] = factors[i]
        model[name]['strike'] = strikes[i]
        model[name]['scale'] = scales[i]
        model[name]['shrink'] = shrinks[i]
        model[name]['hspace'] = hspaces[i]
    return model



models = model_dict(names,geoms,patches,magnitudes,arrow_sizes,nparams,rakes,ramps,factors,strikes,scales,shrinks,hspaces)

def set_stn(c,factor,patch,nrows,ncols,control=0):
      
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn
  
def proj_ysrc_coords(patch,dip):
      proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
      proj_ysrc = np.zeros_like(proj_dysrc)
      for i in range(len(proj_ysrc)):
          proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
      ysrc = np.flip(proj_ysrc)
        
      return ysrc


    

def max_disp_x(dx,dy,dz,patch,nrows,ncols):
    disp = np.sqrt(dx**2+dy**2+dz**2) 
    x,y = set_stn(1,1,patch,nrows,ncols,control=0)
    max_d_patch = np.argmin(abs(disp-np.max(disp)))
    val = max(disp)
    col0_max_d,row0_max_d= max_d_patch%ncols, max_d_patch//ncols
    return col0_max_d, row0_max_d,val
    

def max_disp_t(max_x):
    
    abs_max = np.max(max_x[:,2])
    max_d_patch = np.argmin(abs(max_x[:,2]-abs_max))
    col_abs_max_d = int(max_x[max_d_patch,0])
    row_abs_max_d = int(max_x[max_d_patch,1])
    return col_abs_max_d,row_abs_max_d

names = ['Tohoku','Illapel','Iquique','Pedernales','Gorkha']
ids = 'abcde'

'''
# These lines create displacement time-series along trench-normal cross-section
fig,axes = plt.subplots(len(names),3,figsize=(15,12),dpi=1200,sharex=True)
for n,name in enumerate(names):
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    strike = models[name]['strike']*np.pi/180
    magnitude = models[name]['magnitude']
    folder = name
    nmodels  = 100
    comp = 'N'
    npoints = 256
    
    N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
    E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
    Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
    
    strike = strike - 90.0*np.pi/180.0
    
    X =  -N*np.sin(strike) + E*np.cos(strike)      
    Y =  N*np.cos(strike) + E*np.sin(strike) 
    
    
    mean_X= np.mean(X,axis=0)
    mean_Y= np.mean(Y,axis=0)
    mean_Z= np.mean(Z,axis=0)
    
    std_X = np.zeros_like(mean_X)
    std_Y = np.zeros_like(mean_Y)
    std_Z = np.zeros_like(mean_Z)
    
    for t in range(npoints):
        std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
        std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
        std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
    
    
    
    max_info = np.zeros((npoints,3))
    for t in range(npoints):
        dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
        col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
        max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
        
        
    col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
        
    
    mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
    mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
    mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
    
    mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
    mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
    mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)
    
    
    n_lines = nrows
    cmap = mpl.colormaps['winter']
    
    # Take colors at regular intervals spanning the colormap.
    colors = cmap(np.linspace(0, 1, n_lines))
    
    mean_disp = {'x':mean_X_mesh,'y':mean_Y_mesh,'z':mean_Z_mesh} 
    sigma_disp = {'x':mean_sigma_X_mesh,'y':mean_sigma_Y_mesh,'z':mean_sigma_Z_mesh}
    y_labels = ['Along-strike displacement (m)','Trench-normal displacement (m)','Vertical displacement (m)']
    y_labels = ['$d_{x}$ (m)','$d_{y}$ (m)','$d_{z}$ (m)']
    t = np.arange(0,npoints)
    for i,comp in enumerate(['x','y','z']):
        for j in range(nrows):
            upper_lim = mean_disp[comp][j,col_abs_max_d,:] +  sigma_disp[comp][j,col_abs_max_d,:]
            lower_lim = mean_disp[comp][j,col_abs_max_d,:] -  sigma_disp[comp][j,col_abs_max_d,:]
            axes[n][i].plot(t,mean_disp[comp][j,col_abs_max_d,:],lw=0.4,color=colors[j])
            axes[n][i].fill_between(t,lower_lim,upper_lim,alpha=0.2,facecolor=colors[j])    
            axes[n][0].set_title('$M_w%s$'%(magnitude)+f' {name}',fontweight='bold',loc='left',fontsize=14)    
            axes[n][i].set_ylabel(y_labels[i],fontsize=12)
            axes[n][0].text(-0.2, 1.1, ids[n], fontweight='bold',fontsize=20,transform=axes[n][0].transAxes)
            fig.align_ylabels(fig.get_axes())
            plt.subplots_adjust(wspace=0.2,hspace=0.325)
fig.supxlabel('Time from mainshock (s)',y=0.0725,fontsize=12)            
fig.savefig(f'all_event_evolution_peak_displacement_{nmodels}.pdf')







# These lines create snapshots of 3D displacement
from mpl_toolkits.axes_grid1 import make_axes_locatable
nmodels  = 100

nepochs = 6 #20
fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20),dpi=1200)

for n,name in enumerate(names):
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    strike = models[name]['strike']
    magnitude = models[name]['magnitude']
    scale =  models[name]['scale']
    shrink = models[name]['shrink']
    hspace = models[name]['hspace']
    folder = name
    nmodels  = 100
    comp = 'N'
    npoints = 256
    x = np.arange(patch/2,ncols*patch,patch)
    y = np.arange(patch/2,nrows*patch,patch)
    y = -np.flip(y,axis=0)
    
    Xr,Yr = np.meshgrid(x,y)
    
    strike = strike*np.pi/180

    nsamples = 100
    
    nstations = nrows*ncols
    npoints = 256
    shape = (nmodels,nstations,npoints,3)
    data = np.zeros(shape)
    
    N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
    E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
    Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
    
    strike = strike - 90.0*np.pi/180.0
    
    X =  -N*np.sin(strike) + E*np.cos(strike)      
    Y =  N*np.cos(strike) + E*np.sin(strike) 
    
    
    mean_X= np.mean(X,axis=0)
    mean_Y= np.mean(Y,axis=0)
    mean_Z= np.mean(Z,axis=0)
    
    std_X = np.zeros_like(mean_X)
    std_Y = np.zeros_like(mean_Y)
    std_Z = np.zeros_like(mean_Z)
    
    for t in range(npoints):
        std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
        std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
        std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
    
    
    
    max_info = np.zeros((npoints,3))
    for t in range(npoints):
        dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
        col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
        max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
        
        
    col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
        
    
    mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
    mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
    mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
    
    mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
    mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
    mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)

    

    for k,epoch in enumerate(list(np.logspace(3,np.log2(256),num=nepochs,base=2))):
    #for k,epoch in enumerate(list(np.linspace(0,75,nepochs))):
             
        if epoch==256:
            epoch-=1
        
        Ux = mean_X_mesh[:,:,int(epoch)] 
        Uy = mean_Y_mesh[:,:,int(epoch)] 
        Uz = mean_Z_mesh[:,:,int(epoch)] 
        
        Uxf = mean_X_mesh[:,:,-1] 
        Uyf = mean_Y_mesh[:,:,-1] 
        Uzf = mean_Z_mesh[:,:,-1] 
        U = np.sqrt(Ux**2 + Uy**2)
        
        #cax1 = axes[n][k].pcolormesh(x, y, Uz,edgecolors='k',linewidth=0.25,vmin = min(np.min(Uzf),-np.max(Uzf)),vmax = max(-np.min(Uzf),np.max(Uzf)),cmap='bwr')
        cax1 = axes[n][k].pcolormesh(x, y, Uz,linewidth=0.25,vmin = min(np.min(Uzf),-np.max(Uzf)),vmax = max(-np.min(Uzf),np.max(Uzf)),cmap='bwr')

        q = axes[n][k].quiver(Xr,Yr,Ux.flatten(),Uy.flatten(),scale=scale,scale_units ='x', units='width',width=0.005,headwidth=2.5,headlength=4.0)




        #props = dict(boxstyle='round', edgecolor='black',facecolor='white', alpha=0.5)
        #axes[n][k].text(0.02, 0.12, '             ', transform=axes[n][k].transAxes, fontsize=18,verticalalignment='top', bbox=props)

        axes[n][k].quiverkey(q, 0.86,1.05 , U=np.ceil(U.max()/2),
                      label=' {}m'.format(int(np.ceil(U.max()/2))), labelpos='N',labelsep = 0.05,angle=0,transform=axes[n][k].transAxes)
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

  
        

        axes[0][k].text(0.05, 1.5, '$t = %s s$'%(int(round(epoch,0))), transform=axes[0][k].transAxes,fontsize = 15,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        axes[n][k].set_aspect('equal','box')
        if k==0:
            axes[n][k].set_title('$M_w%s$'%(magnitude)+f' {name}', loc='left',fontweight='bold',fontsize=15)
            axes[n][k].set_ylabel('Trench-normal distance (km)',fontsize=10)
            axes[n][k].set_xlabel('Along-strike distance (km)',fontsize=10)
            axes[n][0].text(-0.15, 1.18, ids[n], fontweight='bold',fontsize=28,transform=axes[n][0].transAxes)
        else:
            axes[n][k].set_xticklabels([])
            axes[n][k].set_yticklabels([])
            axes[n][k].set_xlabel('')
            axes[n][k].set_ylabel('')
            axes[n][k].get_xaxis().set_ticks([])
            axes[n][k].get_yaxis().set_ticks([])
        #fig.align_ylabels(fig.get_axes())
        if k==nepochs -1:
            divider = make_axes_locatable(axes[n][k])
            cax = divider.append_axes("right", size="2%", pad=0.075)

            cbar = fig.colorbar(cax1,cax = cax,shrink=shrink,label='$d_z$ (m)')
            cbar.ax.tick_params(labelsize=10) 
        
        fig.subplots_adjust(wspace=0.1,hspace=hspace)
 
fig.savefig(f'snapshot_xyz_many_{nepochs}.pdf')


'''
#Error in x 


from mpl_toolkits.axes_grid1 import make_axes_locatable
nmodels  = 100

nepochs = 12
fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20))

for n,name in enumerate(names):
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    strike = models[name]['strike']
    magnitude = models[name]['magnitude']
    scale =  models[name]['scale']
    shrink = models[name]['shrink']
    hspace = models[name]['hspace']
    folder = name
    nmodels  = 100
    comp = 'N'
    npoints = 256
    x = np.arange(patch/2,ncols*patch,patch)
    y = np.arange(patch/2,nrows*patch,patch)
    y = -np.flip(y,axis=0)
    
    Xr,Yr = np.meshgrid(x,y)
    
    strike = strike*np.pi/180

    nsamples = 100
    
    nstations = nrows*ncols
    npoints = 256
    shape = (nmodels,nstations,npoints,3)
    data = np.zeros(shape)
    
    N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
    E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
    Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
    
    strike = strike - 90.0*np.pi/180.0
    
    X =  -N*np.sin(strike) + E*np.cos(strike)      
    Y =  N*np.cos(strike) + E*np.sin(strike) 
    
    
    mean_X= np.mean(X,axis=0)
    mean_Y= np.mean(Y,axis=0)
    mean_Z= np.mean(Z,axis=0)
    
    std_X = np.zeros_like(mean_X)
    std_Y = np.zeros_like(mean_Y)
    std_Z = np.zeros_like(mean_Z)
    
    for t in range(npoints):
        std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
        std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
        std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
    
    
    
    max_info = np.zeros((npoints,3))
    for t in range(npoints):
        dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
        col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
        max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
        
        
    col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
        
    
    mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
    mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
    mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
    
    mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
    mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
    mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)

    

    #for k,epoch in enumerate(list(np.logspace(1,np.log2(256),num=nepochs,base=2))):
    for k,epoch in enumerate(list(np.linspace(2,256,nepochs))):

    #for k,epoch in enumerate(list(np.linspace(0,75,nepochs))):
             
        if epoch==256:
            epoch-=1
        
        sigma_Ux = mean_sigma_X_mesh[:,:,int(epoch)] 

        
        sigma_Uxf = mean_sigma_X_mesh[:,:,-1] 

        
        cax1 = axes[n][k].pcolormesh(x, y,  sigma_Ux ,linewidth=0.25,vmin = 0,vmax = np.max(sigma_Uxf),cmap='rainbow')


        
        # minus sign in Y coordinate corrects for y-axis inversion
        

        #coord = q.XY + offsetXY*factor
        # width and height had to be formatted in row-like order (initially in column-like order)
#        ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
#                        width=dictionary['std_U_perp'][i]*factor,
#                        height=dictionary['std_U_parallel'][i]*factor,
#                        angle=0,alpha=0.5,fill=False,edgecolor='k')
#                for i in range(q.N)]
#        for e in ells:
#            ax.add_artist(e)

  
        

        axes[0][k].text(0.05, 1.5, '$t = %s s$'%(int(round(epoch,0))), transform=axes[0][k].transAxes,fontsize = 15,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        axes[n][k].set_aspect('equal','box')
        if k==0:
            axes[n][k].set_title('$M_w%s$'%(magnitude)+f' {name}', loc='left',fontweight='bold',fontsize=15)
            axes[n][k].set_ylabel('Trench-normal distance (km)',fontsize=9)
            axes[n][k].set_xlabel('Along-strike distance (km)',fontsize=9)
            axes[n][0].text(-0.15, 1.18, ids[n], fontweight='bold',fontsize=20,transform=axes[n][0].transAxes)

        else:
            axes[n][k].set_xticklabels([])
            axes[n][k].set_yticklabels([])
            axes[n][k].set_xlabel('')
            axes[n][k].set_ylabel('')
            axes[n][k].get_xaxis().set_ticks([])
            axes[n][k].get_yaxis().set_ticks([])
        #fig.align_ylabels(fig.get_axes())
        if k==nepochs -1:
            divider = make_axes_locatable(axes[n][k])
            cax = divider.append_axes("right", size="2%", pad=0.075)

            fig.colorbar(cax1,cax = cax,shrink=shrink,label='s(m)')
        
        fig.subplots_adjust(wspace=0.1,hspace=hspace)
 
fig.savefig(f'snapshot_sigma_x_{nepochs}.pdf')




#Error in y


from mpl_toolkits.axes_grid1 import make_axes_locatable
nmodels  = 100

nepochs = 12 #20
fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20))

for n,name in enumerate(names):
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    strike = models[name]['strike']
    magnitude = models[name]['magnitude']
    scale =  models[name]['scale']
    shrink = models[name]['shrink']
    hspace = models[name]['hspace']
    folder = name
    nmodels  = 100
    comp = 'N'
    npoints = 256
    x = np.arange(patch/2,ncols*patch,patch)
    y = np.arange(patch/2,nrows*patch,patch)
    y = -np.flip(y,axis=0)
    
    Xr,Yr = np.meshgrid(x,y)
    
    strike = strike*np.pi/180

    nsamples = 100
    
    nstations = nrows*ncols
    npoints = 256
    shape = (nmodels,nstations,npoints,3)
    data = np.zeros(shape)
    
    N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
    E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
    Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
    
    strike = strike - 90.0*np.pi/180.0
    
    X =  -N*np.sin(strike) + E*np.cos(strike)      
    Y =  N*np.cos(strike) + E*np.sin(strike) 
    
    
    mean_X= np.mean(X,axis=0)
    mean_Y= np.mean(Y,axis=0)
    mean_Z= np.mean(Z,axis=0)
    
    std_X = np.zeros_like(mean_X)
    std_Y = np.zeros_like(mean_Y)
    std_Z = np.zeros_like(mean_Z)
    
    for t in range(npoints):
        std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
        std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
        std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
    
    
    
    max_info = np.zeros((npoints,3))
    for t in range(npoints):
        dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
        col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
        max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
        
        
    col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
        
    
    mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
    mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
    mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
    
    mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
    mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
    mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)

    

    #for k,epoch in enumerate(list(np.logspace(4,np.log2(256),num=nepochs,base=2))):
    for k,epoch in enumerate(list(np.linspace(2,256,nepochs))):
    #for k,epoch in enumerate(list(np.linspace(0,75,nepochs))):
             
        if epoch==256:
            epoch-=1
        
        sigma_U = mean_sigma_Y_mesh[:,:,int(epoch)] 

        
        sigma_Uf = mean_sigma_Y_mesh[:,:,-1] 

        
        cax1 = axes[n][k].pcolormesh(x, y,  sigma_U ,linewidth=0.25,vmin = 0,vmax = np.max(sigma_Uf),cmap='rainbow')


        
        # minus sign in Y coordinate corrects for y-axis inversion
        

        #coord = q.XY + offsetXY*factor
        # width and height had to be formatted in row-like order (initially in column-like order)
#        ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
#                        width=dictionary['std_U_perp'][i]*factor,
#                        height=dictionary['std_U_parallel'][i]*factor,
#                        angle=0,alpha=0.5,fill=False,edgecolor='k')
#                for i in range(q.N)]
#        for e in ells:
#            ax.add_artist(e)

  
        

        axes[0][k].text(0.05, 1.5, '$t = %s s$'%(int(round(epoch,0))), transform=axes[0][k].transAxes,fontsize = 15,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        axes[n][k].set_aspect('equal','box')
        if k==0:
            axes[n][k].set_title('$M_w%s$'%(magnitude)+f' {name}', loc='left',fontweight='bold',fontsize=15)
            axes[n][k].set_ylabel('Trench-normal distance (km)',fontsize=9)
            axes[n][k].set_xlabel('Along-strike distance (km)',fontsize=9)
            axes[n][0].text(-0.15, 1.18, ids[n], fontweight='bold',fontsize=20,transform=axes[n][0].transAxes)

        else:
            axes[n][k].set_xticklabels([])
            axes[n][k].set_yticklabels([])
            axes[n][k].set_xlabel('')
            axes[n][k].set_ylabel('')
            axes[n][k].get_xaxis().set_ticks([])
            axes[n][k].get_yaxis().set_ticks([])
        #fig.align_ylabels(fig.get_axes())
        if k==nepochs -1:
            divider = make_axes_locatable(axes[n][k])
            cax = divider.append_axes("right", size="2%", pad=0.075)

            fig.colorbar(cax1,cax = cax,shrink=shrink,label='s(m)')
        
        fig.subplots_adjust(wspace=0.1,hspace=hspace)
 
fig.savefig(f'snapshot_sigma_y_{nepochs}.pdf')




#Error in Z


from mpl_toolkits.axes_grid1 import make_axes_locatable
nmodels  = 100

nepochs = 12 #20
fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20))

for n,name in enumerate(names):
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    strike = models[name]['strike']
    magnitude = models[name]['magnitude']
    scale =  models[name]['scale']
    shrink = models[name]['shrink']
    hspace = models[name]['hspace']
    folder = name
    nmodels  = 100
    npoints = 256
    x = np.arange(patch/2,ncols*patch,patch)
    y = np.arange(patch/2,nrows*patch,patch)
    y = -np.flip(y,axis=0)
    
    Xr,Yr = np.meshgrid(x,y)
    
    strike = strike*np.pi/180

    nsamples = 100
    
    nstations = nrows*ncols
    npoints = 256
    shape = (nmodels,nstations,npoints,3)
    data = np.zeros(shape)
    
    N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
    E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
    Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
    
    strike = strike - 90.0*np.pi/180.0
    
    X =  -N*np.sin(strike) + E*np.cos(strike)      
    Y =  N*np.cos(strike) + E*np.sin(strike) 
    
    
    mean_X= np.mean(X,axis=0)
    mean_Y= np.mean(Y,axis=0)
    mean_Z= np.mean(Z,axis=0)
    
    std_X = np.zeros_like(mean_X)
    std_Y = np.zeros_like(mean_Y)
    std_Z = np.zeros_like(mean_Z)
    
    for t in range(npoints):
        std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
        std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
        std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
    
    
    
    max_info = np.zeros((npoints,3))
    for t in range(npoints):
        dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
        col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
        max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
        
        
    col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
        
    
    mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
    mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
    mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
    
    mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
    mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
    mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)

    

    #for k,epoch in enumerate(list(np.logspace(4,np.log2(256),num=nepochs,base=2))):
    for k,epoch in enumerate(list(np.linspace(2,256,nepochs))):
    #for k,epoch in enumerate(list(np.linspace(0,75,nepochs))):
             
        if epoch==256:
            epoch-=1
        
        sigma_U = mean_sigma_Z_mesh[:,:,int(epoch)] 

        
        sigma_Uf = mean_sigma_Z_mesh[:,:,-1] 

        
        cax1 = axes[n][k].pcolormesh(x, y,  sigma_U ,linewidth=0.25,vmin = 0,vmax = np.max(sigma_Uf),cmap='rainbow')


        
        # minus sign in Y coordinate corrects for y-axis inversion
        

        #coord = q.XY + offsetXY*factor
        # width and height had to be formatted in row-like order (initially in column-like order)
#        ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
#                        width=dictionary['std_U_perp'][i]*factor,
#                        height=dictionary['std_U_parallel'][i]*factor,
#                        angle=0,alpha=0.5,fill=False,edgecolor='k')
#                for i in range(q.N)]
#        for e in ells:
#            ax.add_artist(e)

  
        

        axes[0][k].text(0.05, 1.5, '$t = %s s$'%(int(round(epoch,0))), transform=axes[0][k].transAxes,fontsize = 15,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))

        axes[n][k].set_aspect('equal','box')
        if k==0:
            axes[n][k].set_title('$M_w%s$'%(magnitude)+f' {name}', loc='left',fontweight='bold',fontsize=15)
            axes[n][k].set_ylabel('Trench-normal distance (km)',fontsize=9)
            axes[n][k].set_xlabel('Along-strike distance (km)',fontsize=9)
            axes[n][0].text(-0.15, 1.18, ids[n], fontweight='bold',fontsize=20,transform=axes[n][0].transAxes)

        else:
            axes[n][k].set_xticklabels([])
            axes[n][k].set_yticklabels([])
            axes[n][k].set_xlabel('')
            axes[n][k].set_ylabel('')
            axes[n][k].get_xaxis().set_ticks([])
            axes[n][k].get_yaxis().set_ticks([])
        #fig.align_ylabels(fig.get_axes())
        if k==nepochs -1:
            divider = make_axes_locatable(axes[n][k])
            cax = divider.append_axes("right", size="2%", pad=0.075)

            fig.colorbar(cax1,cax = cax,shrink=shrink,label='s(m)')
        
        fig.subplots_adjust(wspace=0.1,hspace=hspace)
 
 
fig.savefig(f'snapshot_sigma_z_{nepochs}.pdf')

'''

# Relative standard_deviation


#Error in comp


from mpl_toolkits.axes_grid1 import make_axes_locatable
nmodels  = 100

nepochs = 12
for comp in ['x','y','z']:
    fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20),dpi=400)
    
    for n,name in enumerate(names):
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        patch = models[name]['patch']
        strike = models[name]['strike']
        magnitude = models[name]['magnitude']
        scale =  models[name]['scale']
        shrink = models[name]['shrink']
        hspace = models[name]['hspace']
        folder = name
        nmodels  = 100
        npoints = 256
        x = np.arange(patch/2,ncols*patch,patch)
        y = np.arange(patch/2,nrows*patch,patch)
        y = -np.flip(y,axis=0)
        
        Xr,Yr = np.meshgrid(x,y)
        
        strike = strike*np.pi/180
    
        nsamples = 100
        
        nstations = nrows*ncols
        npoints = 256
        shape = (nmodels,nstations,npoints,3)
        data = np.zeros(shape)
        
        N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
        E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
        Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
        
        strike = strike - 90.0*np.pi/180.0
        
        X =  -N*np.sin(strike) + E*np.cos(strike)      
        Y =  N*np.cos(strike) + E*np.sin(strike) 
        
        
        mean_X= np.mean(X,axis=0)
        mean_Y= np.mean(Y,axis=0)
        mean_Z= np.mean(Z,axis=0)
        
        std_X = np.zeros_like(mean_X)
        std_Y = np.zeros_like(mean_Y)
        std_Z = np.zeros_like(mean_Z)
        
        for t in range(npoints):
            std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
            std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
            std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
        
        
        
        max_info = np.zeros((npoints,3))
        for t in range(npoints):
            dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
            col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
            max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
            
            
        col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
            
        
        mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
        mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
        mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
        
        mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
        mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
        mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)
    
        rel_sigma_X_mesh = 100*mean_sigma_X_mesh/abs(mean_X_mesh)
        rel_sigma_Y_mesh = 100*mean_sigma_Y_mesh/abs(mean_Y_mesh)
        rel_sigma_Z_mesh = 100*mean_sigma_Z_mesh/abs(mean_Z_mesh)
    
        #for k,epoch in enumerate(list(np.logspace(4,np.log2(256),num=nepochs,base=2))):
        for k,epoch in enumerate(list(np.linspace(2,256,nepochs))):
        #for k,epoch in enumerate(list(np.linspace(0,75,nepochs))):
                 
            if epoch==256:
                epoch-=1
            
            if comp=='x':
                
                rel_sigma_U = rel_sigma_X_mesh[:,:,int(epoch)] 
            elif comp=='y':
                rel_sigma_U = rel_sigma_Y_mesh[:,:,int(epoch)] 
            else:
                rel_sigma_U = rel_sigma_Z_mesh[:,:,int(epoch)] 
                
    
    
            
            cax1 = axes[n][k].pcolormesh(x, y,  rel_sigma_U,linewidth=0.25,vmin = 0,vmax = 100,cmap='Blues')
    
    
            
            # minus sign in Y coordinate corrects for y-axis inversion
            
    
            #coord = q.XY + offsetXY*factor
            # width and height had to be formatted in row-like order (initially in column-like order)
    #        ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
    #                        width=dictionary['std_U_perp'][i]*factor,
    #                        height=dictionary['std_U_parallel'][i]*factor,
    #                        angle=0,alpha=0.5,fill=False,edgecolor='k')
    #                for i in range(q.N)]
    #        for e in ells:
    #            ax.add_artist(e)
    
      
            
    
            axes[0][k].text(0.05, 1.5, '$t = %s s$'%(int(round(epoch,0))), transform=axes[0][k].transAxes,fontsize = 15,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    
            axes[n][k].set_aspect('equal','box')
            if k==0:
                axes[n][k].set_title('$M_w%s$'%(magnitude)+f' {name}', loc='left',fontweight='bold',fontsize=15)
                axes[n][k].set_ylabel('Trench-normal distance (km)',fontsize=9)
                axes[n][k].set_xlabel('Along-strike distance (km)',fontsize=9)
                axes[n][0].text(-0.15, 1.18, ids[n], fontweight='bold',fontsize=20,transform=axes[n][0].transAxes)

            else:
                axes[n][k].set_xticklabels([])
                axes[n][k].set_yticklabels([])
                axes[n][k].set_xlabel('')
                axes[n][k].set_ylabel('')
                axes[n][k].get_xaxis().set_ticks([])
                axes[n][k].get_yaxis().set_ticks([])
            #fig.align_ylabels(fig.get_axes())
            if k==nepochs -1:
                divider = make_axes_locatable(axes[n][k])
                cax = divider.append_axes("right", size="2%", pad=0.075)
    
                fig.colorbar(cax1,cax = cax,shrink=shrink,label='rel. error (%)',extend='max')
            
            fig.subplots_adjust(wspace=0.1,hspace=hspace)
     
    fig.savefig(f'snapshot_rel_sigma_{comp}_{nepochs}.pdf')
    


# Relative standard_deviation


#Error in comp


from mpl_toolkits.axes_grid1 import make_axes_locatable
nmodels  = 100

nepochs = 12
for comp in ['x','y','z']:
    fig, axes = plt.subplots(5,nepochs,figsize=(int(nepochs*4),20),dpi=400)
    
    for n,name in enumerate(names):
        nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
        patch = models[name]['patch']
        strike = models[name]['strike']
        magnitude = models[name]['magnitude']
        scale =  models[name]['scale']
        shrink = models[name]['shrink']
        hspace = models[name]['hspace']
        folder = name
        nmodels  = 100
        npoints = 256
        x = np.arange(patch/2,ncols*patch,patch)
        y = np.arange(patch/2,nrows*patch,patch)
        y = -np.flip(y,axis=0)
        
        Xr,Yr = np.meshgrid(x,y)
        
        strike = strike*np.pi/180
    
        nsamples = 100
        
        nstations = nrows*ncols
        npoints = 256
        shape = (nmodels,nstations,npoints,3)
        data = np.zeros(shape)
        
        N  = np.array(h5py.File(os.path.join(folder,f'{name}_N_nmodels_{nmodels}.h5'),'r')['N'])
        E = np.array(h5py.File(os.path.join(folder,f'{name}_E_nmodels_{nmodels}.h5'),'r')['E'])
        Z = np.array(h5py.File(os.path.join(folder,f'{name}_Z_nmodels_{nmodels}.h5'),'r')['Z'])
        
        strike = strike - 90.0*np.pi/180.0
        
        X =  -N*np.sin(strike) + E*np.cos(strike)      
        Y =  N*np.cos(strike) + E*np.sin(strike) 
        
        
        mean_X= np.mean(X,axis=0)
        mean_Y= np.mean(Y,axis=0)
        mean_Z= np.mean(Z,axis=0)
        
        std_X = np.zeros_like(mean_X)
        std_Y = np.zeros_like(mean_Y)
        std_Z = np.zeros_like(mean_Z)
        
        for t in range(npoints):
            std_X[:,t] = np.sqrt(np.diag(np.cov(X[:,:,t].T)))
            std_Y[:,t] = np.sqrt(np.diag(np.cov(Y[:,:,t].T)))
            std_Z[:,t] = np.sqrt(np.diag(np.cov(Z[:,:,t].T)))
        
        
        
        max_info = np.zeros((npoints,3))
        for t in range(npoints):
            dx,dy,dz = mean_X[:,t],mean_Y[:,t],mean_Z[:,t]
            col0_max_d,row0_max_d,val = max_disp_x(dx,dy,dz,patch,nrows,ncols)
            max_info[t,:] = np.array([col0_max_d,row0_max_d,val])
            
            
        col_abs_max_d,row_abs_max_d = max_disp_t(max_info)
            
        
        mean_X_mesh = mean_X.reshape(nrows,ncols,npoints)
        mean_Y_mesh = mean_Y.reshape(nrows,ncols,npoints)
        mean_Z_mesh = mean_Z.reshape(nrows,ncols,npoints)
        
        mean_sigma_X_mesh = std_X.reshape(nrows,ncols,npoints)
        mean_sigma_Y_mesh = std_Y.reshape(nrows,ncols,npoints)
        mean_sigma_Z_mesh = std_Z.reshape(nrows,ncols,npoints)
    
        rel_sigma_X_mesh = 100*mean_sigma_X_mesh/abs(mean_X_mesh)
        rel_sigma_Y_mesh = 100*mean_sigma_Y_mesh/abs(mean_Y_mesh)
        rel_sigma_Z_mesh = 100*mean_sigma_Z_mesh/abs(mean_Z_mesh)
    
        for k,epoch in enumerate(list(np.logspace(4,np.log2(256),num=nepochs,base=2))):
        #for k,epoch in enumerate(list(np.linspace(0,75,nepochs))):
                 
            if epoch==256:
                epoch-=1
            
            if comp=='x':
                
                U = mean_X_mesh[:,:,int(epoch)] 
                Ut = mean_X_mesh[:,:,-1]
            elif comp=='y':
                U = mean_Y_mesh[:,:,int(epoch)] 
                Ut = mean_Y_mesh[:,:,-1]
            else:
                U = mean_Z_mesh[:,:,int(epoch)] 
                Ut = mean_Z_mesh[:,:,-1]
                
                
    
    
            
            cax1 = axes[n][k].pcolormesh(x, y, U ,edgecolors='k',linewidth=0.25,vmin = min(np.min(Ut),-np.max(Ut)),vmax =  max(-np.min(Ut),np.max(Ut)),cmap='bwr')
    
    
            
            # minus sign in Y coordinate corrects for y-axis inversion
            
    
            #coord = q.XY + offsetXY*factor
            # width and height had to be formatted in row-like order (initially in column-like order)
    #        ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
    #                        width=dictionary['std_U_perp'][i]*factor,
    #                        height=dictionary['std_U_parallel'][i]*factor,
    #                        angle=0,alpha=0.5,fill=False,edgecolor='k')
    #                for i in range(q.N)]
    #        for e in ells:
    #            ax.add_artist(e)
    
      
            
    
            axes[0][k].text(0.05, 1.5, '$t = %s s$'%(int(round(epoch,0))), transform=axes[0][k].transAxes,fontsize = 15,bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    
            axes[n][k].set_aspect('equal','box')
            if k==0:
                axes[n][k].set_title('$M_w%s$'%(magnitude)+f' {name}', loc='left',fontweight='bold',fontsize=15)
                axes[n][k].set_ylabel('Trench-normal distance (km)',fontsize=9)
                axes[n][k].set_xlabel('Along-strike distance (km)',fontsize=9)
                axes[n][0].text(-0.15, 1.18, ids[n], fontweight='bold',fontsize=20,transform=axes[n][0].transAxes)

            else:
                axes[n][k].set_xticklabels([])
                axes[n][k].set_yticklabels([])
                axes[n][k].set_xlabel('')
                axes[n][k].set_ylabel('')
                axes[n][k].get_xaxis().set_ticks([])
                axes[n][k].get_yaxis().set_ticks([])
            #fig.align_ylabels(fig.get_axes())
            if k==nepochs -1:
                divider = make_axes_locatable(axes[n][k])
                cax = divider.append_axes("right", size="2%", pad=0.075)
    
                fig.colorbar(cax1,cax = cax,shrink=shrink,label='$ %s $(m)'%(comp))
            
            fig.subplots_adjust(wspace=0.1,hspace=hspace)
     
    fig.savefig(f'snapshot_mean_{comp}_{nepochs}.pdf')
    

'''
