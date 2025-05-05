# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 15:58:33 2024

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 14:57:55 2024

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
scales = [0.5,0.8e-1,1e-1,1.25e-1,0.70e-1]
shrinks = [0.2,0.4,0.4,0.3,0.4]
hspaces = [0.01,0.3,0.3,0.3,0.1]
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


markers = ['s','o','^','D','v']
msizes = [8,6,6,4,4]
colors = ['blue','','']

n_lines = 5
cmap = mpl.colormaps['hsv']

# Take colors at regular intervals spanning the colormap.
colors = cmap(np.linspace(0, 1, n_lines))

colors = ['red','yellow','skyblue','blue','grey']
names = ['Tohoku','Illapel','Iquique','Pedernales','Gorkha']
fig,axes = plt.subplots(2,1,figsize=(7,6))
for i,name in enumerate(names):
    df = pd.read_csv(f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')
    nsamples = 100
    
    working_dir = os.getcwd()
    #edks_file = os.path.join(working_dir,f'{name}/{nsamples}_samples/EDKS/EDKS_{name}_displacement_nsamples_{nsamples}.h5')
    #okada_file = os.path.join(working_dir,f'{name}/{nsamples}_samples/Okada/Okada_{name}_displacement_nsamples_{nsamples}.h5')

    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    strike = models[name]['strike']
    magnitude = models[name]['magnitude']
    scale =  models[name]['scale']
    shrink = models[name]['shrink']
    hspace = models[name]['hspace']
    folder = os.path.join(working_dir,'H5_60')
    folder = os.path.join(folder,name)
    nmodels  = 60
    npoints = 256
    x = np.arange(patch/2,ncols*patch,patch)
    y = np.arange(patch/2,nrows*patch,patch)
    y = -np.flip(y,axis=0)
    
    Xr,Yr = np.meshgrid(x,y)
    
    strike = strike*np.pi/180
    strike_corr = strike*(180/np.pi) - 90
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
    
    GD = np.sqrt(N**2 + E**2 + Z**2)
    
    std_GD = np.zeros_like(mean_X)

    for t in range(npoints):
        std_GD[:,t] = np.sqrt(np.diag(np.cov(GD[:,:,t].T)))

    peak_ground = np.zeros(nstations)
    std_peak_ground =np.zeros(nstations)
    for l in range(nstations):
        dx,dy,dz = mean_X[l,:],mean_Y[l,:],mean_Z[l,:]
        gds = np.sqrt(dx**2+dy**2+dz**2) 
        peak_ground[l] = np.max(gds)
        tindex = np.argmin(abs(gds - np.max(gds)))
        std_peak_ground[l] = std_GD[l,tindex] 


    
    
        
    
    
    
    
    # if name=='Pedernales':
        
    #     RotAngle = 360-99
    #     RotAngle2 = RotAngle*((np.pi) / 180.)
         
    #     rotation = np.arctan2(np.tan(strike) - np.tan(RotAngle2), 
    #                          np.cos(dip)*(1.+np.tan(RotAngle2)*np.tan(strike)))
        
        
    #      # If RotAngle within ]90, 270], change rotation
    #     if RotAngle > 90. and RotAngle<=270:
    #          rotation += np.pi
        
        
    #     hypo_dd= hypo_as *np.cos(rotation) - hypo_dd*np.sin(rotation)
    #     hypo_as= hypo_as *np.sin(rotation) + hypo_dd*np.cos(rotation)
    
    hypo_dd = -df['Hypo_dd'].values[0]
    hypo_as = df['Hypo_as'].values[0]
    
    xr = np.arange(patch/2 , ncols*patch,patch)
    yr = np.arange(-(nrows-1/2)*patch,0,patch)
    X,Y = np.meshgrid(xr,yr)
       # shift accordingly at surface
    #yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
    
    
    DIP = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0)
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0) 
    
    
    method = 'linear'
    
    DEPTH_interp = RegularGridInterpolator((yr, xr), DEPTH,bounds_error=False, fill_value=None,method = method)
    
    hypo_dep = DEPTH_interp((hypo_dd,hypo_as))
    
    r_target = np.array([hypo_as,hypo_dd,hypo_dep])
    
    
    
    
    xflat = Xr.flatten()
    yflat = Yr.flatten()
    
    
    xflat = Xr.flatten()
    yflat = Yr.flatten()
    
    
    r_all = np.column_stack((xflat,yflat))
    
    drx = r_all[:,0]  - r_target[0] 
    dry = r_all[:,1]  - r_target[1]
    azimuths = np.arctan2(drx,dry)*(180/np.pi)
    
    # for n,az in enumerate(azimuths):
    #     if az<0:
    #         az+=360
    #         azimuths[n] = az
    # azimuths  += strike_corr
    
    for l,az in enumerate(azimuths):
        if az>360:
            az-=360
            azimuths[l] = az
      
    for n,az in enumerate(azimuths):
        if az<0:
            az+=360
            azimuths[n] = az    

    azimuths -= 90
    azimuth_sorted = np.sort(azimuths)
    ind_azimuth_sorted = np.argsort(azimuths)

    r_all = np.column_stack((xflat,yflat))
    
    dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2 + (r_target[2])**2)
    
    dr_sorted = np.sort(dr)
    ind_dr_sorted = np.argsort(dr)
    
    #plt.scatter(dr_sorted,peak_ground[ind_dr_sorted]*100, marker=markers[i],ec='black',lw=0.2,s=msizes[i],label=f'{name}',color=colors[i])
    #plt.errorbar(dr_sorted,peak_ground[ind_dr_sorted]*100, yerr =std_peak_ground[ind_dr_sorted]*100,mec='black',mew=0.5,lw=0.01,ms=msizes[i],fmt='o',label=f'{name}',color=colors[i])
    axes[0].errorbar(dr_sorted,peak_ground[ind_dr_sorted]*100, yerr =std_peak_ground[ind_dr_sorted]*100,mec='k',mew=0.2,elinewidth =0.4,capsize=0.75,fmt='.',label='$M_w$'+f'{magnitude} {name}',ms=msizes[i],color=colors[i])
    Mw =[6,7,8,9]
    r = [30,125,320,700]
    for l,m in enumerate(Mw):
        
        x = np.arange(1,2000,100)
        y = 10**(-4.434 + 1.047*m - 0.138*m*np.log10(x))
        axes[0].plot(x,y,c='k',lw=0.25)
        axes[0].text(r[l],10**(-4.434 + 1.047*m - 0.138*m*np.log10(r[l])),f'Mw{m}',fontsize=6)
        
    axes[0].set_xscale('log')
    axes[0].set_yscale('log')
    axes[0].set_ylim(3,1e4)
    axes[0].set_xlim(3,1e3)

    axes[0].set_xlabel('Hypocentral distance (km)')
    axes[0].set_ylabel('PGD (cm)')
    axes[0].legend(loc='lower left',handletextpad=0.02,borderaxespad=0.1,fontsize=8)
    
    axes[0].text(-0.1, 1.2, 'a',transform=axes[0].transAxes,fontsize = 18)

    axes[1].errorbar(azimuth_sorted,peak_ground[ind_azimuth_sorted]*100, yerr =std_peak_ground[ind_azimuth_sorted]*100,mec='k',mew=0.2,elinewidth =0.4,capsize=0.75,fmt='.',ms=msizes[i],color=colors[i])
    
    #plt.scatter(dr_sorted,peak_ground[ind_dr_sorted]*100, marker=markers[i],ec='black',lw=0.2,s=msizes[i],label=f'{name}',color=colors[i])
    #plt.errorbar(dr_sorted,peak_ground[ind_dr_sorted]*100, yerr =std_peak_ground[ind_dr_sorted]*100,mec='black',mew=0.5,lw=0.01,ms=msizes[i],fmt='o',label=f'{name}',color=colors[i])
    #plt.scatter(azimuth_sorted,peak_ground[ind_azimuth_sorted]*100)
    #plt.errorbar(azimuth_sorted,peak_ground[ind_azimuth_sorted]*100, yerr =std_peak_ground[ind_azimuth_sorted]*100,mec='k',mew=0.2,elinewidth =0.4,capsize=0.75,fmt='.',label='$M_w$'+f'{magnitude} {name}',ms=msizes[i],color=colors[i])

    axes[1].text(-0.1, 1.2, 'b',transform=axes[1].transAxes,fontsize = 18) 
    axes[1].set_yscale('log')
    axes[1].set_ylim(3,1e4)
    axes[1].set_xticks(list(np.arange(-90,271,45)))
    axes[1].grid(axis='x',linestyle='--',linewidth=0.4)
    axes[1].set_xlabel(r'Angle measured from strike line')
    axes[1].set_ylabel('PGD (cm)')
    fig.subplots_adjust(hspace=0.4)    
plt.savefig('PGD_wrt_distance_and_strike.pdf')
plt.close()

#### Azimuth dependence ###


