


import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator
import h5py
import matplotlib.pyplot as plt
import os

import matplotlib as mpl

names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [0,0.5,0.5,0.5,0.4]
wspaces = [0.30,0.5,0.5,0.5,0.3]
sizes = [(14,6),(10,8),(10,6.0),(10,8),(14,7)]
shrinks = [0.5,0.7,0.7,0.75,0.75]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['hspace'] = hspaces[i]
        model[name]['wspace'] = wspaces[i]
        model[name]['size'] = sizes[i]
        model[name]['shrink'] = shrinks[i]
    return model

def read_h5file(file_name,key):
       # self.Multiple_Okada_displacement()

       f = h5py.File(file_name,'r')
       dset = np.array(f[key])
       return dset
    
def cov(file_name,key):
    dset  = read_h5file(file_name,key)
    covariance = np.cov(dset.transpose())
    
    nparameters = dset.shape[1]
    
    cov1 = covariance[:nparameters//3,:nparameters//3]
    cov2 = covariance[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    cov3 = covariance[2*nparameters//3:,2*nparameters//3:]

    # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
    # standard deviation (= square root of variance)

    std1 = np.sqrt(cov1.diagonal())
    std2 = np.sqrt(cov2.diagonal())
    std3 = np.sqrt(cov3.diagonal())
    
    return std1,std2,std3

def average(file_name,key):
    dset  = read_h5file(file_name,key)
    mean = np.mean(dset.transpose(),axis=1)
    
    nparameters = dset.shape[1]
    
    x = mean[:nparameters//3]
    y = mean[nparameters//3:2*nparameters//3]
    z = mean[2*nparameters//3:]

    # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
    # standard deviation (= square root of variance)

    
    return x,y,z

def corr(file_name,key):
    dset  = read_h5file()
    correlation = np.corrcoef(dset.transpose())
    
    nparameters = dset.shape[1]
    corr1 = correlation[:nparameters//3,:nparameters//3]
    corr2 = correlation[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    corr3 = correlation[2*nparameters//3:,2*nparameters//3:]

    return corr1, corr2, corr3

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
  
    
models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks)

markers = ['s','o','^','D','v']
msizes = [20,15,12,5,5]
colors = ['blue','','']

n_lines = 5
cmap = mpl.colormaps['hsv']

# Take colors at regular intervals spanning the colormap.
colors = cmap(np.linspace(0, 1, n_lines))
colors = ['red','yellow','skyblue','blue','grey']
names = ['Tohoku','Illapel','Iquique','Pedernales','Gorkha']
for i,name in enumerate(names):
    df = pd.read_csv(f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')
    nsamples = 100
    
    working_dir = os.getcwd()
    #edks_file = os.path.join(working_dir,f'{name}/{nsamples}_samples/EDKS/EDKS_{name}_displacement_nsamples_{nsamples}.h5')
    #okada_file = os.path.join(working_dir,f'{name}/{nsamples}_samples/Okada/Okada_{name}_displacement_nsamples_{nsamples}.h5')
    
    edks_file =f'FK_displacement_100_samples/EDKS_{name}_displacement_nsamples_100.h5'
    key = 'displacement'
    
    
    std1,std2,std3 = cov(edks_file,key)
    
    
    std = {'Along-strike':std1,
                'Trench-normal':std2,
                'Vertical':std3}
        
    patch = models[name]['patch']
    hspace = models[name]['hspace']
    
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(2,4,patch,nrows,ncols,control=1)
    
    Xr,Yr = np.meshgrid(x,y)
    
        
    
    
    
    
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
    
    xS = np.arange(patch/2 , ncols*patch,patch)
    yS = np.arange(-(nrows-1/2)*patch,0,patch)
    X,Y = np.meshgrid(xS,yS)
       # shift accordingly at surface
    #yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
    
    
    DIP = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0)
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0) 
    
    
    method = 'linear'
    
    DEPTH_interp = RegularGridInterpolator((yS, xS), DEPTH,bounds_error=False, fill_value=None,method = method)
    
    hypo_dep = DEPTH_interp((hypo_dd,hypo_as))
    
    r_target = np.array([hypo_as,hypo_dd,hypo_dep])
    
    
    
    
    dx,dy,dz = average(edks_file,key)
    
    d = np.sqrt(dx**2 + dy**2 + dz**2)
    xflat = Xr.flatten()
    yflat = Yr.flatten()
    
    
    r_all = np.column_stack((xflat,yflat))
    
    dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2 + (r_target[2])**2)
    
    dr_sorted = np.sort(dr)
    ind_dr_sorted = np.argsort(dr)


    plt.scatter(dr_sorted,d[ind_dr_sorted]*100,marker=markers[i],ec='black',lw=0.2,s=msizes[i],label=f'{name}',color=colors[i])
    Mw =[6,7,8,9]
    for m in Mw:
        
        x = np.arange(1,2000,100)
        y = 10**(-4.434 + 1.047*m - 0.138*m*np.log10(x))
        plt.plot(x,y,c='k',lw=0.25)
        plt.text(10,10**(-4.434 + 1.047*m - 0.138*m*np.log10(10)),f'Mw{m}',fontsize=7)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1,1e3)
    plt.ylim(1e-1,2e4)
    plt.xlabel('Hypocentral distance (km)')
    plt.ylabel('Permanent Deformation (cm)')
    plt.legend(handletextpad=0.02,borderaxespad=0.1,fontsize=8)
plt.savefig('Permament_EDKS_deformation_scaling.png', dpi=700)

