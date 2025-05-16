# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 18:56:39 2024

@author: joanv
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.colors import TwoSlopeNorm
import pandas as pd

names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
#hspaces = [0,0.5,0.5,0.5,0.4]
hspaces = [-0.25,-0.2,-0.2,-0.1,0.05]
wspaces = [0.30,0.4,0.4,0.5,0.3]
sizes = [(14,6),(10,8),(10,6.0),(10,8),(14,7)]
shrinks = [0.5,0.6,0.6,0.6,0.75]
trench_name = ['kuriljapan','southamerica','southamerica','southamerica','MHT']
lon_lim = [(138,145),(-74,-68),(-75,-69),(-81,-80),(84,87)]
lat_lim = [(34,44),(-22,-17),(-34,-28),(-1.,1.),(26,30)]


def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,trench_name,lon_lim,lat_lim):
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
        model[name]['trench'] = trench_name[i]
        model[name]['lon_lim'] = lon_lim[i]
        model[name]['lat_lim'] = lat_lim[i]
    return model

models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,trench_name,lon_lim,lat_lim)

nsamples = 100
name = str(sys.argv[1])
working_dir = os.getcwd()
edks_file = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS/{nsamples}_samples/EDKS_{name}_displacement_nsamples_{nsamples}_parallel.h5')
okada_file = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/Okada/{nsamples}_samples/Okada_{name}_displacement_nsamples_{nsamples}_parallel.h5')
key = 'displacement'

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
  
def plot_cov(name,file_name,key,method = 'EDKS'):
    std1,std2,std3 = cov(file_name,key)


    std = {'Along-strike':std1,
            'Trench-normal':std2,
            'Vertical':std3}
    
    patch = models[name]['patch']
    hspace = models[name]['hspace']

    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    
    x,y = set_stn(2,4,patch,nrows,ncols,control=1)
    nrows = len(y)
    ncols = len(x)
    X,Y = np.meshgrid(x,y)



    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=1400)
    for i,parameter_id in enumerate(std.keys()):
        parameter = std[parameter_id]

        title = f"{parameter_id} Displacement "

        im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow',linewidth=0.25)


        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)

        axes[i].margins(0)        
        fig.colorbar(im0, ax=axes[i],shrink=0.65,label='Uncertainty (m)')
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=12)
    
    plt.subplots_adjust(hspace=hspace)
    fig.suptitle(f'{name } uncertainty in {method} Surface Displacement')
    os.makedirs(os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/{method}/cov/{nsamples}_samples'),exist_ok=True)
    figure_dir = os.path.join(os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/{method}/cov/{nsamples}_samples'),f'{method}_{name}_uncertainty_parallel.png')	
    fig.savefig(figure_dir,dpi=700)
        
	
plot_cov(name,edks_file,key)
plot_cov(name,okada_file,key,method='Okada')


def plot_mean_cov(name,file_name,key,method = 'EDKS'):
    std1,std2,std3 = cov(file_name,key)
    d1,d2,d3 = average(file_name,key)
    
    d = {'Along-strike':d1,
            'Trench-normal':d2,
            'Vertical':d3}
    

    std = {'Along-strike':std1,
            'Trench-normal':std2,
            'Vertical':std3}
    
    
    rel_std = {'Along-strike':100*std1/np.abs(d1),
            'Trench-normal':100*std2/np.abs(d2),
            'Vertical':100*std3/np.abs(d3)}
    
    patch = models[name]['patch']
    hspace = models[name]['hspace']
    wspace = models[name]['wspace']
    size = models[name]['size']
    shrink = models[name]['shrink']
    trench_name = models[name]['trench']
    lat_lim = models[name]['lat_lim']
    lon_lim = models[name]['lon_lim']
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    mean_data_dir = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')

    trench_dir = os.path.join(working_dir,f'Trench_TPGA/{trench_name}.lonlat')
    df = pd.read_csv(mean_data_dir)

    LON = np.flip(df['lon'].values.reshape(nrows,ncols,order='F'),axis=0)
    LAT = np.flip(df['lat'].values.reshape(nrows,ncols,order='F'),axis=0)
    lon = LON.flatten()
    lat = LAT.flatten()
    trench = np.loadtxt(trench_dir)
    trench_lon = trench[:,0]
    trench_lat = trench[:,1]

    ind_lat = np.where((trench_lat<lat_lim[1]) & (trench_lat>lat_lim[0]))
    trench_lat = trench_lat[ind_lat]
    trench_lon = trench_lon[ind_lat]
    
    yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
    
    xs,ys = set_stn(1,1,patch,nrows,ncols,control=0)
    X,Y = np.meshgrid(xs,yS)

    X[:,0] = X[:,0] - patch/2
    X[:,-1] = X[:,-1] + patch/2
    Y[0,:] = Y[0,:] - patch/2
    Y[-1,:] = Y[-1,:] + patch/2

     
    
    coast_lon = []
    coast_lat = []
    import cartopy.feature as cfeature
    for geom in cfeature.COASTLINE.geometries():
        lons = geom.coords._coords[:,0]
        lats = geom.coords._coords[:,1]
        for lon0,lat0 in zip(lons,lats):
            if (lat0<lat_lim[1]) & (lat0>lat_lim[0]) & (lon0<lon_lim[1]) & (lon0>lon_lim[0]):
                coast_lon.append(lon0)
                coast_lat.append(lat0)

    if name=='Gorkha':
        
        import cartopy.io.shapereader as shpreader
        shpfilename = shpreader.natural_earth(resolution='110m',
                                          category='cultural',
                                          name='admin_0_countries')
        reader = shpreader.Reader(shpfilename)
        countries = reader.records()
        for country in countries:
            if country.attributes['ADM0_A3'] == 'NPL':
                l = country.geometry
                
        coast_lon = np.array(l.boundary.xy[0])  
        coast_lat = np.array(l.boundary.xy[1])    
    # ind_coast_lat = np.where((coast_lat<lat_lim[1]) & (coast_lat>lat_lim[0]))
    # coast_lat = coast_lat[ind_coast_lat]
    # coast_lon = coast_lon[ind_coast_lat]

    # ind_coast_lon = np.where((coast_lon<lon_lim[1]) & (coast_lon>lon_lim[0]))   
    # coast_lat = coast_lat[ind_coast_lon]
    # coast_lon = coast_lon[ind_coast_lon]
    methodInt = 'linear'

    x,y = set_stn(2,4,patch,nrows,ncols,control=1)
    xl,yl = set_stn(40,4,patch,nrows,ncols,control=1)
    Xl,Yl = np.meshgrid(xl,yl)
    
    
    from scipy.interpolate import RegularGridInterpolator,NearestNDInterpolator
    
    lat_interp = RegularGridInterpolator((ys, xs), LAT,bounds_error=False, fill_value=None,method = methodInt)
    lon_interp = RegularGridInterpolator((ys, xs), LON,bounds_error=False, fill_value=None,method = methodInt)
    
    LON2 = lon_interp((Yl,Xl))
    LAT2 = lat_interp((Yl,Xl))
    lon2 = LON2.flatten()
    lat2 = LAT2.flatten()
    interp_x = NearestNDInterpolator(list(zip(lon2,lat2)), Xl.flatten())
    interp_y = NearestNDInterpolator(list(zip(lon2,lat2)), Yl.flatten())

    trench_x = interp_x(trench_lon,trench_lat)
    trench_y = interp_y(trench_lon,trench_lat)
    if name=='Tohoku':
        
        coast1 = np.loadtxt(os.path.join(os.getcwd(),f'boundaries/{name}_coast1.lonlat'))
        coast_lon1  = coast1[:,0]
        coast_lat1  = coast1[:,1]
        coast_x1 = interp_x(coast_lon1,coast_lat1)
        coast_y1 = interp_y(coast_lon1,coast_lat1)
        
        coast2 = np.loadtxt(os.path.join(os.getcwd(),f'boundaries/{name}_coast2.lonlat'))
        coast_lon2  = coast2[:,0]
        coast_lat2  = coast2[:,1]
        coast_x2 = interp_x(coast_lon2,coast_lat2)
        coast_y2 = interp_y(coast_lon2,coast_lat2)
    
    else:
        coast = np.loadtxt(os.path.join(os.getcwd(),f'boundaries/{name}_coast.lonlat'))
        coast_lon  = coast[:,0]
        coast_lat  = coast[:,1]
        coast_x = interp_x(coast_lon,coast_lat)
        coast_y = interp_y(coast_lon,coast_lat)
        
    nrows = len(y)
    ncols = len(x)


    
    # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    fig, axes = plt.subplots(3,3,figsize=size,sharex=True,dpi=700)
    symbol_d = ['$d_x$ (m)','$d_y$ (m)','$d_z$ (m)']
    symbol_sigma = ['$s$ (m)','$s$ (m)','$s$ (m)']
    for i,parameter_id in enumerate(std.keys()):
        sigma = std[parameter_id]
        val = d[parameter_id]
        rel_sigma = rel_std[parameter_id]
        title = f"{parameter_id}"

        im0 = axes[i][0].pcolormesh(x,y,val.reshape(nrows,ncols),cmap='bwr',norm=TwoSlopeNorm(0,vmin=(min(val.min(),-val.max())),vmax=(max(-val.min(),val.max()))))
        im1 = axes[i][1].pcolormesh(x,y,sigma.reshape(nrows,ncols),cmap='rainbow')
        im2 = axes[i][2].pcolormesh(x,y,rel_sigma.reshape(nrows,ncols),cmap='Blues',vmin=0,vmax=100)
        axes[i][0].plot(X[0,:],Y[0,:],color='k',ls='--',lw=0.5)
        axes[i][0].plot(X[-1,:],Y[-1,:],color='k',ls='--',lw=0.5)
        axes[i][0].plot(X[:,0],Y[:,0],color='k',ls='--',lw=0.5)
        axes[i][0].plot(X[:,-1],Y[:,-1],color='k',ls='--',lw=0.5)
        axes[i][1].plot(X[0,:],Y[0,:],color='k',ls='--',lw=0.5)
        axes[i][1].plot(X[-1,:],Y[-1,:],color='k',ls='--',lw=0.5)
        axes[i][1].plot(X[:,0],Y[:,0],color='k',ls='--',lw=0.5)
        axes[i][1].plot(X[:,-1],Y[:,-1],color='k',ls='--',lw=0.5)
        axes[i][2].plot(X[0,:],Y[0,:],color='k',ls='--',lw=0.5)
        axes[i][2].plot(X[-1,:],Y[-1,:],color='k',ls='--',lw=0.5)
        axes[i][2].plot(X[:,0],Y[:,0],color='k',ls='--',lw=0.5)
        axes[i][2].plot(X[:,-1],Y[:,-1],color='k',ls='--',lw=0.5)

        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        axes[i][0].plot(trench_x,trench_y,color='black',lw=0.8)
        axes[i][1].plot(trench_x,trench_y,color='black',lw=0.8)

        axes[i][2].plot(trench_x,trench_y,color='black',lw=0.8)

        if name=='Tohoku':
            axes[i][0].plot(coast_x1,coast_y1,color='gray',lw=0.8)
            axes[i][1].plot(coast_x1,coast_y1,color='gray',lw=0.8)
            axes[i][2].plot(coast_x1,coast_y1,color='gray',lw=0.8)
            
            axes[i][0].plot(coast_x2,coast_y2,color='gray',lw=0.8)
            axes[i][1].plot(coast_x2,coast_y2,color='gray',lw=0.8)
            axes[i][2].plot(coast_x2,coast_y2,color='gray',lw=0.8)
        else:
    
            axes[i][0].plot(coast_x,coast_y,color='gray',lw=0.8)
            axes[i][1].plot(coast_x,coast_y,color='gray',lw=0.8)
            axes[i][2].plot(coast_x,coast_y,color='gray',lw=0.8)

       # axes[i][0].margins(0)
        #axes[i][1].margins(0)
        #axes[i][2].margins(0)
        fig.colorbar(im0, ax=axes[i][0],shrink=shrink,label= symbol_d[i])
        fig.colorbar(im1, ax=axes[i][1],shrink=shrink,label= symbol_sigma[i])
        fig.colorbar(im2, ax=axes[i][2],shrink=shrink,label='rel. error (%)',extend='max')
        

        # modified [1,0]->[i,0] size 9 -> 7.5    
        axes[1][0].set_ylabel('Trench-normal distance (km)',fontsize=7)
        axes[2][i].set_xlabel('Along-strike distance (km)',fontsize=7)
        axes[i][0].set_aspect('equal', 'box')
        axes[i][0].set_title(title,fontweight='bold')
        axes[i][0].tick_params(labelsize=7)

        axes[0][i].text(-0.1,1.25,'abc'[i],fontweight='bold',fontsize=20,transform=axes[0][i].transAxes)
        
        # axes[i][1].set_ylabel('Trench-normal distance (km)',fontsize=10)
        # axes[i][1].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i][1].set_aspect('equal', 'box')
        axes[i][1].set_title(title,fontweight='bold')
        axes[i][1].tick_params(labelsize=7)
        
        # axes[i][2].set_ylabel('Trench-normal distance (km)',fontsize=10)
        # axes[i][2].set_xlabel('Along-strike distance (km)',fontsize=10)
        axes[i][2].set_aspect('equal', 'box')
        axes[i][2].set_title(title,fontweight='bold')
        axes[i][2].tick_params(labelsize=7)
    
    plt.subplots_adjust(wspace=wspace,hspace=hspace)
    #fig.suptitle(f'{name}',fontsize=13,fontweight='bold')
    os.makedirs(os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/{method}/{nsamples}_samples'),exist_ok=True)
    figure_dir = os.path.join(os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/{method}/{nsamples}_samples'),f'{method}_{name}_mean_errors_parallel.pdf')
    fig.savefig(figure_dir)

plot_mean_cov(name,edks_file,key)
plot_mean_cov(name,okada_file,key,method='Okada')
# def cov(self):
#     dset  = self.read_h5file()
#     self.covariance = np.cov(dset.transpose())
    
#     nparameters = dset.shape[1]
    
#     self.cov1 = self.covariance[:nparameters//3,:nparameters//3]
#     self.cov2 = self.covariance[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
#     self.cov3 = self.covariance[2*nparameters//3:,2*nparameters//3:]

#     # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
#     # standard deviation (= square root of variance)

#     self.std1 = np.sqrt(self.cov1.diagonal())
#     self.std2 = np.sqrt(self.cov2.diagonal())
#     self.std3 = np.sqrt(self.cov3.diagonal())
