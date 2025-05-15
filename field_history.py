# -*- coding: utf-8 -*-
"""
Created on Tue Sep 24 12:03:44 2024

@author: vite0005
"""
import re
import obspy as obs
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
import sys
from matplotlib.patheffects import Stroke
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapely.geometry as sgeom
import numpy as np
import cartopy.io.img_tiles as cimgt

from scipy.interpolate import RegularGridInterpolator
from matplotlib import pyplot as plt, animation

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']

nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patches = [29,17,18,10,15]
geoms = [(nrows[i],ncols[i]) for i in range(len(names))]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,90,90,107,360-99]
ramps = [0,3,0,0,9]
factors =[3,3,3,2,2]
strikes = [194,-13.58,4,293,27.05]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,nramps,factors,strikes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i]
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i]
        model[name]['nramp'] = nramps[i]
        model[name]['factor'] = factors[i]
        model[name]['factor'] = factors[i]
        model[name]['strike'] = strikes[i]
    return model


models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,ramps,factors,strikes)
eq_name = 'Tohoku'


patch = models[eq_name]['patch']
nrows,ncols = models[eq_name]['geom']
strike  = models[eq_name]['strike']

x = np.arange(patch/2,ncols*patch,patch)
y = np.arange(patch/2,nrows*patch,patch)
y = -np.flip(y,axis=0)


strike = strike*np.pi/180
strike = strike - 90.0*np.pi/180.0

def set_stn(c,factor,patch,geom,control=0):
      nrows,ncols =geom[0],geom[1]
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn

def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)

def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()

def proj_ysrc_coords(patch,dip):
      proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
      proj_ysrc = np.zeros_like(proj_dysrc)
      for i in range(len(proj_ysrc)):
          proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i]
      ysrc = np.flip(proj_ysrc)

      return ysrc


nmodels = 10
nsamples = 100

nstations = nrows*ncols
npoints = 256
shape = (nmodels,nstations,npoints,3)
data = np.zeros(shape)

parent_dir = '/home/josevilo/Dynamic/MudPy/examples/%s/%s_samples'%(eq_name,nsamples)
for i in range(1,nmodels+1):
   imodel_dir = parent_dir+f'/{i}/output/waveforms/{eq_name}'
   for m,comp in enumerate(['N','E','Z']):
       for k in range(nstations):
           file_dir =imodel_dir+'/%.4d.LY%s.sac'%(k,comp)
           stream = obs.read(file_dir)
     
           trace = stream[0].data

           data[i-1,k,:,m] = stream[0].data 
           


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [0,0.5,0.5,0.5,0.5]
wspaces = [0.2,0.3,0.4,0.3,0.4]
sizes = [(14,16),(10,11),(10,12),(10,11),(14,16)]
shrinks = [0.5,0.5,0.5,0.5,0.2]
limits_model = [[138,145,35,44],[-73,-69,-21,-18],[-74,-70,-32,-29],[-82,-79,-1,1],[84,87,27,29]]
limits_disp = [[138-1,145+2,35-2,44+2],[-72-1,-69+1,-21-1,-18+1],[-72-2,-70+1,-33-1,-29+1],[-81-1,-80+1,-1-1,1+1],[84-0.25,87,27-0.25,29+0.25]]
nsamples = 100
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,limits_model,limits_disp):
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
        model[name]['limits_model'] = limits_model[i]
        model[name]['limits_disp'] = limits_disp[i]
    return model

models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,limits_model,limits_disp)

N = data[:,:,:,0]  
E = data[:,:,:,1]  
name = 'Tohoku'

df = pd.read_csv(f'/home/josevilo/Dynamic/MudPy/examples/{name}/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')
X = df['lon'].values
Y = df['lat'].values

print(X)

# geographical projection used: 'PlateCarree'
proj = ccrs.PlateCarree()

    



U = data[0,:,0,0]
V = data[0,:,0,1]
Z = data[0,:,0,2]

limits = models[name]['limits_disp']

fig = plt.figure(dpi=600)
axtext = fig.add_axes([0.0,0.95,0.1,0.05])
axtext.axis("off")

ax = plt.axes(projection=proj)

#google_terrain = cimgt.GoogleTiles(style="satellite")
#ax.add_image(google_terrain, 10)

min_lon,max_lon,min_lat,max_lat = limits

  # set grid geographical boundaries for aoi ie. [min_lon, max_Lon,min_lat,max_lat]
ax.set_extent([min_lon, max_lon, min_lat, max_lat])

#ax.coastlines(linewidth =0.4)
#ax.add_feature(cfeature.BORDERS,linewidth=0.8)
LON = one_array_formatter(X,shape)
LAT = one_array_formatter(Y,shape)
scale = 30
q = ax.quiver(LON,LAT,V,U,scale=scale,scale_units ='x', color='black',units='width',width=0.002)


sep_map = 0.25
cmap = plt.get_cmap('bwr') # Choose your desired colormap
norm = TwoSlopeNorm(0,vmin=(min(Z.min(),-Z.max())),vmax=(max(-Z.min(),Z.max())))
cax1 = ax.tricontourf(LON,LAT,Z,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,levels = np.arange(round(min(Z.min(),-Z.max()),1),round(max(-Z.min(),Z.max()),1) + sep_map,sep_map))
ax.set_title(f'Dynamic displacement')

time = axtext.text(0.5,-0.25, '$t=$'+str(0)+' s', ha="left", va="top")
bbox = ax.get_position()
Ukey = 20

ax.quiverkey(q, X=bbox.x0 + 0.88*bbox.x1, Y=bbox.y0 - 0.8*bbox.y0  , U=Ukey,
                label=' {}m'.format(Ukey),labelpos='N',angle=0,labelcolor='white',fc='white',labelsep = 0.015,fontproperties={'size': 8,'weight':'bold'})
factor = 1/q.scale    # scaling factor from UV unit to XY units
fig.savefig('test')
'''
  # Set the colorbar limit
  # plot GPS stations using corresponding lon and lat coordinates
  #ax.scatter(X,Y,s=3,ec='k', c='k',linewidths=0.25,marker='.',
  #         transform=ccrs.PlateCarree(),label='GPS station')
  #ax.contourf(LONGITUDE,LATITUDE,E,transform=ccrs.Geodetic(),cmap=cmap)
  # set gridlines and outer geographical labels
gl = ax.gridlines(crs = ccrs.PlateCarree(),draw_labels=True, linewidth=0.1,linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xlabel_style = {'size': 7}
gl.ylabel_style = {'size': 7}

x_coord,y_coord = [0.3,0.7]
sub_ax = fig.add_axes([x_coord, y_coord, 0.10, 0.10],
                            projection=ccrs.PlateCarree())
extent = [min_lon -15 ,max_lon + 12,min_lat - 20 ,max_lat + 5]

sub_ax.set_extent(extent, geodetic)

  # Make a nice border around the inset axes.
#effect = Stroke(linewidth=1.5)
#sub_ax.spines['geo'].set_path_effects([effect])

  # Add the land, coastlines and the extent of the Solomon Islands.
sub_ax.add_feature(cfeature.LAND,color='lightgray')
sub_ax.add_feature(cfeature.BORDERS,linewidth=0.4)
sub_ax.coastlines(linewidth =0.4)
extent_box = sgeom.box(min_lon, min_lat, max_lon, max_lat)
sub_ax.add_geometries([extent_box], ccrs.PlateCarree(), facecolor='none',
                        edgecolor='blue', linewidth=1.25)

nframes = 10
def animate(i):
    q.set_UVC(data[0,:,i,1],data[0,:,i,0])
    cax1.set_array(data[0,:,i,2])

    # arr = G0[:,:,i]
    # vmax     = np.max(arr)
    # vmin     = np.min(arr)

    # cf = ax.pcolormesh(x,y,arr, vmax=vmax, vmin=vmin)
    # fig.colorbar(cf, ax=ax)
    time.set_text('$t=$'+str(i)+' s')
plt.subplots_adjust(hspace=-0.2)
anim = animation.FuncAnimation(fig, animate,frames =nframes,interval =2)
anim.save(filename=f"{eq_name}_disp_field_{nmodels}.gif", writer="pillow")
'''
