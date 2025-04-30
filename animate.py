# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 20:44:37 2024

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
from matplotlib.patheffects import Stroke
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import shapely.geometry as sgeom
import numpy as np
import cartopy.io.img_tiles as cimgt
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from scipy.interpolate import RegularGridInterpolator
from matplotlib import pyplot as plt, animation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']




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
limits_disp = [[139.5-1,144.5+1,35-2,43+1],[-72-1,-69+1,-21-1,-18+1],[-72-2,-70+1,-33-1,-29+1],[-81-1,-80+1,-1-1,1+1],[84-0.25,87,27-0.25,29+0.25]]
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



nsamples = 100
npoints = 256
def animate_displacement(name,sampling,Ukey_horizontal,scale_horizontal,Ukey_vertical,scale_vertical,insetbox_x,insetbox_y,textbox_y,loc='upper right'):
  eq_name = name


  patch = models[eq_name]['patch']
  nrows,ncols = models[eq_name]['geom']


  x = np.arange(patch/2,ncols*patch,patch)
  y = np.arange(patch/2,nrows*patch,patch)
  y = -np.flip(y,axis=0)
  df = pd.read_csv(f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')
  U_parallel = df['U_parallel'].values
  U_perp = df['U_perp'].values
  Slip = np.sqrt(U_parallel**2 + U_perp**2)
  strike = np.mean(df['strike'].values)
  strike = strike*np.pi/180
  strike = strike - 90.0*np.pi/180.0
  nmodels = 1
  nstations = nrows*ncols
  shape = (nmodels,nstations,npoints,3)
  data = np.zeros(shape)

  parent_dir = os.path.join(os.getcwd(),f'Dynamic_Simulations/{name}/{nsamples}_samples/1/output/waveforms/{name}')
  for i in range(1,nmodels+1):
    imodel_dir = parent_dir
    for m,comp in enumerate(['N','E','Z']):
        for k in range(nstations):
            file_dir =imodel_dir+'/%.4d.LY%s.sac'%(k,comp)
            stream = obs.read(file_dir)

            data[i-1,k,:,m] = stream[0].data






  X = df['lon'].values
  Y = df['lat'].values


  geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
  # geographical projection used: 'PlateCarree'
  proj = ccrs.PlateCarree()




  m = sampling

  U = data[0,:,0,0]
  V = data[0,:,0,1]
  Z = data[0,:,0,2]

  Uts = data[0,:,::m,0]
  Vts = data[0,:,::m,1]
  Zts = data[0,:,::m,2]
  Zt = data[0,:,-1,2]
  limits = models[name]['limits_disp']


  fig,axes = plt.subplots(
    1, 2, figsize=(9, 6),dpi=300,sharey=True,
    subplot_kw={'projection': proj, "aspect": 2},
    gridspec_kw = {'wspace':0.14, 'hspace':0.007},
)

  axtext = fig.add_axes([0.08,textbox_y,0.1,0.05])
  axtext.axis("off")
  google_terrain = cimgt.GoogleTiles(style="satellite")
  axes[0].add_image(google_terrain, 10)

  min_lon,max_lon,min_lat,max_lat = limits

    # set grid geographical boundaries for aoi ie. [min_lon, max_Lon,min_lat,max_lat]
  axes[0].set_extent([min_lon, max_lon, min_lat, max_lat])

  axes[0].coastlines(linewidth =0.4)
  axes[0].add_feature(cfeature.BORDERS,linewidth=0.8)
  LON = one_array_formatter(X,(nrows,ncols))
  LAT = one_array_formatter(Y,(nrows,ncols))
  scale = scale_horizontal

  sep_map = 0.2
  #cmap = plt.get_cmap('bwr') # Choose your desired colormap
  #norm = TwoSlopeNorm(0,vmin=(min(Zt.min(),-Zt.max())),vmax=(max(-Zt.min(),Zt.max())))
  #cax1 = ax.tricontourf(LON,LAT,Z,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,levels = np.arange(round(min(Zt.min(),-Zt.max()),1),round(max(-Zt.min(),Zt.max()),1) + sep_map,sep_map))
  axes[0].set_title(' Horizontal Deformation Field',fontweight='bold',fontsize=10.5)


  #cbar = fig.colorbar(cax1,ax=ax,shrink = 0.4,cmap=cmap,pad = 0.075,norm=norm,orientation ='horizontal')
  #cax1 = ax.tricontourf(LON,LAT,Zt,transform=ccrs.PlateCarree(),cmap=cmap)
  #ax.set_title(f'Dynamic displacement')
  #cbar = fig.colorbar(cax1,ax=ax,shrink = 0.4,cmap=cmap,pad = 0.075,norm=norm,orientation ='horizontal')

  #cbar.set_label(label='$d_z$ (m)')
  #cbar.ax.tick_params(labelsize=7.5)
  time = axtext.text(0.5,-0.25, '$t=$'+str(0)+' s', ha="left", va="top")
  bbox = axes[0].get_position()
  Ukey = Ukey_horizontal
  im1 = axes[0].tricontourf(X,Y,Slip,alpha=1.0,cmap='YlOrRd')
  
  axins = inset_axes(axes[0], width="18%", height="1.5%", loc=loc)
# Add the colorbar
  ticks = [int(round(i,0)) for i in list(np.arange(0,int(max(Slip))  + round(int(max(Slip))/2,1) ,round(int(max(Slip))/2,1))) ]  
    
  cbar = plt.colorbar(im1, cax=axins,orientation='horizontal',ticks=ticks)
  cbar.ax.tick_params(labelsize=5,color='white')
  cbar.ax.set_xticklabels(ticks,color='white')  # horizontal colorbar
  cbar.ax.set_title('U(m)',color='white',pad = 0.5,fontsize=4,fontweight='bold')

  q = axes[0].quiver(LON,LAT,V,U,scale=scale,scale_units ='x', color='cyan',units='width',width=0.004)


  axes[0].quiverkey(q, X=bbox.x0 + 0.8*bbox.x1, Y=bbox.y0 - 0.2*bbox.y0  , U=Ukey_horizontal,
                  label=' {}m'.format(Ukey),labelpos='N',angle=0,labelcolor='cyan',fc='cyan',labelsep = 0.015,fontproperties={'size': 9.5,'weight':'bold'})
   # scaling factor from UV unit to XY units

    # Set the colorbar limit
    # plot GPS stations using corresponding lon and lat coordinates
    #ax.scatter(X,Y,s=3,ec='k', c='k',linewidths=0.25,marker='.',
    #         transform=ccrs.PlateCarree(),label='GPS station')
    #ax.contourf(LONGITUDE,LATITUDE,E,transform=ccrs.Geodetic(),cmap=cmap)
    # set gridlines and outer geographical labels

  cardinal_labels = {"east": "", "west": "", "north": "", "south": ""}
  latitude_formatter = LatitudeFormatter(cardinal_labels=cardinal_labels)
  longitude_formatter = LongitudeFormatter(cardinal_labels=cardinal_labels)
  gl = axes[0].gridlines(
      draw_labels=["left", "bottom"],
      linestyle="--",
      xformatter=longitude_formatter,
      yformatter=latitude_formatter,
      ylabel_style={"size": 8},
      xlabel_style={"size": 8},
      linewidth=0.3
  )

  x_coord,y_coord = [insetbox_x,insetbox_y]
  sub_ax = fig.add_axes([x_coord, y_coord, 0.15, 0.15],
                              projection=ccrs.PlateCarree())
  extent = [min_lon -15 ,max_lon + 15,min_lat - 25 ,max_lat + 15]

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


  google_terrain = cimgt.GoogleTiles(style="satellite")
  axes[1].add_image(google_terrain, 10)

  min_lon,max_lon,min_lat,max_lat = limits

    # set grid geographical boundaries for aoi ie. [min_lon, max_Lon,min_lat,max_lat]
  axes[1].set_extent([min_lon, max_lon, min_lat, max_lat])

  axes[1].coastlines(linewidth =0.4)
  axes[1].add_feature(cfeature.BORDERS,linewidth=0.8)
  LON = one_array_formatter(X,(nrows,ncols))
  LAT = one_array_formatter(Y,(nrows,ncols))
  scale = scale_vertical

  sep_map = 0.2
  #cmap = plt.get_cmap('bwr') # Choose your desired colormap
  #norm = TwoSlopeNorm(0,vmin=(min(Zt.min(),-Zt.max())),vmax=(max(-Zt.min(),Zt.max())))
  #cax1 = ax.tricontourf(LON,LAT,Z,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,levels = np.arange(round(min(Zt.min(),-Zt.max()),1),round(max(-Zt.min(),Zt.max()),1) + sep_map,sep_map))
  axes[1].set_title('Vertical Deformation Field',fontweight='bold',fontsize=10.5)


  #cbar = fig.colorbar(cax1,ax=ax,shrink = 0.4,cmap=cmap,pad = 0.075,norm=norm,orientation ='horizontal')
  #cax1 = ax.tricontourf(LON,LAT,Zt,transform=ccrs.PlateCarree(),cmap=cmap)
  #ax.set_title(f'Dynamic displacement')
  #cbar = fig.colorbar(cax1,ax=ax,shrink = 0.4,cmap=cmap,pad = 0.075,norm=norm,orientation ='horizontal')

  #cbar.set_label(label='$d_z$ (m)')
  #cbar.ax.tick_params(labelsize=7.5)

  bbox = axes[0].get_position()
  Ukey = Ukey_vertical
  im2 = axes[1].tricontourf(X,Y,Slip,alpha=1.0,cmap='YlOrRd')
  axins = inset_axes(axes[1], width="18%", height="1%", loc=loc)
# Add the colorbar
  cbar = plt.colorbar(im2, cax=axins,orientation='horizontal',ticks=ticks)
  cbar.ax.tick_params(labelsize=5,color='white')
  cbar.ax.set_xticklabels(ticks,color='white')  # horizontal colorbar
  cbar.ax.set_title('U(m)',color='white',pad = 0.5,fontsize=4,fontweight='bold')

  q2 = axes[1].quiver(LON,LAT,np.zeros_like(Z),Z,scale=scale,scale_units ='x', color='cyan',units='width',width=0.004)


  axes[1].quiverkey(q2, X=bbox.x0 + 1.5*bbox.x1, Y=bbox.y0 - 0.2*bbox.y0  , U=Ukey_vertical,
                  label=' {}m'.format(Ukey),labelpos='N',angle=90,labelcolor='cyan',fc='cyan',labelsep = 0.19,fontproperties={'size': 9.5,'weight':'bold'})


    # Set the colorbar limit
    # plot GPS stations using corresponding lon and lat coordinates
    #ax.scatter(X,Y,s=3,ec='k', c='k',linewidths=0.25,marker='.',
    #         transform=ccrs.PlateCarree(),label='GPS station')
    #ax.contourf(LONGITUDE,LATITUDE,E,transform=ccrs.Geodetic(),cmap=cmap)
    # set gridlines and outer geographical labels
  #gl = axes[1].gridlines(crs = ccrs.PlateCarree(),draw_labels=True, linewidth=0.25,linestyle='--')
  #gl.xlabels_top = False
  #gl.ylabels_right = False
  #gl.xlabel_style = {'size': 6}
  #gl.ylabel_style = {'size': 6}

  cardinal_labels = {"east": "", "west": "", "north": "", "south": ""}
  latitude_formatter = LatitudeFormatter(cardinal_labels=cardinal_labels)
  longitude_formatter = LongitudeFormatter(cardinal_labels=cardinal_labels)
  gl = axes[1].gridlines(
      draw_labels=["left", "bottom"],
      linestyle="--",
      xformatter=longitude_formatter,
      yformatter=latitude_formatter,
      ylabel_style={"size": 8},
      xlabel_style={"size": 8},
      linewidth=0.3
  )
  # x_coord,y_coord = [0.53,insetbox_y]
  # sub_ax = fig.add_axes([x_coord, y_coord, 0.15, 0.15],
  #                             projection=ccrs.PlateCarree())
  # extent = [min_lon -15 ,max_lon + 12,min_lat - 25 ,max_lat + 12]

  # sub_ax.set_extent(extent, geodetic)

  #   # Make a nice border around the inset axes.
  # #effect = Stroke(linewidth=1.5)
  # #sub_ax.spines['geo'].set_path_effects([effect])

  #   # Add the land, coastlines and the extent of the Solomon Islands.
  # sub_ax.add_feature(cfeature.LAND,color='lightgray')
  # sub_ax.add_feature(cfeature.BORDERS,linewidth=0.4)
  # sub_ax.coastlines(linewidth =0.4)
  # extent_box = sgeom.box(min_lon, min_lat, max_lon, max_lat)
  # sub_ax.add_geometries([extent_box], ccrs.PlateCarree(), facecolor='none',
  #                         edgecolor='blue', linewidth=1.25)

  nframes = (npoints//m)
  times = np.arange(0,256,m)
  def animate(i):
      q.set_UVC(Vts[:,i],Uts[:,i])
      q2.set_UVC(np.zeros_like(Zts[:,i]),Zts[:,i])
      #cax1.set_array(Zts[:,i])

      # arr = G0[:,:,i]
      # vmax     = np.max(arr)
      # vmin     = np.min(arr)

      # cf = ax.pcolormesh(x,y,arr, vmax=vmax, vmin=vmin)
      # fig.colorbar(cf, ax=ax)
      time.set_text('$t=$'+str(times[i])+' s')
  fig.suptitle(f'{name}',fontweight='bold',fontsize=14)
  plt.subplots_adjust()
  anim = animation.FuncAnimation(fig, animate,frames = nframes,interval=100)
  anim.save(filename=f'Dynamic_Simulations/{name}/{nsamples}_samples/{eq_name}_combined_disp_field_{nmodels}_with_slip.gif', writer="pillow")
  plt.close()  

m = 5
#animate_displacement('Pedernales',m,1,2.0,1,2.5,0.11,0.70,0.95)
#animate_displacement('Gorkha',m,2,4.65,1,4.25,0.365,0.55,0.95,loc='upper left')
animate_displacement('Iquique',m,1,3,1,2.5,0.11,0.6,0.95)
#animate_displacement('Illapel',m,4,5,1,2,0.11,0.65,0.95)
#animate_displacement('Tohoku',m,30,20,5,8,0.12,0.72,0.95)
