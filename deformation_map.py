# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 23:34:16 2025

@author: joanv
"""


from matplotlib.patheffects import Stroke
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import shapely.geometry as sgeom
import numpy as np
import cartopy.io.img_tiles as cimgt
import pandas as pd
import h5py
from scipy.interpolate import RegularGridInterpolator
from matplotlib.colors import TwoSlopeNorm
import cartopy.io.shapereader as shpreader
import os

NATIONS = {'Tohoku':'JPN',
           'Iquique':'CHL',
           'Illapel':'CHL',
           'Pedernales':'ECU',
           'Gorkha':'NPL'}
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
trench_name = ['kuriljapan','southamerica','southamerica','southamerica','MHT']
limits_model = [[138,145,35,44],[-73,-69,-21,-18],[-74,-70,-32,-29],[-82,-79,-1,1],[84,87,27,29]]
limits_disp = [[138-1,145+2,35-2,44+2],[-72-1,-69+1,-21-1,-18+1],[-72-2,-70+1,-33-1,-29+1],[-81-1,-80+1,-1-0.5,1+0.5],[84-0.25,87,27-0.25,29+0.25]]
nsamples = 100
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,trench_name,limits_model,limits_disp):
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
        model[name]['limits_model'] = limits_model[i]
        model[name]['limits_disp'] = limits_disp[i]
    return model

models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,trench_name,limits_model,limits_disp)

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
  
    


geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
# geographical projection used: 'PlateCarree'
proj = ccrs.PlateCarree()

working_dir = os.getcwd()

names = ['Tohoku','Illapel','Iquique','Pedernales','Gorkha']

positions = {'Illapel':(1,0),'Iquique':(2,0),'Pedernales':(1,1),'Gorkha':(2,1),'Tohoku':(0,0)}
Ukeys = {'Iquique':1,'Illapel':3,'Pedernales':2,'Gorkha':2,'Tohoku':30}
inset_coords = {'Iquique':[0.672,0.748],'Illapel':[0.393,0.772],'Pedernales':[0.395,0.364],'Gorkha':[0.78,0.335],'Tohoku':[0.115,0.59]}
scales = {'Iquique':4,'Illapel':7,'Pedernales':4,'Gorkha':9,'Tohoku':18}
shrinks = {'Iquique':0.6,'Illapel':0.6,'Pedernales':0.6,'Gorkha':0.5,'Tohoku':0.3}
sep_contours = {'Iquique':0.25,'Illapel':0.5,'Pedernales':0.25,'Gorkha':0.5,'Tohoku':2}
sep_maps = {'Iquique':0.1,'Illapel':0.1,'Pedernales':0.05,'Gorkha':0.05,'Tohoku':0.25}


fig = plt.figure(figsize=(15,10))
gs0 = fig.add_gridspec(1, 3,wspace=0.25)

gs00 = gs0[0].subgridspec(1, 1,height_ratios=[10],hspace=0.05)
gs01 = gs0[1].subgridspec(2,1,hspace=0.05)
gs02 = gs0[2].subgridspec(2,1,hspace=0.05)
gs = [gs00,gs01,gs02]
for l,name in enumerate(names):
  position = positions[name]
  Ukey = Ukeys[name]
  inset_coord = inset_coords[name]
  scale = scales[name]
  shrink = shrinks[name]
  sep_contour = sep_contours[name]
  sep_map = sep_maps[name]

  if name == 'Tohoku':
    ax = fig.add_subplot(gs[0][0],projection = proj)
  else:
    ax = fig.add_subplot(gs[positions[name][0]][positions[name][1],0],projection = proj)
  limits = models[name]['limits_disp']
  method = 'linear'
  nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
  patch = models[name]['patch']
  print(patch)
  df = pd.read_csv(f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')
  strike = np.mean(df['strike'].values)
  X = df['lon'].values
  Y = df['lat'].values
  U_parallel = df['U_parallel'].values
  U_perp = df['U_perp'].values
  Slip = np.sqrt(U_parallel**2 + U_perp**2)
  std_U_parallel = df['std_U_parallel'].values
  std_U_perp = df['std_U_perp'].values
  hypo_as = df['Hypo_as'].values
  hypo_dd = df['Hypo_dd'].values
  dictionary = {'U_parallel':U_parallel,'U_perp':U_perp,'std_U_parallel':std_U_parallel,'std_U_perp':std_U_perp}
  trench_name = models[name]['trench']
  shape = (nrows,ncols)


  km = 1000
  xstn,ystn = set_stn(2,4,patch,shape,control=1)

  Xstn,Ystn = np.meshgrid(xstn,ystn)
  xstn_flat,ystn_flat =  Xstn.flatten()*km,Ystn.flatten()*km


  LON = two_array_formatter(X,shape)
  LAT = two_array_formatter(Y,shape)
  print(LON.shape)



  patch_interp = patch/3


  xS_interp = xstn
  yS_interp = ystn
  XS_interp,YS_interp  = np.meshgrid(xS_interp,yS_interp)

  xS = np.arange(patch/2 , ncols*patch,patch)
  yS_unproj = np.arange(-(nrows-1/2)*patch,0,patch)
    # shift accordingly at surface
  yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
  XS, YS = np.meshgrid(xS,yS)



  lat_interp = RegularGridInterpolator((yS, xS), LAT,bounds_error=False, fill_value=None,method = method)
  lon_interp = RegularGridInterpolator((yS, xS), LON,bounds_error=False, fill_value=None,method = method)

  LON = lon_interp((YS_interp,XS_interp))
  LAT = lat_interp((YS_interp,XS_interp))



  LON = LON.flatten()
  LAT = LAT.flatten()



  h5file_dir = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS/{nsamples}_samples')
  h5f_name_edks = os.path.join(h5file_dir,f'EDKS_{name}_displacement_nsamples_{nsamples}_parallel.h5')



  h5file_edks = h5py.File(h5f_name_edks,'r')



  size = h5file_edks['displacement'].shape[1]
  E_edks = np.mean(h5file_edks['displacement'][:,:int(size/3)],axis=0)
  N_edks = np.mean(h5file_edks['displacement'][:,int(size/3):2*int(size/3)],axis=0)
  Z_edks = np.mean(h5file_edks['displacement'][:,2*int(size/3):],axis=0)


  trench_dir = os.path.join(working_dir,f'Trench_TPGA/{trench_name}.lonlat')
  trench = np.loadtxt(trench_dir)
  trench_lon = trench[:,0]
  trench_lat = trench[:,1]
  Z=Z_edks

  # next two lines give error   
  #google_terrain = cimgt.GoogleTiles(style="satellite")
  #ax.add_image(google_terrain, '10')

  min_lon,max_lon,min_lat,max_lat = limits

  # set grid geographical boundaries for aoi ie. [min_lon, max_Lon,min_lat,max_lat]
  ax.set_extent([min_lon, max_lon, min_lat, max_lat])

  ax.coastlines(linewidth =0.8)
  ax.add_feature(cfeature.BORDERS,linewidth=0.7)
  # set ocean and land features
  ax.add_feature(cfeature.LAND, color='lightgray')
  ax.add_feature(cfeature.OCEAN, color='skyblue')
  if name=='Gorkha':
      markevery = 5
  else:
      markevery = 25
  ax.plot(trench_lon,trench_lat,color='black',lw=0.4,marker=">",markevery=markevery,ms=0.8)


  # set US state boudaries
  #ax.add_feature(cfeature.STATES, linewidth=0.75)
  cmap = plt.get_cmap('bwr') # Choose your desired colormap
  norm = TwoSlopeNorm(0,vmin=(min(Z.min(),-Z.max())),vmax=(max(-Z.min(),Z.max())))
  #ax.tricontour(LON, LAT, Z ,levels=np.arange(np.floor(min(Z.min(),-Z.max())),np.ceil(max(-Z.min(),Z.max())) + sep_contour,sep_contour),inline=True, linewidths=0.15, colors='k')
  im = ax.tricontourf(LON,LAT,Z,transform=ccrs.PlateCarree(),cmap=cmap,norm=norm,levels = np.arange(round(min(Z.min(),-Z.max()),1),round(max(-Z.min(),Z.max()),1) + sep_map,sep_map))
  # Add a colorbar
  cbar = fig.colorbar(im,ax=ax,shrink = shrink,cmap=cmap,pad = 0.15,norm=norm,orientation ='vertical')
  cbar.set_label(label='$d_z$ (m)')
  cbar.ax.tick_params(labelsize=8)
  U = N_edks
  V = E_edks   # strike-slip displacement (see p.4 Minson et al. (2013) Part II)
  displacement = np.sqrt(U**2 + V**2)
  strike  = strike*np.pi/180
  strike = (strike - 90*np.pi/180)
  Ux = V*np.cos(strike) + U*np.sin(strike)
  Uy = -V*np.sin(strike) + U*np.cos(strike)
  #Ux = V*np.cos(strike - np.pi/2) + U*np.sin(np.pi - strike)
  #Uy = V*np.cos(strike) + U*np.cos(strike - np.pi/2)



  q = ax.quiver(LON,LAT,Ux,Uy,scale=scale,scale_units ='x', color='black',units='width',width=0.00225)

  bbox = ax.get_position()
  ax.quiverkey(q, X=bbox.x0 + 0.1*bbox.x1, Y=bbox.y0 - 0.79*bbox.y0  , U=Ukey,
                label=' {}m'.format(Ukey),labelpos='N',angle=0,labelcolor='black',fc='black',labelsep = 0.015,fontproperties={'size': 8,'weight':'bold'})
  factor = 1/q.scale    # scaling factor from UV unit to XY units

  ax.text(0.34,1.15,f'{name}',fontsize=14,fontweight='bold',transform=ax.transAxes)
  ax.text(0.010,1.15,'abcde'[l],fontsize=14,fontweight='bold',transform=ax.transAxes)
  # Set the colorbar limit
  # plot GPS stations using corresponding lon and lat coordinates
  #ax.scatter(X,Y,s=3,ec='k', c='k',linewidths=0.25,marker='.',
  #         transform=ccrs.PlateCarree(),label='GPS station')
  #ax.contourf(LONGITUDE,LATITUDE,E,transform=ccrs.Geodetic(),cmap=cmap)
  # set gridlines and outer geographical labels
  gl = ax.gridlines(crs = ccrs.PlateCarree(),draw_labels=True, linewidth=0.1,linestyle='--')
  gl.xlabels_top = False
  gl.ylabels_right = False
  gl.xlabel_style = {'size': 6}
  gl.ylabel_style = {'size': 6}

  x_coord,y_coord = inset_coord
  sub_ax = fig.add_axes([x_coord, y_coord, 0.06, 0.06],
                            projection=ccrs.PlateCarree())
  extent = [min_lon -15 ,max_lon + 12,min_lat - 20 ,max_lat + 5]

  sub_ax.set_extent(extent, geodetic)

  # Make a nice border around the inset axes.
  effect = Stroke(linewidth=1.0)
  sub_ax.spines['geo'].set_path_effects([effect])

  # Add the land, coastlines and the extent of the Solomon Islands.
  sub_ax.add_feature(cfeature.LAND,color='lightgray')
  sub_ax.add_feature(cfeature.BORDERS,linewidth=0.4)
  sub_ax.coastlines(linewidth =0.4)
  extent_box = sgeom.box(min_lon, min_lat, max_lon, max_lat)
  shpfilename = shpreader.natural_earth(resolution='110m',
                                    category='cultural',
                                    name='admin_0_countries')
  reader = shpreader.Reader(shpfilename)
  countries = reader.records()
  print(NATIONS[name])
  for country in countries:
      if country.attributes['ADM0_A3'] in NATIONS[name]:
          sub_ax.add_geometries(country.geometry, ccrs.PlateCarree(), facecolor='red')
  sub_ax.add_geometries([extent_box], ccrs.PlateCarree(), facecolor='none',
                        edgecolor='blue', linewidth=0.7)

plt.savefig('all_deformation_edks_countries_red.pdf')
  