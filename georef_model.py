# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:04:12 2025

@author: joanv
"""

from scipy.interpolate import RegularGridInterpolator
import pandas as pd

import os
import numpy as np

from matplotlib.patheffects import Stroke
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import shapely.geometry as sgeom

import matplotlib as mpl
from matplotlib.patches import Ellipse
import matplotlib.patches as patches
from matplotlib.path import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable  
import cartopy.io.shapereader as shpreader
working_dir = os.getcwd()

def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)
    
def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()
    
def set_stn(c,factor,patch,control=0):
      
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



km = 1000


NATIONS = {'Tohoku':'JPN',
           'Iquique':'CHL',
           'Illapel':'CHL',
           'Pedernales':'ECU',
           'Gorkha':'NPL'}
names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
magnitudes = ['9.0','8.1','8.3','7.8','7.8']
refs = ['(Minson et al. 2014)','(Duputel et al. 2015)','(Caballero et al. 2023)','(Gombert et al. 2018)','(Yue et al. 2017)']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patch_sizes = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [0,0.5,0.5,0.5,0.5]
wspaces = [0.2,0.3,0.4,0.3,0.4]
sizes = [(14,16),(10,11),(10,12),(10,11),(14,16)]
shrinks = [0.5,0.5,0.5,0.5,0.2]
limits_model = [[139.5,145.5,35,43],[-73.5,-69,-21.5,-18.5],[-74,-70,-33,-29],[-81.5,-79,-1,1],[84,86.5,27,29.25]]# [[137,146,34,45],[-74,-68,-22,-17],[-72,-68,-33,-28],[-82,-79,-1,1],[83,88,26,30]]
limits_disp = [[138-2,145+2,35-2,44+2],[-73-1,-69+1,-21-1,-18+1],[-74-1,-70-1,-32-1,-29+1],[-82-1,-79+1,-1-1,1+1],[84-0.25,87,27-0.25,29+0.25]]
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

models = model_dict(names,geoms,patch_sizes,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,limits_model,limits_disp)





positions = {'Illapel':(1,0),'Iquique':(2,0),'Pedernales':(1,1),'Gorkha':(2,1),'Tohoku':(0,0)}
Ukeys = {'Iquique':8,'Illapel':10,'Pedernales':4,'Gorkha':4,'Tohoku':30}
inset_coords = {'Iquique':[0.675,0.726],'Illapel':[0.395,0.777],'Pedernales':[0.395,0.36],'Gorkha':[0.84,0.38],'Tohoku':[0.115,0.64]}
scales = {'Iquique':20,'Illapel':40,'Pedernales':20,'Gorkha':20,'Tohoku':60}
shrinks = {'Iquique':0.6,'Illapel':0.6,'Pedernales':0.6,'Gorkha':0.5,'Tohoku':0.3}
sep_contours = {'Iquique':0.25,'Illapel':0.5,'Pedernales':0.25,'Gorkha':0.5,'Tohoku':2}
sep_maps = {'Iquique':0.1,'Illapel':0.1,'Pedernales':0.05,'Gorkha':0.05,'Tohoku':0.25}
key_positions = {'Iquique':[0.1,0.1],'Illapel':[0.1,0.1],'Pedernales':[0.1,0.1],'Gorkha':[0.9,0.1],'Tohoku':[0.85,0.05]}
percentage = {'Iquique':0.93,'Illapel':0.92,'Pedernales':0.94,'Gorkha':0.97,'Tohoku':0.94}
fig = plt.figure(figsize=(15,10))
gs0 = fig.add_gridspec(1, 3,wspace=0.3)

gs00 = gs0[0].subgridspec(1, 1,height_ratios=[15],hspace=0.02)
gs01 = gs0[1].subgridspec(2,1,hspace=0.02)
gs02 = gs0[2].subgridspec(2,1,hspace=0.02)
gs = [gs00,gs01,gs02]
proj = ccrs.PlateCarree()

for l,name in enumerate(names):




    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    nparam = models[name]['nparam']

    shape = (nrows,ncols)
    patch = models[name]['patch']





        
    df = pd.read_csv(f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')



    LON = one_array_formatter(df['lon'].values,shape)
    LAT = one_array_formatter(df['lat'].values,shape)

    ONSET = np.loadtxt(f'RuptureTime_{name}_all_kinematic_model.txt')



    ONSET = ONSET.reshape(shape)
    # df.insert(8, "onset", onset, True)
    # df.insert(5,"type",rise_type,True)

    # L = np.ones_like(onset)*patchsize*1e3
    # W = np.ones_like(onset)*patchsize*1e3

    # df.insert(10, "L(m)", L, True)
    # df.insert(11,"W(m)",W,True)
    percent = percentage[name]

    Uperp = two_array_formatter(df['U_perp'].values,shape)
    Uparallel =  two_array_formatter(df['U_parallel'].values,shape)
    SLIP = two_array_formatter(df['Slip'].values,shape)
    RAKE = two_array_formatter(np.arctan2(Uparallel,Uperp)*(180/np.pi),shape)
    DIP = two_array_formatter(df['dip'].values,shape)
    # vertical shift for seafloor
    DEPTH = two_array_formatter(df['depth'].values,shape)
    # DEPTH = array_formatter(df['depth'].values,shape) 
    DURATION = two_array_formatter(df['Tr'].values,shape)

    LON = two_array_formatter(df['lon'].values,shape)
    LAT = two_array_formatter(df['lat'].values,shape)
    STRIKE = two_array_formatter(df['strike'].values,shape)



    patch_interp = patch/10

    xS_interp = np.arange(patch_interp , ncols*patch,patch_interp)
    yS_interp = np.arange(-(nrows*percent*patch)-patch_interp, patch_interp/8,patch_interp)
    XS_interp,YS_interp  = np.meshgrid(xS_interp,yS_interp)

    nrows_interp,ncols_interp = len(yS_interp),len(xS_interp)

    xS = np.arange(patch/2 , ncols*patch,patch)
    yS = np.arange(-(nrows-1/2)*patch,0,patch)
       # shift accordingly at surface
    yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
    XS, YS = np.meshgrid(xS,yS)


    DIP = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0)
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0) 
    method = 'linear'
    dip_interp = RegularGridInterpolator((yS, xS), DIP,bounds_error=False, fill_value=None,method = method)
    strike_interp = RegularGridInterpolator((yS, xS), STRIKE,bounds_error=False, fill_value=None,method = method)

    depth_interp = RegularGridInterpolator((yS, xS), DEPTH,bounds_error=False, fill_value=None,method = method)

    lat_interp = RegularGridInterpolator((yS, xS), LAT,bounds_error=False, fill_value=None,method = method)
    lon_interp = RegularGridInterpolator((yS, xS), LON,bounds_error=False, fill_value=None,method = method)
    onset_interp = RegularGridInterpolator((yS, xS), ONSET,bounds_error=False, fill_value=None,method = method)

    DIP = dip_interp((YS_interp,XS_interp))
    DEPTH = depth_interp((YS_interp,XS_interp))
    LON = lon_interp((YS_interp,XS_interp))
    LAT = lat_interp((YS_interp,XS_interp))
    ONSET = onset_interp((YS_interp,XS_interp))


    lat_hyp = lat_interp((-df['Hypo_dd'].values[0],df['Hypo_as'].values[0]))
    lon_hyp = lon_interp((-df['Hypo_dd'].values[0],df['Hypo_as'].values[0]))
    position = positions[name]
    Ukey = Ukeys[name]
    inset_coord = inset_coords[name]
    scale = scales[name]
    shrink = shrinks[name]
    sep_contour = sep_contours[name]
    sep_map = sep_maps[name]
    ukey_pos = key_positions[name]
    if name == 'Tohoku':
        ax = fig.add_subplot(gs[0][0],projection = proj)
    else:
        ax = fig.add_subplot(gs[positions[name][0]][positions[name][1],0],projection = proj)
    data = df 

    geodetic = ccrs.Geodetic(globe=ccrs.Globe(datum='WGS84'))
    # geographical projection used: 'PlateCarree'
    
    
    
    min_lon,max_lon,min_lat,max_lat = models[name]['limits_model']
    print(models[name]['limits_model'])
    # set grid geographical boundaries for aoi ie. [min_lon, max_Lon,min_lat,max_lat]
    ax.set_extent([min_lon, max_lon, min_lat, max_lat])
    
    ax.coastlines(linewidth =0.4)
    ax.add_feature(cfeature.BORDERS,linewidth=0.8)
    patch_size = models[name]['patch']
    
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    X = data['lon'].values
    Y = data['lat'].values
    
    X = data['lon'].values
    Y = data['lat'].values
    U_parallel = data['U_parallel'].values
    U_perp = data['U_perp'].values
    Slip = np.sqrt(U_parallel**2 + U_perp**2)
    
    std_U_parallel = data['std_U_parallel'].values
    std_U_perp = data['std_U_perp'].values
    hypo_as = data['Hypo_as'].values
    hypo_dd = data['Hypo_dd'].values
    
    dictionary = {'U_parallel':U_parallel,'U_perp':U_perp,'std_U_parallel':std_U_parallel,'std_U_perp':std_U_perp}
    strike = np.mean(data['strike'].values)
    along_strike  = strike*np.pi/180
    along_dip = (along_strike - 90*np.pi/180)
    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]
    
    km_2_degree = 1/111.
    
    move_ad = (patch_size/2)*km_2_degree*np.array([np.sin(along_dip),np.cos(along_dip)])
    move_as = (patch_size/2)*km_2_degree*np.array([np.sin(along_strike),np.cos(along_strike)])
    
    
    cmap = plt.get_cmap('YlOrRd')
    
    # Normalize data to [0, 1] range
    norm = plt.Normalize(Slip.min(), Slip.max())
    
    # Get colors for each value in the array
    colors = cmap(norm(Slip))
   
    vertices =  {i:[np.array([X[i],Y[i]]) + move_ad + move_as,np.array([X[i],Y[i]]) - move_ad + move_as,np.array([X[i],Y[i]]) - move_ad - move_as, np.array([X[i],Y[i]]) + move_ad - move_as,np.array([X[i],Y[i]]) + move_ad + move_as] for i in np.arange(len(X))}
    for i in np.arange(len(X)):
      vertex = []
      for k in np.arange(5):
        vertex.append((vertices[i][k][0],vertices[i][k][1]))
      path = Path(vertex)
      #patch = patches.PathPatch(path,facecolor=colors[i],edgecolor='black',lw=1e-3)
      patch = patches.PathPatch(path,facecolor=colors[i],edgecolor=colors[i],lw=None)
      ax.add_patch(patch)
      ax.set_xlabel('Longitude')
      ax.set_ylabel('Latitude')
      ax.set_aspect('equal')

    
    im = ax.tricontourf(LON.flatten(),LAT.flatten(),ONSET.flatten(),alpha=0,colors='blue',linewidths=5,levels = 12)
    
    for line in im.collections:
        paths = line.get_paths()
        for path in paths:
          #patch = patches.PathPatch(path, edgecolor = 'blue',facecolor='none', lw=0.1,ls='-')
          patch = patches.PathPatch(path, edgecolor = 'blue',facecolor='none', lw=0.3,ls='-')
          ax.add_patch(patch)
    
          #vertices = path.vertices
          #x_coords, y_coords = vertices[:, 0], vertices[:, 1]
    
            # Now you have the x and y coordinates of each contour point
            # Do something with them, e.g., print them
        #ax.plot(x_coords,y_coords,color='k',lw=0.25)
    
    #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
    ax.margins(0)
    U = dictionary['U_parallel'] # dip-slip displacement
    V = dictionary['U_perp']    # strike-slip displacement (see p.4 Minson et al. (2013) Part II)
    
    strike = along_dip
    Ux = V*np.cos(strike) + U*np.sin(strike)
    Uy = -V*np.sin(strike) + U*np.cos(strike)
    q = ax.quiver(X,Y,Ux,Uy,scale=scale,scale_units ='x', color='black',units='width',edgecolors='black',width=0.005,headwidth=3,headlength=5.0)
    if name =='Pedernales':
        off_x = 0.05
    else:
        off_x = 0
        
    ax.text(0.3 - off_x,1.17,' $M_{w}$ ' + f'{magnitudes[l]} ' + f'{name}',fontsize=13.5,fontweight='bold',transform=ax.transAxes)
    #ax.text(0.21,1.1, f' {refs[l]}',fontsize=13,fontweight='bold',transform=ax.transAxes)

    ax.text(-0.025,1.1,'acbde'[l],fontsize=15,fontweight='bold',transform=ax.transAxes)
    bbox = ax.get_position()
    
    ax.quiverkey(q, X=ukey_pos[0],Y=ukey_pos[1], U=Ukey,
                  label=' {} m'.format(Ukey), color='black',labelpos='N',angle=0)
    factor = 1/q.scale    # scaling factor from UV unit to XY units
    
    offsetXY = np.column_stack((Ux,Uy))
    
    
    ax.coastlines(linewidth =0.8)
    ax.add_feature(cfeature.BORDERS,linewidth=0.7)
    # set ocean and land features
    ax.add_feature(cfeature.LAND, color='lightgray')
    ax.add_feature(cfeature.OCEAN, color='skyblue')
    
    offsetXY = np.array([[xy[0],xy[1]] for xy in offsetXY])
    coord = q.XY + offsetXY*factor
    # width and height had to be formatted in row-like order (initially in column-like order)
    ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
                    width=dictionary['std_U_perp'][i]*factor,
                    height=dictionary['std_U_parallel'][i]*factor,
                    angle=0,alpha=0.5,fill=False,edgecolor='black',lw=0.5)
            for i in range(q.N)]
    for e in ells:
        ax.add_artist(e)
    ax.plot(lon_hyp,lat_hyp,marker='*',ms=18,c='blue',markeredgecolor='k')
    
    sub_ax = fig.add_axes([inset_coord[0], inset_coord[1], 0.07, 0.07],
                              projection=ccrs.PlateCarree())
    extent = [min_lon -15 ,max_lon + 12,min_lat - 20 ,max_lat + 5]
    gl = ax.gridlines(crs = ccrs.PlateCarree(),draw_labels=True, linewidth=0.1,linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabel_style = {'size': 7.5}
    gl.ylabel_style = {'size': 7.5}

    sub_ax.set_extent(extent, geodetic)
    
    # Make a nice border around the inset axes.
    effect = Stroke(linewidth=1.5)
    sub_ax.spines['geo'].set_path_effects([effect])
    
    # Add the land, coastlines and the extent of the Solomon Islands.
    sub_ax.add_feature(cfeature.LAND,color='lightgray')
    sub_ax.add_feature(cfeature.BORDERS,linewidth=0.4)
    sub_ax.coastlines(linewidth =0.4)
    shpfilename = shpreader.natural_earth(resolution='110m',
                                      category='cultural',
                                      name='admin_0_countries')
    reader = shpreader.Reader(shpfilename)
    countries = reader.records()
    
    for country in countries:
        if country.attributes['ADM0_A3'] in NATIONS[name]:
            sub_ax.add_geometries(country.geometry, ccrs.PlateCarree(), facecolor='red')
    extent_box = sgeom.box(min_lon, min_lat, max_lon, max_lat)
    sub_ax.add_geometries([extent_box], ccrs.PlateCarree(), facecolor='none',
                          edgecolor='blue', linewidth=0.8)
    
   
    
    divider = make_axes_locatable(sub_ax)  

     

    cax = divider.append_axes('bottom', '10%', pad='12%', axes_class=mpl.pyplot.Axes)
    sub_ax.get_figure().colorbar(mpl.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax,label='Slip (m)',orientation='horizontal')

    cax.xaxis.set_tick_params(labelsize=6.5)
plt.savefig('all_models_countries_no_contour.pdf')