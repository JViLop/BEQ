# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 10:39:38 2024

@author: joanv
"""

import numpy as np
import pandas as pd
from scipy.interpolate import RegularGridInterpolator


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
hspaces = [-0.25,0.10,0.10,0.10,0.10]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['hspace'] = hspaces[i]
    return model

km =1000

name = 'Tohoku'


models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces)

for name in names:
    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    patch = models[name]['patch']
    df = pd.read_csv(f'INPUT/{name}/model/kinematic/500_samples/mean/{name}_mean_kinematic_model.csv')
    hypo_dd = -df['Hypo_dd'].values[0]
    hypo_as = df['Hypo_as'].values[0]
    strike = np.mean(df['strike'].values)
    dip = df['dip'].values[:nrows][2]
    def proj_ysrc_coords(patch,dip):
          proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
          proj_ysrc = np.zeros_like(proj_dysrc)
          for i in range(len(proj_ysrc)):
              proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
          ysrc = np.flip(proj_ysrc)
            
          return ysrc
      
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
    
    xS = np.arange(patch/2 , ncols*patch,patch)
    yS = np.arange(-(nrows-1/2)*patch,0,patch)
       # shift accordingly at surface
    #yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
    XS, YS = np.meshgrid(xS,yS)
    
    
    DIP = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0)
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0) 
    
    LON = np.flip(df['lon'].values.reshape(nrows,ncols,order='F'),axis=0)
    LAT = np.flip(df['lat'].values.reshape(nrows,ncols,order='F'),axis=0)
    
    method = 'linear'
    LON_interp = RegularGridInterpolator((yS, xS), LON,bounds_error=False, fill_value=None,method = method)
    LAT_interp = RegularGridInterpolator((yS, xS), LAT,bounds_error=False, fill_value=None,method = method)
    
    lon_hyp = LON_interp((hypo_dd,hypo_as))
    lat_hyp = LAT_interp((hypo_dd,hypo_as))
    
    print(f'{name}',lon_hyp,lat_hyp)
    

