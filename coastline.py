# -*- coding: utf-8 -*-
"""
Created on Thu May 15 13:02:37 2025

@author: jvilo
"""



import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
import re

import cartopy.io.shapereader as shpreader

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
lon_lim = [(139,144.5),(-74,-68),(-75,-69),(-83,-79),(80,89)]
lat_lim = [(33,45.5),(-22,-17),(-33.5,-28),(-1.2,1.2),(26,32)]

country = [['JPN'],['CHL','PER'],['CHL'],['ECU'],['NPL']]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,trench_name,lon_lim,lat_lim,country):
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
        model[name]['country'] = country[i]
    return model

models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,hspaces,wspaces,sizes,shrinks,trench_name,lon_lim,lat_lim,country)
name = 'Tohoku'

shpfilename = shpreader.natural_earth(resolution='10m',
                                  category='cultural',
                                  name='admin_0_countries')
reader = shpreader.Reader(shpfilename)
countries = reader.records()
target_countries = models[name]['country'] 
geoms = []
for country in countries:
    for target_country in target_countries:
        if country.attributes['ADM0_A3'] == target_country :
            geoms.append(country.geometry)
latitude_coast = np.array([])
longitude_coast =  np.array([])           
for geom in geoms:        
    poly = str(geom)
    matches = re.findall('-?\d+.\d{14}',poly)
    
    coords = np.array(matches,dtype='float')
    
    if name in ('Tohoku'):
        lon = np.array([coords[i] for i in range(0,len(coords),2)])[:-1]  
        lat = np.array([coords[i] for i in range(1,len(coords),2)])

    else:
        lon = np.array([coords[i] for i in range(0,len(coords),2)])   
        lat = np.array([coords[i] for i in range(1,len(coords),2)])
    ids_problem = []
    
    for i, latc in enumerate(lat):
        if latc<models[name]['lat_lim'][0] or latc>models[name]['lat_lim'][1]:
            ids_problem.append(i)
    
    lon_copy = np.copy(lon)
    lat_copy = np.copy(lat)
    
    lon[ids_problem] = lat_copy[ids_problem]
    lat[ids_problem] = lon_copy[ids_problem]
    
    
    ids_bounds_lon = np.where((lon<models[name]['lon_lim'][1]) & (lon>models[name]['lon_lim'][0]))                
    lon_coast = lon[ids_bounds_lon]
    lat_coast = lat[ids_bounds_lon]
    
    ids_bounds_lat = np.where((lat_coast<models[name]['lat_lim'][1]) & (lat_coast>models[name]['lat_lim'][0]))                
    lon_coast = lon_coast[ids_bounds_lat]
    lat_coast = lat_coast[ids_bounds_lat]
    
    latitude_coast = np.concatenate((latitude_coast,lat_coast))
    longitude_coast = np.concatenate((longitude_coast,lon_coast))
   
boundaries = np.column_stack((longitude_coast,latitude_coast))
plt.plot(longitude_coast,latitude_coast)
# for i in range(len(boundaries)):
#     plt.annotate(f'{i}',(boundaries[i][0],boundaries[i][1]),fontsize=1.5)
# plt.savefig('tohoku.pdf')
os.makedirs(os.path.join(os.getcwd(),'boundaries'),exist_ok=True)

if name=='Tohoku':
    # verified manually
    longitude_coast1 = longitude_coast[:616]
    longitude_coast2 = longitude_coast[624:850]
    
    latitude_coast1 = latitude_coast[:616]
    latitude_coast2 = latitude_coast[624:850]
    
    boundaries1 = np.column_stack((longitude_coast1,latitude_coast1))
    boundaries2 = np.column_stack((longitude_coast2,latitude_coast2))
    

    np.savetxt(f'boundaries/{name}_coast1.lonlat',boundaries1)
    np.savetxt(f'boundaries/{name}_coast2.lonlat',boundaries2)
else:
    np.savetxt(f'boundaries/{name}_coast.lonlat',boundaries)
    

    

        
        
        
        
        
        