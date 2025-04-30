# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 13:30:17 2023

@author: joanv
"""



import numpy as np
import pandas as pd 

import os
from matplotlib.patches import Ellipse


def geom_reader(name:str,header=[0,1,4,5,6],nrows_skip = 1):
    '''
    

    Parameters
    ----------
    name : str
        DESCRIPTION.
    header : list
        contains list of header indeces in following order
        lon, lat, depth, strike, dip    
    Returns
    -------
    None.

    '''
    EQs_dir = os.path.join(os.getcwd(),'EQ')
    EQ_dir =  os.path.join(EQs_dir,name)
    geometry_dir = os.path.join(EQ_dir,'geometry')
    geometry_formatted_dir = os.path.join(geometry_dir,'geometry_formatted')
    geometry_files = os.listdir(geometry_dir)
    geometry_file_name = [ i for i in geometry_files if not i.endswith("formatted")][0]
    geometry_file_dir = os.path.join(geometry_dir, geometry_file_name) # must contain one file only
    
    
    ### Reading Geometry  ###
    try: 
        
        geometry  = pd.read_csv(geometry_file_dir,skiprows=nrows_skip,sep=' ') # subject to change            
        geometry = geometry[['lat','lon','depth','strike','dip']] # taking only columns of interest
        geometry_d = geometry.to_dict(orient='list')
    
    except: 
        
        geometry  = np.loadtxt(geometry_file_dir,skiprows=nrows_skip) # subject to change      
        i_lat = header[0]
        i_lon = header[1]
        i_depth = header[2]
        i_strike = header[3]
        i_dip = header[4]
        geometry_d =dict()
        geometry_d['lat'] =geometry[:,i_lat]
        geometry_d['lon'] = geometry[:,i_lon]
        geometry_d['depth'] = geometry[:,i_depth]
        geometry_d['strike'] = geometry[:,i_strike]
        geometry_d['dip'] = geometry[:,i_dip]
    
    csv_file_dir = os.path.join(geometry_formatted_dir,f'{name}_geometry.csv')
    geom_df = pd.DataFrame.from_dict(geometry_d, orient='index') # handle `Value error exception` 
    geom_df = geom_df.transpose()
    
    try:
        os.makedirs(geometry_formatted_dir)
    except FileExistsError:
        print(f'Directory {geometry_formatted_dir} already exist')
    
    geom_df.to_csv(csv_file_dir)
   
