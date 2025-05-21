# -*- coding: utf-8 -*-
"""
Created on Wed May 14 11:23:52 2025

@author: jvilo
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

for name in ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']:
# for name in ['Gorkha']:
    os.system(f'python utils/run_edks_mean_errors_parallel.py {name}')
    

