import os
import numpy as np
import os
from os.path import join,exists,basename
import math
from glob import glob

# UTM Projection
lon0 = -80.0
lat0 =   0.0

# Hypocenter (USGS)
# lat_hypo = 0.382
# lon_hypo = -79.922
# dep_hypo =  20.6

# Hypocenter (Nocquet at al. 2016)
lat_hypo = 0.35
lon_hypo = -80.17
dep_hypo =  17

# Green's functions calculation using EDKS
plot_figs   = True
comp_GFs    = False
comp_KinGFs = True 
comp_bigG   = True
comp_kinCp  = False

# EDKS 
ref_model = 'model2.edks'

# Fault
faultGeometry = os.path.expanduser('~/Projects/Ecuador/geometry4/Ecuador_8x10.lonlat')
Nstrike = 10
Ndip    = 8
pSize   = 15.
Np = Nstrike*Ndip
RotAngle = -99.
if RotAngle < 0.:
    RotAngle += 360.

# InSAR
varresSDes=os.path.expanduser('~/Projects/Ecuador/data/insar/descending_S1/Descending_Sentinel')
varresADes=os.path.expanduser('~/Projects/Ecuador/data/insar/descending_ALOSv4/Descending_ALOS')
varresAAsc=os.path.expanduser('~/Projects/Ecuador/data/insar/ascending/THRE0.45/Ascending_ALOS')

# Tsunami
dart_dir = os.path.expanduser('~/Projects/Ecuador/data/tsunami3_8x10')
dart_wav = join(dart_dir,'tsunami_cut')
dart_info = join(dart_dir,'station_info.txt')
dart_nobs_per_trace = 66

# Static GPS
data_cgps = os.path.expanduser('~/Projects/Ecuador/data/cgps/cgps.dat')
data_hrgps = os.path.expanduser('~/Projects/Ecuador/data/cgps/hrgps.dat')
data_gps = os.path.expanduser('~/Projects/Ecuador/data/gps/gps.dat')

# GFs
staticGFdir = 'staticGFs'
GFdir = 'KinGFs'

# Strong motion data
seis_dir_sm       = os.path.expanduser('~/Projects/Ecuador/data/SM/Meter_filt_decim/')
sac_lst_sm        = join(seis_dir_sm,'o_inv_sac_file_lst4')

# HRGPS data
seis_dir_gps       = os.path.expanduser('~/Projects/Ecuador/data/HRGPS/Meter_filt_decim/')
sac_lst_gps        = join(seis_dir_gps,'o_inv_sac_file_lst')

# GF names
rakes_key = ['RP','AR']


# Green's functions for seismic data
earth_model = os.path.expanduser('~/Projects/Ecuador/Kernels_cinematique/ecuador.model2')
gf_db_dir   = os.path.expanduser('~/Projects/Ecuador/Kernels_cinematique/GF_PEDERNALES_run3')
gf_db_scale = 1.0e-28 # 1.0e-26 in cm -> 1.0e-28 in m
Ntriangles  = 70
Dtriangles  = 1.
Npt         = 4
Nmesh       = 4

