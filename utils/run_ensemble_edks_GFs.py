# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 09:59:17 2024

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 16:57:56 2024

@author: joanv
"""

# if __name__ =='__main__':
#from mpi4py import MPI
import numpy as np 
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
import sys
import os
import h5py
working_dir = os.getcwd()
# Okada
import okada4py as ok92

# EDKS
from layered_disloc import layered_disloc



        

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
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,nramps,factors):
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
    return model

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
   



    
models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,ramps,factors)

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()
#nprocs = comm.Get_size()

name = str(sys.argv[1])
km = 1000
nsamples = 100
factor = models[name]['factor']
    # units
# in meters
# definitions
edks = f'{name}.edks'
prefix = name
#edks = 'tohoku_nied_2_shifted2.8km.edks'
#prefix = 'tohoku_nied_2_shifted2.8km'
#edks = 'tohoku_layered_removed.edks'
#prefix = 'tohoku_layered_removed'
BIN_EDKS = '${HOME}/EDKS_py/bin'
input_dir  = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')

nrows, ncols = models[name]['geom'][0], models[name]['geom'][1]
patch = models[name]['patch']
nparam = models[name]['nparam']

   

df = pd.read_csv(input_dir)
df= df.drop(df.columns[0],axis=1)

xstn,ystn = set_stn(2,4,patch,control=1)

Xstn,Ystn = np.meshgrid(xstn,ystn)
xstn_flat,ystn_flat =  Xstn.flatten()*km,Ystn.flatten()*km
zstn_flat = np.zeros(xstn_flat.shape)
model_type = 'kinematic'


file_folder = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/bin_data')
file_dir = os.path.join(file_folder,f'{name}_{model_type}_n_{nsamples}.dat')
bin_data = np.fromfile(file_dir,'float').reshape((nparam,nsamples))
nsteps = bin_data.shape[1]

h5f_name_edks =  f'EDKS_{name}_displacement_nsamples_{nsamples}.h5'
h5f_name_okada =  f'Okada_{name}_displacement_nsamples_{nsamples}.h5'
h5f_folder_edks = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS/{nsamples}_samples')
h5f_folder_okada = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/Okada/{nsamples}_samples')
os.makedirs(h5f_folder_edks,exist_ok=True)
os.makedirs(h5f_folder_okada,exist_ok=True)
h5file_dir_edks = os.path.join(h5f_folder_edks,h5f_name_edks)
h5file_dir_okada = os.path.join(h5f_folder_okada,h5f_name_okada)   
#h5file_edks = h5py.File(h5file_dir_edks,'w')
#h5file_okada = h5py.File(h5file_dir_okada,'w')        
#edks_displacement_dset = h5file_edks.create_dataset('displacement',shape=(nsamples,len(xstn)*len(ystn)*3))
#okada_displacement_dset = h5file_okada.create_dataset('displacement',shape=(nsamples,len(xstn)*len(ystn)*3))




patch_interp = patch/factor

xS_interp = np.arange(patch_interp/2 , ncols*patch,patch_interp)
yS_interp = np.arange(-(nrows*patch) + patch_interp/2,0,patch_interp)
XS_interp,YS_interp  = np.meshgrid(xS_interp,yS_interp)

xS = np.arange(patch/2 , ncols*patch,patch)
yS = np.arange(-(nrows-1/2)*patch,0,patch)
    # shift accordingly at surface
yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
XS, YS = np.meshgrid(xS,yS)


DIP = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0)
STRIKE0 = np.flip(df['strike'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0) 
method = 'linear'
dip_interp = RegularGridInterpolator((yS, xS), DIP,bounds_error=False, fill_value=None,method = method)
depth_interp = RegularGridInterpolator((yS, xS), DEPTH,bounds_error=False, fill_value=None,method = method)

dip = dip_interp((YS_interp,XS_interp))
depth = depth_interp((YS_interp,XS_interp))




xsrc = xS_interp*km
ysrc = yS_interp*km
zsrc = depth.flatten()*km
Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
xsrc_flat,ysrc_flat = Xsrc.flatten(),Ysrc.flatten()
width = np.ones_like(xsrc_flat)*patch_interp*km
length = np.ones_like(xsrc_flat)*patch_interp*km
strike = np.ones_like(xsrc_flat)*90
strike_rad = strike*(np.pi/180)
dip_rad= dip.flatten()*(np.pi/180)




#print(xstn_flat.shape,ystn_flat.shape,zstn_flat.shape)
#print(xsrc_flat.shape,ysrc_flat.shape,zsrc.shape,width.shape,dip.shape,depth.shape)
strike0_rad= STRIKE0*(np.pi/180)
dip0_rad= DIP.flatten()*(np.pi/180) 

nu = 0.25


aS = np.ones_like(dip.flatten())*(width[0]*length[0])

ones = np.ones_like(dip.flatten())

# calculate the GF's


rake = 0 * np.ones(aS.shape)
slip = 1 * np.ones(aS.shape)

GFeSS, GFnSS, GFzSS = layered_disloc(xsrc_flat, ysrc_flat, zsrc, strike, dip.flatten(), rake, slip, aS, xstn_flat, ystn_flat, edks, prefix, BIN_EDKS)
rake = 90 * np.ones(aS.shape)
slip = 1 * np.ones(aS.shape)
GFeDS, GFnDS, GFzDS = layered_disloc(xsrc_flat, ysrc_flat, zsrc, strike, dip.flatten(), rake, slip, aS, xstn_flat, ystn_flat, edks, prefix, BIN_EDKS)

file_dir  = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS_GFs/{nsamples}_samples')
os.makedirs(file_dir,exist_ok =True)
file_name = os.path.join(file_dir,'GFeSS.txt')
np.savetxt(file_name,GFeSS)
file_name = os.path.join(file_dir,'GFnSS.txt')
np.savetxt(file_name,GFnSS)
file_name = os.path.join(file_dir,'GFzSS.txt')
np.savetxt(file_name,GFzSS)
file_name = os.path.join(file_dir,'GFeDS.txt')
np.savetxt(file_name,GFeDS)
file_name = os.path.join(file_dir,'GFnDS.txt')
np.savetxt(file_name,GFnDS)
file_name = os.path.join(file_dir,'GFzDS.txt')
np.savetxt(file_name,GFzDS)


