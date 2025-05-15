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
import time
# EDKS
from layered_disloc import layered_disloc


# Check inputs
def ArraySizes(*args):
    '''
    Only requirement is that each arguments has the same size and can be converted to a numpy array
    Returns : Numpy arrays
    '''
    
    # Create a list of sizes
    Sizes = []
    Arrays = []

    # Check class
    for arg in args:
        if arg.__class__ in (list, tuple):
            arg = np.array(arg)
        elif arg.__class__ in (float, np.float64, int):
            arg = np.array([arg])
        Arrays.append(arg)
        Sizes.append(arg.shape)
    
    # Assert sizes
    assert (len(np.unique(Sizes))==1), 'The {} provided arrays are not the same size'.format(len(args))

    # All done
    return Arrays

#--------------------------------------------------
# Displacements only
def displacement(xs, ys, zs, xc, yc, zc, width, length, strike, dip, ss, ds, ts, nu=0.25):
    '''
    Returns the displacements at the stations located on (xs, ys, zs) for patches
        with centers on (xc, yc, zc). All arguments can be float, list or array.
    '''

    # Here Mu can be anything. RJ tested it and the displacement is not-sensitive to Mu as it should be.
    # Although, it does not work with Mu = 0.0 GPa... So we take a random value of 30GPa
    mu = 30e9

    # Nu does matter here, and it is by default 0.25

    # Check 
    xs, ys, zs = ArraySizes(xs, ys, zs)
    xc, yc, zc, width, length, strike, dip, ss, ds, ts = ArraySizes(xc, yc, zc, width, length, strike, dip, ss, ds, ts)

    # Normally, StaticInv does angles in Radians
    dip = dip*180./np.pi
    strike = strike*180./np.pi

    # Run okada
    u, d, s, flag, flag2 = ok92.okada92(xs, ys, zs, xc, yc, zc, length, width, dip, strike, ss, ds, ts, mu, nu)

    # Check if things went well
    if not (flag==0).all():
        if not np.where(flag!=0)==[]:
            print(' Error: {}'.format(tuple(np.where(flag!=0))))
            print('Something went wrong in okada4py... You should check...')

    # Reshape the displacement
    u = u.reshape((len(xs), 3))

    # All Done
    return u
        


        

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
BIN_EDKS = '${HOME}/static/edks/bin'
input_dir  = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{method}/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')

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


file_folder = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{method}/{nsamples}_samples/bin_data')
file_dir = os.path.join(file_folder,f'{name}_{model_type}_n_{nsamples}.dat')
bin_data = np.fromfile(file_dir,'float').reshape((nparam,nsamples))
	
h5f_name_edks =  f'EDKS_{name}_displacement_nsamples_{nsamples}.h5'
h5f_name_okada =  f'Okada_{name}_displacement_nsamples_{nsamples}.h5'
h5f_folder_edks = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS/{nsamples}_samples')
h5f_folder_okada = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/Okada/{nsamples}_samples')
os.makedirs(h5f_folder_edks,exist_ok=True)
os.makedirs(h5f_folder_okada,exist_ok=True)
h5file_dir_edks = os.path.join(h5f_folder_edks,h5f_name_edks)
h5file_dir_okada = os.path.join(h5f_folder_okada,h5f_name_okada)   
h5file_edks = h5py.File(h5file_dir_edks,'w')
h5file_okada = h5py.File(h5file_dir_okada,'w')        
edks_displacement_dset = h5file_edks.create_dataset('displacement',shape=(nsamples,len(xstn)*len(ystn)*3))
okada_displacement_dset = h5file_okada.create_dataset('displacement',shape=(nsamples,len(xstn)*len(ystn)*3))




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

nu = 0.3400051746245543

start = time.time()
for i in range(bin_data.shape[1]):
#for i in range(5):    
    Np = nrows*ncols
    data = bin_data[:,i] 
    data0 = np.copy(data)
    if name =='Gorkha':
          data[:Np],data[Np:2*Np] = data0[Np:2*Np],data0[:Np]
    if name == 'Pedernales':
        ss = np.zeros(Np)
        ds = np.zeros(Np)
        RotAngle = 360-99
        for p in range(Np):
            RotAngle2 =RotAngle*((np.pi) / 180.)
            
            rotation = np.arctan2(np.tan(strike0_rad[p]) - np.tan(RotAngle2), 
                                np.cos(dip0_rad[p])*(1.+np.tan(RotAngle2)*np.tan(strike0_rad[p])))

            # If RotAngle within ]90, 270], change rotation
            if RotAngle > 90. and RotAngle<=270.:
                rotation += np.pi

            rp = data[p]    # rake-perpendicular
            ar = data[p+Np] # rake-parallel

            ss[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
            ds[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
        UPERP = ss
        UPARALLEL = ds    
    else:
        rake = models[name]['rake']
        rake0 = (rake - 90)*(np.pi/180)
        UPERP =  data[:Np]*np.cos(rake0)  - data[Np:2*Np]*np.sin(rake0) 
        UPARALLEL = data[:Np]*np.sin(rake0)  + data[Np:2*Np]*np.cos(rake0)   

    
    UPERP = np.flip(UPERP.reshape(nrows,ncols,order='F'),axis=0)
   
    UPARALLEL = np.flip(UPARALLEL.reshape(nrows,ncols,order='F'),axis=0)
    
    Uperp_interp = RegularGridInterpolator((yS, xS), UPERP,bounds_error=False, fill_value=None,method = method)
    Uparallel_interp = RegularGridInterpolator((yS, xS), UPARALLEL,bounds_error=False, fill_value=None,method = method)
    
    Uperp = Uperp_interp((YS_interp,XS_interp))
    Uparallel = Uparallel_interp((YS_interp,XS_interp))
    #Slip_interp = np.sqrt(Uperp**2 + Uparallel**2)
    #print(Uperp.shape,Uparallel.shape)		


    
    
    ss = Uperp.flatten()
    ds = Uparallel.flatten()
    ts = np.zeros(ss.shape)
    Displacement = displacement(
                        xstn_flat,
                        ystn_flat,
                        zstn_flat, 
                        xsrc_flat, 
                        ysrc_flat, 
                        zsrc, 
                        width, 
                        length, 
                        strike_rad, 
                        dip_rad, 
                        ss, 
                        ds, 
                        ts, nu=nu)
    
    E_okada = Displacement[:,0]
    N_okada = Displacement[:,1]
    Z_okada = Displacement[:,2]
    
    okada_displacement_dset[i,:] = Displacement.flatten(order='F')
 
     
   
        
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
    
    ss = Uperp.flatten()
    ds = Uparallel.flatten()
    ts = np.zeros_like(Uperp.flatten())
    
    nu = 0.3400051746245543
    
    
    
    
    
    
    # EDKS initialization
    
    # reshape to form vectors.
    # source coordinates
    
    aS = np.ones_like(dip.flatten())*(width[0]*length[0])
    
    SSlip = ss
    DSlip = ds
    
    # calculate the GF's
    
    if i==0:
        rake = 0 * np.ones(SSlip.shape)
        slip = 1 * np.ones(SSlip.shape)
        GFeSS, GFnSS, GFzSS = layered_disloc(xsrc_flat, ysrc_flat, zsrc, strike, dip.flatten(), rake, slip, aS, xstn_flat, ystn_flat,\
                                  edks, prefix, BIN_EDKS)
        rake = 90* np.ones(SSlip.shape)
        slip = 1 * np.ones(SSlip.shape)
        GFeDS, GFnDS, GFzDS = layered_disloc(xsrc_flat, ysrc_flat, zsrc, strike, dip.flatten(), rake, slip, aS, xstn_flat, ystn_flat,\
                                  edks, prefix, BIN_EDKS)
  
     # compute forward displacement calculation.
    E_edks = np.dot(GFeSS,SSlip) + np.dot(GFeDS, DSlip)
    N_edks = np.dot(GFnSS,SSlip) + np.dot(GFnDS, DSlip)
    Z_edks = np.dot(GFzSS,SSlip) + np.dot(GFzDS, DSlip)
    edks_displacement_dset[i,:] = np.array([E_edks,N_edks,Z_edks]).flatten()
    print(f'model: {i} done!')
    
h5file_okada.close()
h5file_edks.close()
end = time.time()
print(end - start)
    
     
        
 
    
 
    
 
    




dict_okada = {'x':E_okada,'y':N_okada,'z':Z_okada}
df_okada = pd.DataFrame(dict_okada)

dict_edks = {'x':E_edks,'y':N_edks,'z':Z_edks}
df_edks = pd.DataFrame(dict_edks)

delta_E =  np.abs(E_okada - E_edks)
delta_N =   np.abs(N_okada - N_edks)
delta_Z =   np.abs(Z_okada - Z_edks)
dict_residual = {'x':delta_E ,'y':delta_N,'z':delta_Z} 
df_residual = pd.DataFrame(dict_residual)

dfs = [df_okada,df_edks,df_residual]

shape_stn = (len(ystn),len(xstn))
cols = ['Okada','EDKS','Residual']
fig, axes = plt.subplots(3,3,figsize=(14,12),dpi=600)
for j, df in enumerate(dfs):
       
       for i,parameter_id in enumerate(df.keys()):
           test = 'with interpolation (Okada halfspace and EDKS layered)'
           parameter = df[parameter_id].values
           title = f"Test {test}"
           file_name = f"Test {test}"
           # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
           if j!=2:
               im = axes[i][j].pcolormesh(xstn,ystn,parameter.reshape(shape_stn), cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
               fig.colorbar(im, ax=axes[i][j],shrink=0.18,label='{}'.format(f'{parameter_id} (m)'))
           else:
               im = axes[i][j].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),cmap = 'viridis')
               fig.colorbar(im, ax=axes[i][j],shrink=0.18,label='{}'.format('abs. residual (m)'))
#ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
           
                   
           #fig.colorbar(im, ax=axes[i][j],shrink=0.6,label='{}'.format(f'{parameter_id} (m)'))
       
       
               
           axes[i][j].set_ylabel('Trench-normal distance (km)',fontsize=7)
           axes[i][j].set_xlabel('Along-strike distance (km)',fontsize=7)
           axes[i][j].set_aspect('equal', 'box')
           axes[i][j].tick_params(labelsize=7)

           axes[i][j].set_title(f'{cols[j]}' ,fontweight='bold',fontsize=11)
           plt.subplots_adjust(wspace=0.5,hspace=-0.6)

fig.suptitle(title,x=0.5,y=0.85,fontweight='bold') # uncomment later
#plt.tight_layout()
#plt.subplots_adjust(wspace=0,hspace=0.75)

fig_dir  = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS_Okada_plot/{nsamples}_samples')
os.makedirs(fig_dir,exist_ok =True)
fig_name = os.path.join(fig_dir,file_name.replace(" ","_") + '.png')
fig.savefig(fig_name)
plt.close()
