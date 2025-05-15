
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 13:22:10 2024

@author: joanv
"""



import re
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
#import h5py as h5
import io
import os
import sys
from ruptureTime import FastSweep, Hypocenter, Grid
# working_dir = f'/home/josevilo/dynamic/MudPy/'
working_dir = os.getcwd()

name = str(sys.argv[1])
nmodels = int(sys.argv[2])

factor_interp = int(sys.argv[3])

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


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,90,90,0,107]
z_offsets = [2.8,0,0,0,0]
nramps = [0,3,0,9,0]
def build_model(names,geoms,patches,arrow_sizes,nparams,rakes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['z_offset'] = z_offsets[i]
        model[name]['nramp'] = nramps[i]
    return model


def get_mu(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VP'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        return mu

models = build_model(names,geoms,patches,arrow_sizes,nparams,rakes)

km = 1000


model_type = 'kinematic'
nsamples = 100
z_offset = models[name]['z_offset']

earth_model_dir = os.path.join(working_dir,f'EQ/{name}/earth_model/model_{name}')
f = open(earth_model_dir, "r")
content = f.readlines()

#nlayers =int(sys.argv[4])

nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
nparam = models[name]['nparam']

shape = (nrows,ncols)
patch = models[name]['patch']
nramp = models[name]['nramp']
earth_model = content[12:]
keys = re.findall('\w+',content[11],flags=re.DOTALL)
fmt = re.compile('\d+.\d+')
layers = [list(map(float,re.findall(fmt,m))) for m in earth_model]
#print(layers)
layers = np.array(layers).T
print(layers)
layers0 = np.zeros(layers.shape)
layers0[0,:] = layers[0,:]
layers0[1,:] = layers[2,:]
layers0[2,:] = layers[1,:]
layers0[3,:] = layers[3,:]
layers0[4,:] = layers[5,:]
layers0[5,:] = layers[4,:]
model_dict = dict(zip(keys,layers0))
#print(layers0)

model_df = pd.DataFrame(model_dict)
print(model_df)



    
    #logs_folder = os.path.join(working_dir,f'{name}/{nsamples}_samples/{j}/logs/')
    #os.makedirs(logs_folder,exist_ok=True)
for j in range(1,nmodels+1):
    GF_folder = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j}/GFs/')
    os.makedirs(GF_folder,exist_ok=True)
    dynamic_folder = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j}/GFs/dynamic')
    matrices_folder = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j}/GFs/matrices')
    
    os.makedirs(dynamic_folder,exist_ok=True)
    os.makedirs(matrices_folder,exist_ok=True)

j=1
mod_folder = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j}/structure/')
os.makedirs(mod_folder,exist_ok=True)
mod_file_dir = os.path.join(mod_folder,f'{name}.mod')

f = io.open(mod_file_dir,'w',newline='\n')
for i in range(layers.shape[1]):
    line = list(model_df.iloc[i].values)
    string = ''
    for k,item in enumerate(line):
        if k==0:
            string = string + str(item)
        else:
            string = string + '\t'+ str(item)
    string+='\n'
    
    f.write(string)
f.close()

file_folder = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/bin_data')
file_dir = os.path.join(file_folder,f'{name}_{model_type}_n_{nsamples}.dat')
bin_data = np.fromfile(file_dir,'float').reshape((nparam,nsamples))



mean_data_dir = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_kinematic_model.csv')
df = pd.read_csv(mean_data_dir)



LON = one_array_formatter(df['lon'].values,shape)
LAT = one_array_formatter(df['lat'].values,shape)

j = 1    
header = '#station	lon    lat	"static,disp,vel,tsun,insar"									"Static file, displacement file, velocity file,tsunami file,strain file"															"static sigmas(n,e,u),displacement sigmas(n,e,u),velocity sigmas(n,e,u),tsunami sigma,strain sigmas(5 components?)"\n'
common_text ='0	1	0	1	0	/foo/bar	/foo/bar	/Users/dmelgar/Slip_inv/Nepal_forward_test/data/waveforms/3086.vel		/foo/bar	/foo/bar	1	1	1	1	1	1	1	1	1	0	0	0	0	1	0' 

stn_folder = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j}/data/station_info')
os.makedirs(stn_folder,exist_ok=True)
stn0_file_dir = os.path.join(stn_folder,f'{name}.gflist')
f = open(stn0_file_dir, 'w')
f.write(header)
for i in  range(len(LON)):
    f.write('%.4d'%(i)+'\t'+'%s'%(LON[i])+'\t'+'%s'%(LAT[i])+'\t'+common_text +'\n')
f.close()


header = '#station	lon    lat	\n'
stn1_file_dir = os.path.join(stn_folder,f'{name}.sta')
f = open(stn1_file_dir, 'w')
f.write(header)
for i in  range(len(LON)):
    f.write('%.4d'%(i)+'\t'+'%s'%(LON[i])+'\t'+'%s'%(LAT[i])+'\n')
f.close()





# df.insert(8, "onset", onset, True)
# df.insert(5,"type",rise_type,True)

# L = np.ones_like(onset)*patchsize*1e3
# W = np.ones_like(onset)*patchsize*1e3

# df.insert(10, "L(m)", L, True)
# df.insert(11,"W(m)",W,True)


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



patch_interp = patch/factor_interp

xS_interp = np.arange(patch_interp/2 , ncols*patch,patch_interp)
yS_interp = np.arange(-(nrows*patch) + patch_interp/2,0,patch_interp)
XS_interp,YS_interp  = np.meshgrid(xS_interp,yS_interp)

xS = np.arange(patch/2 , ncols*patch,patch)
yS_unproj = np.arange(-(nrows-1/2)*patch,0,patch)
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

DIP = dip_interp((YS_interp,XS_interp))
DEPTH = depth_interp((YS_interp,XS_interp))
LON = lon_interp((YS_interp,XS_interp))
LAT = lat_interp((YS_interp,XS_interp))
STRIKE = strike_interp((YS_interp,XS_interp))

DIP  = DIP.flatten()
STRIKE = STRIKE.flatten()
DEPTH = DEPTH.flatten()
LON = LON.flatten()
LAT = LAT.flatten()


rise_type = np.ones_like(STRIKE)*0.5

xsrc = xS_interp*km
ysrc = yS_interp*km
zsrc = DEPTH*km
Xsrc,Ysrc = np.meshgrid(xsrc,ysrc)
xsrc_flat,ysrc_flat = Xsrc.flatten(),Ysrc.flatten()
W = np.ones_like(xsrc_flat)*patch_interp*km
L = np.ones_like(xsrc_flat)*patch_interp*km
strike = np.ones_like(xsrc_flat)*90
strike_rad = strike*(np.pi/180)
dip_rad= DIP*(np.pi/180)

strike0_rad=  np.flip(df['strike'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()*(np.pi/180)
dip0_rad= DIP.flatten()*(np.pi/180) 
M0_total = 0
MU = np.zeros_like(STRIKE)
for i in range(len(MU)):
    mu = get_mu(dict(model_df),DEPTH[i])
    MU[i] = mu
    # A = (patch_interp*1e3)**2
    # M0 = A*mu*SLIP[i]
    # # print(f'depth = {round(DEPTH[i],3)}, mu = {round(mu*1e-9,3)}, M0 = {round(M0,3)}')
    # M0_total += M0
    

#for i in range(bin_data.shape[1]):
for j in range(0,nmodels):    
    Np = nrows*ncols
    model = bin_data[:,j] 
    model0 = np.copy(model)
    if name =='Gorkha':
          model[:Np],model[Np:2*Np] = model0[Np:2*Np],model0[:Np]
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

            rp = model[p]    # rake-perpendicular
            ar = model[p+Np] # rake-parallel

            ss[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
            ds[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
        UPERP = ss
        UPARALLEL = ds    
    else:
        rake = models[name]['rake']
        rake0 = (rake - 90)*(np.pi/180)
        UPERP =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0) 
        UPARALLEL = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)   
 

    #rake = models[name]['rake']
    #rake0 = (rake - 90)*(np.pi/180)
    
    #UPERP =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0)
    UPERP = two_array_formatter(UPERP,shape)
    
    #UPARALLEL = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)
    UPARALLEL = two_array_formatter(UPARALLEL,shape)
    
    
    Uperp_interp = RegularGridInterpolator((yS, xS), UPERP,bounds_error=False, fill_value=None,method = method)
    Uparallel_interp = RegularGridInterpolator((yS, xS), UPARALLEL,bounds_error=False, fill_value=None,method = method)
    
    Uperp = Uperp_interp((YS_interp,XS_interp))
    Uparallel = Uparallel_interp((YS_interp,XS_interp))
    
    Uperp = Uperp.flatten()
    Uparallel = Uparallel.flatten()
    
       
    DURATION = model[2*Np+nramp:3*Np+nramp]
    DURATION = two_array_formatter(DURATION,shape)
    Duration_interp = RegularGridInterpolator((yS, xS), DURATION,bounds_error=False, fill_value=None,method = method)
    
    Duration = Duration_interp((YS_interp,XS_interp))
    DURATION = Duration.flatten()
    
    
    
    VR = model[3*Np+nramp:4*Np+nramp]
    VR  = two_array_formatter(VR,shape)
    #Vr_interp = RegularGridInterpolator((yS, xS), VR,bounds_error=False, fill_value=None,method = method)
    #VR = Vr_interp((YS_interp,XS_interp)) # no need to flatten
    
    
    
    hypo_as = float(model[4*Np+nramp:4*Np+nramp+1])
    hypo_ad = -float(model[4*Np+nramp+1:4*Np+nramp+2]) # must be negative in my convention
    
    
    solver = FastSweep()
    solver.setGrid(xS,yS_unproj,VR) # approximation
    solver.setHypo(hypo_as,hypo_ad)
    solver.fastSweep(verbose=False)
    t0 = solver.t0[1:-1,1:-1].flatten()
    T0 = t0.reshape(shape)
    t0_interp = RegularGridInterpolator((yS, xS), T0,bounds_error=False, fill_value=None,method = method)
    t0 = t0_interp((YS_interp,XS_interp))
    t0 = t0.flatten()
    
    ONSET = t0
    
        
    N = [int(n)+1 for n in range(len(Uperp))]
    fault_dict = {'No':N,
                  'Lon':LON,
                  'Lat':LAT,
                  'depth(km)':DEPTH -z_offset,
                  'strike':STRIKE,
                  'dip':DIP,
                  'rise_type':rise_type,
                  'rise_time(s)':DURATION,
                  'length(m)':L,
                  'Width(m)':W
                  }
    
    ## halfspace ##
    rupt_dict = {'No.':N,
                  'lon':LON,
                  'lat':LAT,
                  'z(km)':DEPTH - z_offset,
                  'str':STRIKE,
                  'dip':DIP,
                  'type':rise_type,
                  'rise(s)':DURATION,
                  'ss(m)':Uperp,
                  'ds(m)':Uparallel,
                  'l(m)':L,
                  'W(m)':W,
                  'onset(s)':ONSET,
                  'mu(Pa)':MU
                  }
    
    
        
    
    
    df_fault = pd.DataFrame(fault_dict)
    
    fault_file_dir = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j+1}/data/model_info/{name}.fault')
    os.makedirs(os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j+1}/data/model_info/'),exist_ok=True)
    f = io.open(fault_file_dir, 'w', newline='\n')
    init_l = '# No,Lon, Lat,depth(km),strike,dip,rise_type,rise_time(s),length(m),width(m)\n'
    f.write(init_l)
    for i in range(len(ONSET)):
        line = list(df_fault.iloc[i].values)
        string = ''
        for k,item in enumerate(line):
            if k==0:
                string = string + str(item)
            else:
                string = string + ' '+ str(item)
        string+='\n'
        
        f.write(string)
            
    f.close()
    
    
    df_rupt = pd.DataFrame(rupt_dict)
    rupt_file_dir = os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j+1}/forward_models/{name}.rupt')
    os.makedirs(os.path.join(working_dir,f'Dynamic_Simulations/{name}/{nsamples}_samples/{j+1}/forward_models/'),exist_ok=True)
    
    f = io.open(rupt_file_dir, 'w', newline='\n')
    init_l = '#No.  lon     lat   z(km)   str dip type rise(s) ss(m) ds(m) L(m)    W(m) onset(s)    mu(Pa)\n'
    f.write(init_l)
    for i in range(len(ONSET)):
        line = list(df_rupt.iloc[i].values)
        string = ''
        for k,item in enumerate(line):
            if k==0:
                string = string + str(item)
            else:
                string = string + ' '+ str(item)
        string+='\n'
        
        f.write(string)
            
    f.close()
    
    

