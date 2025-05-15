# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:52:56 2023

@author: joanv
"""

import h5py
import numpy as np
import pandas as pd 

import struct
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Ellipse
from scipy.stats import skew

# Np=132 # Number of subfaults
# file_dir = './Data_Iquique_model/step_056.h5'
# f = h5py.File(file_dir,'r')
# samples = f['Sample Set']

# strike_slip = samples[:,:Np]
# dip_slip = samples[:,Np:2*Np]
# nuisance_p = samples[:,2*Np:2*Np+3]
# rise_time = samples[:,2*Np+3:3*Np+3]
# rupture_velocity = samples[:,3*Np+3:4*Np+3] #must change see 
# hypocenter = samples[:,4*Np+3:]



def stacking(d,keys,n):
    if n!=0:
        x = np.hstack((stacking(d,keys,n-1),d[keys[n+1]]))
    else:
        x = np.hstack((d[keys[n]],d[keys[n+1]]))
    return x
    

def hstacking(d):   
    keys = ['strikeslip','dipslip','risetime','rupturevelocity','hypo_as','hypo_dd']
    n = len(keys)
    return stacking(d,keys,n-2)

def dict_constructor(model_type,sampled_data,Np,nramp=0):
    # mean of each parameter at each node 
    mean = np.mean(sampled_data,axis=1) 

    # standard deviation of each parameter at each node 
    std = np.std(sampled_data,axis=1)
    
    # skewness of each parameter at each node
    
    skewness = skew(sampled_data,axis=1)
    
    if model_type=='static':
            mean_d = {
                    'U_perp':mean[:Np],         # along-strike
                    'U_parallel':mean[Np:2*Np], # along-dip
                    'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2)
                            }
            std_d = {
                   'std_U_perp':std[:Np],         # along-strike
                   'std_U_parallel':std[Np:2*Np], # along-dip
        
                   }
            skew_d = {
                'skew_U_perp':skewness[:Np],
                'skew_U_parallel':skewness[Np:2*Np]
                    }
    else:
        mean_d = {
                'U_perp':mean[:Np],         # along-strike
                'U_parallel':mean[Np:2*Np], # along-dip
                'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2),
                'Tr':mean[2*Np+nramp:3*Np+nramp],
                'Vr':mean[3*Np+nramp:4*Np+nramp],
                'Hypo_as': mean[4*Np+nramp:4*Np+nramp+1], # along-strike
                'Hypo_dd': mean[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                        }
        std_d = {
               'std_U_perp':std[:Np],         # along-strike
               'std_U_parallel':std[Np:2*Np], # along-dip
               'std_Tr':std[2*Np+nramp:3*Np+nramp],
               'std_Vr':std[3*Np+nramp:4*Np+nramp],
               'std_Hypo_as': std[4*Np+nramp:4*Np+nramp+1], # along-strike
               'std_Hypo_dd': std[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                       }
        skew_d = {
            'skew_U_perp':skewness[:Np],
            'skew_U_parallel':skewness[Np:2*Np],
            'skew_Tr':skewness[2*Np+nramp:3*Np+nramp],
            'skew_Vr':skewness[3*Np+nramp:4*Np+nramp],
            'skew_Hypo_as':skewness[4*Np+nramp:4*Np+nramp+1], # along-strike
            'skew_Hypo_dd':skewness[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                }

        return mean_d,std_d,skew_d     


def model_reader(name:str,
               model_type:str,
               file_type:str,
               geom_file_shape:tuple,
               patchsize:int,
               header = [0,1,4,5,6],
               nrows_skip=1,
               model_file_shape=None,
               precision='double',
               sampling=False,
               nsamples=20000):
    

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
    geom_df0 = geom_df.transpose()
    try:
        os.makedirs(geometry_formatted_dir)
    except FileExistsError:
        print(f'Directory {geometry_formatted_dir} already exist')
    
    geom_df0.to_csv(csv_file_dir)
       
    
    nrows = geom_file_shape[0]
    ncols = geom_file_shape[1]
    Np= int(nrows*ncols) # number of patches
    
    EQs_dir = os.path.join(os.getcwd(),'EQ')
    EQ_dir =  os.path.join(EQs_dir,name)
    model_dir = os.path.join(EQ_dir,'model')
    model_type_dir = os.path.join(model_dir,model_type)
    model_file_name = os.listdir(model_type_dir)[0]
    model_file_dir = os.path.join(model_type_dir,model_file_name)
    
    
    output_EQs_dir = os.path.join(os.getcwd(),'output')
    output_EQ_dir =  os.path.join(output_EQs_dir,name)
    output_model_dir = os.path.join(output_EQ_dir,'model')
    output_model_type_dir = os.path.join(output_model_dir,model_type)
    output_mean_model_dir = os.path.join(output_model_type_dir,'mean')
    output_plots_model_dir = os.path.join(output_model_type_dir,'figures')
    
    if file_type=='h5':
        f = h5py.File(model_file_dir,'r')
        try:
            
            data = f['Sample Set']
            data = np.array(data).transpose()
        
        except: 
            data = hstacking(f['ParameterSets']).transpose()
            

    ### Reading Source Model file ###
    elif file_type=='binary':  
        data = np.fromfile(model_file_dir,precision).reshape(model_file_shape)
    

    
    rng = np.random.default_rng()    # uniform sampler
    # uniform sampling `nsamples` from `data`
    
    if sampling==True:
        sampled_data = rng.choice(data,nsamples,axis=1) 
    else:
        sampled_data = data
        
    # mean of each parameter at each node 
    mean = np.mean(sampled_data,axis=1) 

    # standard deviation of each parameter at each node 
    std = np.std(sampled_data,axis=1)
    
    # skewness of each parameter at each node
    
    skewness = skew(sampled_data,axis=1)
    
    # covariance matrix 
    cov_m = np.cov(sampled_data)
    
    
    if name in ('Tohoku','Illapel'):
        if model_type=='static':
            mean_d = {
                    'U_perp':mean[:Np],         # along-strike
                    'U_parallel':mean[Np:2*Np], # along-dip
                    'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2)
                            }
            std_d = {
                   'std_U_perp':std[:Np],         # along-strike
                   'std_U_parallel':std[Np:2*Np], # along-dip
        
                   }
            skew_d = {
                'skew_U_perp':skewness[:Np],
                'skew_U_parallel':skewness[Np:2*Np]
                    }

        else:
            mean_d = {
                    'U_perp':mean[:Np],         # along-strike
                    'U_parallel':mean[Np:2*Np], # along-dip
                    'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2),
                    'Tr':mean[2*Np:3*Np],
                    'Vr':mean[3*Np:4*Np],
                    'Hypo_as': mean[4*Np:4*Np+1], # along-strike
                    'Hypo_dd': mean[4*Np+1:4*Np+2] # down-dip
                            }
            std_d = {
                   'std_U_perp':std[:Np],         # along-strike
                   'std_U_parallel':std[Np:2*Np], # along-dip
                   'std_Tr':std[2*Np:3*Np],
                   'std_Vr':std[3*Np:4*Np],
                   'std_Hypo_as': std[4*Np:4*Np+1], # along-strike
                   'std_Hypo_dd': std[4*Np+1:4*Np+2] # down-dip
                           }
            skew_d = {
                'skew_U_perp':skewness[:Np],
                'skew_U_parallel':skewness[Np:2*Np],
                'skew_Tr':skewness[2*Np:3*Np],
                'skew_Vr':skewness[3*Np:4*Np],
                'skew_Hypo_as':skewness[4*Np:4*Np+1], # along-strike
                'skew_Hypo_dd':skewness[4*Np+1:4*Np+2] # down-dip
                    }
    
    
    elif name=='Iquique':
        nramp = 3
        mean_d = {
                'U_perp':mean[:Np],         # along-strike
                'U_parallel':mean[Np:2*Np], # along-dip
                'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2),
                'Tr':mean[2*Np+nramp:3*Np+nramp],
                'Vr':mean[3*Np+nramp:4*Np+nramp],
                'Hypo_as': mean[4*Np+nramp:4*Np+nramp+1], # along-strike
                'Hypo_dd': mean[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                        }
        std_d = {
                'std_U_perp':std[:Np],
                'std_U_parallel':std[Np:2*Np],
                'std_Tr':std[2*Np+nramp:3*Np+nramp],
                'std_Vr':std[3*Np+nramp:4*Np+nramp],
                'std_Hypo_as':std[4*Np+nramp:4*Np+nramp+1], # along-strike
                'std_Hypo_dd':std[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                    }
        skew_d = {
            'skew_U_perp':skewness[:Np],
            'skew_U_parallel':skewness[Np:2*Np],
            'skew_Tr':skewness[2*Np+3:3*Np],
            'skew_Vr':skewness[3*Np:4*Np],
            'skew_Hypo_as':skewness[4*Np:4*Np+1], # along-strike
            'skew_Hypo_dd':skewness[4*Np+1:4*Np+2] # down-dip
                }
    
    elif name=='Pedernales':
        nramp = 9
        # If RotAngle within ]90, 270], change rotation
        # if RotAngle > 90. and RotAngle<=270.:
        #     rotation += np.pi

        # rp = alt.parameter('strike',samples[:,p]   ,bins=100,dtype=np.float32,color='red')
        # ar = alt.parameter('dip'   ,samples[:,p+Np],bins=100,dtype=np.float32,color='red')

        # ss = copy.deepcopy(rp)
        # ds = copy.deepcopy(ar)
        # ss.values = ar.values*np.cos(rotation) - rp.values*np.sin(rotation)
        # ds.values = ar.values*np.sin(rotation) + rp.values*np.cos(rotation)
        
        mean_d = {
                'U_perp':mean[:Np],         # along-strike
                'U_parallel':mean[Np:2*Np], # along-dip
                'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2),
                'Tr':mean[2*Np+nramp:3*Np+nramp],
                'Vr':mean[3*Np+nramp:4*Np+nramp],
                'Hypo_as': mean[4*Np+nramp:4*Np+nramp+1], # along-strike
                'Hypo_dd': mean[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                        }
        std_d = {
                'std_U_perp':std[:Np],
                'std_U_parallel':std[Np:2*Np],
                'std_Tr':std[2*Np+nramp:3*Np+nramp],
                'std_Vr':std[3*Np+nramp:4*Np+nramp],
                'std_Hypo_as':std[4*Np+nramp:4*Np+nramp+1], # along-strike
                'std_Hypo_dd':std[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                    }
    
    elif name=='Gorkha':
        mean_d = {
                'U_parallel':mean[:Np],         # along-dip
                'U_perp':mean[Np:2*Np], # along-strike
                'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2),
                'Tr':mean[2*Np:3*Np],
                'Vr':mean[3*Np:4*Np]
                        }
        std_d = {
                'std_U_parallel':std[:Np],
                'std_U_perp':std[Np:2*Np],
                'std_Tr':std[2*Np:3*Np],
                'std_Vr':std[3*Np:4*Np]
                    }
    
        
    
    full_m = mean_d | std_d    # full model (mean + std)
    full_m = full_m | geometry_d  # full model (mean +std + geom)
    try: 
        os.makedirs(output_mean_model_dir)    
        print(f'Directory {output_mean_model_dir}  just created')
    except FileExistsError:
        print(f'Directory {output_mean_model_dir}  just created')    
    
    full_model_csv_file_dir = os.path.join(output_mean_model_dir,f'{name}_mean_{model_type}_model.csv')
    full_m_df = pd.DataFrame.from_dict(full_m, orient='index') # handle `Value error exception` 
    full_m_df = full_m_df.transpose()
    full_m_df.to_csv(full_model_csv_file_dir)
    
    def distribution(parameter,color_bar_label):
        
        
        x = np.arange(patchsize/2,ncols*patchsize,patchsize)
        y = np.arange(patchsize/2,nrows*patchsize,patchsize)
        fig, ax = plt.subplots(figsize=(12,10),dpi=600)
        param = mean_d[parameter].reshape(ncols,nrows).transpose()
        X, Y = np.meshgrid(x, y)
        try:
            hypocenter = (mean_d['Hypo_as'],mean_d['Hypo_dd'])
            ax.plot(hypocenter[0],hypocenter[1],marker='*',color='blue',markeredgecolor='black',markersize=25)
        except:
            pass
            #hypocenter = (145,35) # for Gorkha
            #ax.plot(hypocenter[0],hypocenter[1],marker='*',color='blue',markeredgecolor='black',markersize=25)
        im = ax.pcolormesh(x,y,param,edgecolors='k',cmap='YlOrRd')

        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        ax.margins(0)
        if parameter=='Slip':
            U = mean_d['U_parallel'].reshape(ncols,nrows).transpose() # dip-slip displacement
            V = mean_d['U_perp'].reshape(ncols,nrows).transpose()     # strike-slip displacement (see p.4 Minson et al. (2013) Part II)
            q = ax.quiver(X,Y,V,U,scale=0.8,scale_units ='x', units='width',width=0.002,headwidth=4.5,headlength=6)
            ax.quiverkey(q, X=0.1, Y=1.1, U=round(param.max(),-1),
                          label=' {} m'.format(round(param.max(),-1)), labelpos='E')
            factor = 1/q.scale    # scaling factor from UV unit to XY units
            offsetXY = np.stack((V,U),axis=2).flatten().reshape(nrows*ncols,2)   
            
            # minus sign in Y coordinate corrects for y-axis inversion
            
            offsetXY = np.array([[xy[0],-xy[1]] for xy in offsetXY])  
            coord = q.XY + offsetXY*factor
            # width and height had to be formatted in row-like order (initially in column-like order)
            ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
                            width=full_m['std_U_perp'].reshape((nrows,ncols),order='F').flatten()[i]*factor, 
                            height=full_m['std_U_parallel'].reshape((nrows,ncols),order='F').flatten()[i]*factor,
                            angle=0,alpha=0.5,fill=False,edgecolor='k')
                    for i in range(q.N)]
            for e in ells:
                ax.add_artist(e)  
                
        fig.colorbar(im, ax=ax,shrink=0.45,label='{}'.format(color_bar_label))
        im.figure.axes[1].tick_params(axis='both', labelsize=12) 
        
        ax.set_ylabel('Down-dip distance (km)',fontsize=12)
        ax.set_xlabel('Along-strike distance (km)',fontsize=12)
        ax.set_aspect('equal', 'box')
        ax.tick_params(labelsize=12)
        ax.invert_yaxis()
        try: 
            os.makedirs(output_plots_model_dir)  
            print(f'Directory {output_plots_model_dir}  just created')
        except FileExistsError:
            print(f'Directory {output_plots_model_dir} already created')
        plot_path = os.path.join(output_plots_model_dir,f'{name}_mean_{model_type}_{parameter}_{name}.jpg')
        fig.suptitle(f'{name}_{model_type}_model'.replace("_"," "),x=0.18,y=0.9,fontweight='bold')
        fig.savefig(plot_path)

    
    keys = list(mean_d.keys())
    try:
        keys.remove('Hypo_as')
        keys.remove('Hypo_dd')
    except:
        keys = list(mean_d.keys())
        
    for parameter in keys:    
        distribution(parameter,parameter.replace('_',' '))
        

 
    
    