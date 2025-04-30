# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 10:51:35 2023

@author: joanv
"""

import sys
 
# adding csi functions to the system path
sys.path.insert(0, '/home/josevilo/csi/csi')

from okadafull import displacement, stress, strain
import pandas as pd
import numpy as np
import os
import matplolib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def dir_creator(parent_dir,name,model_type,parameter):
    
    output_dir = os.path.join(parent_dir,'ouput')
    EQ_dir = os.path.join(output_dir,name)
    model_dir =  os.path.join(EQ_dir,'model')
    model_type_dir =  os.path.join(model_dir,model_type)
    parameter_dir = os.path.join(model_type_dir ,parameter)
    parameter_data_dir = os.path.join(parameter_dir, 'data')
    parameter_figures_dir = os.path.join(parameter_dir,'figures')
    
    child_dir = {'data':parameter_data_dir,'figures':parameter_figures_dir}
    
    return child_dir

working_dir = os.getcwd()
input_dir = os.path.join(working_dir,'input')
name,model_type = '',''


def dir_finder(parent_dir,name,model_type,sample_type='mean'):
    input_dir = os.path.join(parent_dir,'input')
    EQ_dir = os.path.join(input_dir,name)
    model_dir =  os.path.join(EQ_dir,'model')
    model_type_dir =  os.path.join(model_dir,model_type)
    sample_type_dir = os.path.join(model_type_dir,sample_type)
    input_file_dir = os.path.join(sample_type_dir,f'{name}_{sample_type}_{model_type}_model.csv')
    
    return input_file_dir

def okada_across_fault_calculator(name,model_type,patchsize,shape_geom_file,sample_type='mean',custom_strike=False,strike_value='Default',mu=30.0*1e9,nu=0.25):
    working_dir = os.getcwd()
    input_file_dir = dir_finder(working_dir,name,model_type)    
     
    # Accesing data with mean model ; need to include sampled model case
    df = pd.read_csv(input_file_dir)
    df= df.drop(df.columns[0],axis=1)
    
    nrows = shape_geom_file[0] # before for Tohoku = 9 , Iquique = 11
    ncols = shape_geom_file[1] # before for Tohoku = 24, Iquique = 12
    Uparallel = df['U_parallel'].values.reshape(ncols,nrows).transpose().flatten()
    Uperp = df['U_perp'].values.reshape(ncols,nrows).transpose().flatten()
    slip =  df['Slip'].values.reshape(ncols,nrows).transpose().flatten()
    dip = df['dip'].values.reshape(ncols,nrows).transpose().flatten()
    depth = df['depth'].values.reshape(ncols,nrows).transpose().flatten()
    strike = df['strike'].values.reshape(ncols,nrows).transpose().flatten()
    
    try:
        hypocenter = [df['Hypo_as'][0],df['Hypo_dd'][0]]
    except:
        hypocenter = [10,10]
    #input_file_dir = os.path.join(input_data_directory,'n_10000_sampled_models_kinematic_tohoku.dat')
    #data = np.fromfile(input_file_dir,'float').reshape((866,10000))
    
    #-----------------------------------
    # 1. Dislocation and receiver
    #-----------------------------------
    
    patchsize  = patchsize*1e3  # Tohoku =29, Iquique = 17
    
    
    # DISLOCATION
    x = np.arange(patchsize/2,ncols*patchsize,patchsize)
    y = np.arange(patchsize/2,nrows*patchsize,patchsize)
    x,y = np.meshgrid(x,y)
    x, y = x.flatten(),y.flatten()
    xc = x
    yc = y
    zc  = depth*1e3
    
    
    # STATION
    
    xs = xc
    ys = yc
    zs = -zc # before it was at zero dept h=  np.zeros(xc.shape)
    
    
    length_c = patchsize*np.ones(xc.shape)
    width_c = patchsize*np.ones(yc.shape)
    dip_c = dip*(np.pi)/180
    
    if custom_strike:
        strike_c = strike_value*(np.pi)/180*np.ones(xc.shape)
    else:    
        strike_c = strike*(np.pi)/180 
    
        
    # DISLOCATION 
    ts = np.zeros(xc.shape)
    ss = Uperp
    ds = Uparallel
    
    
    # Calling Okada module for parameter calculations 
    DISPLACEMENT= displacement(xs, ys,zs, xc, yc, zc, width_c, length_c, strike_c, dip_c, ss, ds, ts, nu=0.25)
    STRAIN = strain(xs, ys,zs, xc, yc, zc, width_c, length_c, strike_c, dip_c, ss, ds, ts, nu=0.25)
    STRESS = stress(xs, ys,zs, xc, yc, zc, width_c, length_c, strike_c, dip_c, ss, ds, ts, nu=0.25)
    
    
    # Directory placeholders
    
    # local function to save csv files
    def output_file_writer(parameter,arr):
        parameter_dirs = dir_creator(working_dir,name,model_type,parameter)
        parameter_file = os.path.join(parameter_dirs['data'],f'{name}_{sample_type}_{model_type}_model_{parameter}_strike_{strike_value}.txt')
 
        os.makedirs(parameter_dirs['data'],exist_ok=True)
        np.savetxt(parameter_file,arr)
        
    def output_figure(parameter='stress'):
        parameter_dirs = dir_creator(working_dir,name,model_type,parameter)
        parameter_file = os.path.join(parameter_dirs['figures'],f'{name}_{sample_type}_{model_type}_model_{parameter}_strike_{strike_value}.jpg')
        os.makedirs(parameter_dirs['figures'],exist_ok=True)
        
        return parameter_file
        
    output_file_writer('displacement',DISPLACEMENT)    
    output_file_writer('stress',STRESS[0])
    output_file_writer('strain',STRAIN)
    
    
    disp_pd = pd.DataFrame(DISPLACEMENT,columns=['x','y','z'])
    stress_pd = pd.DataFrame(STRESS[0],columns=['xx','xy','xz','yy','yz','zz'])
   
    displacement_dict = {'displacement':disp_pd}
    stress_dict = {'stress':stress_pd}
    
    # set fault geometry
    
    
    # create placeholder matrices
    
    # Create tnn, ts1, ts2 
    npatches = int(shape_geom_file[0]*shape_geom_file[1])
    tnn = np.zeros(npatches)
    ts1 = np.zeros(npatches)
    ts2 = np.zeros(npatches)
    # compute stress for each patch (total number: ns)
    for k in range(npatches):
        
        Sxx = stress_dict['stress']['xx'].loc[k]*1e-6
        Sxy = stress_dict['stress']['xy'].loc[k]*1e-6
        Sxz = stress_dict['stress']['xz'].loc[k]*1e-6
        Syy = stress_dict['stress']['yy'].loc[k]*1e-6
        Syz = stress_dict['stress']['yz'].loc[k]*1e-6
        Szz = stress_dict['stress']['zz'].loc[k]*1e-6
        
        # assemble stress tensor
        
        tau = [[Sxx,Sxy,Sxz],[Sxy,Syy,Syz],[Sxz,Syz,Szz]]
        
        #tau = np.array([[S[k,0], S[k,1], S[k,2]],[S[k,1], S[k,3], S[k,4]], [S[k,2], S[k,4], S[k,5]]])
        
        # set cosine/sine constants
        CS, SS = np.cos(strike), np.sin(strike)
        CD, SD = np.cos(dip), np.sin(dip)
    
        # set direction vectors
        NVec = [SD*CS, -SS*SD, CD] # fault-normal
        Ss = [SS, -CS, 0] # along-strike
        Sd = [-CD*CS, CS*SS, SD] # along-dip
        
        # compute traction vector
        trn = np.matmul(tau, NVec) 
        tnn[k] = np.dot(trn, NVec) # fault-normal
        ts1[k] = np.dot(trn, Ss) # along-strike 
        ts2[k] = np.dot(trn, Sd) # along-dip
        
      
        
   
    fig, axes = plt.subplots(3,1,figsize=(12,10),dpi=600)
    traction_d = {'Normal':tnn,'Along-strike Shear':ts1,'Along-dip Shear':ts2}
    
    X, Y = np.meshgrid(x, y)
    for i,parameter_id in enumerate(traction_d.keys()):
        
        parameter = traction_d[parameter_id]
        title = f"{parameter_id} Stress across fault $\phi={strike_value}$"
        # parameter Stress-related is already in correct order ie. row-wise as obtained from Okada
        im = axes[i].pcolormesh(x,y,parameter.reshape(ncols,nrows),edgecolors='k',
                          cmap='bwr',
                          norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
        
        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        axes[i].margins(0)
                
        fig.colorbar(im, ax=axes[i],shrink=0.9,label='{}'.format('Stress (MPa)'))
        try:
            
            axes[i].plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markersize=18,markeredgecolor='black')
        except:
            pass
        
        CS = axes[i].contour(X, Y, slip, [10],linewidths=3)

        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
    
        #U = ts1.reshape(nrows,ncols)     # along-strike 
        #V = ts2.reshape(nrows,ncols)     # along-dip 
        # need to invert (U,V)=(Along-strike stress,Along-dip stress) as well by setting `angles`=`xy` 
        #q = ax.quiver(X,Y,U,V,angles='xy',scale=1.00,scale_units ='x', units='width',width=0.002,headwidth=4.5,headlength=6)
        #ax.quiverkey(q, X=0.1, Y=1.1, U=50,
         #             label='50MPa', labelpos='E')
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=12)
        axes[i].invert_yaxis()
    fig.suptitle(f'{name}_{model_type}_model'.replace("_"," "),x=0.18,y=0.9,fontweight='bold')
    plt.tight_layout()
    figure_dir = output_figure()
    fig.savefig(figure_dir)



okada_calculator('Tohoku','kinematic',29, (9,24))
okada_calculator('Tohoku','static',29, (9,24))
okada_calculator('Iquique','kinematic',17, (11,12))
okada_calculator('Illapel','kinematic',18, (10,17))
okada_calculator('Gorkha','kinematic',10,(9,18))

