# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 16:28:55 2023

@author: joanv
"""


import random
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
import os 
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

def okada_calculator(name,model_type,patchsize,shape_geom_file,sample_type='mean',custom_strike=False,strike_value='Default',mu=30.0*1e9,nu=0.25):
    working_dir = os.getcwd()
    input_file_dir = dir_finder(working_dir,name,model_type)    
     
    # Accesing data with mean model ; need to include sampled model case
    df = pd.read_csv(input_file_dir)
    df= df.drop(df.columns[0],axis=1)
    
    nrows = shape_geom_file[0] # before for Tohoku = 9 , Iquique = 11
    ncols = shape_geom_file[1] # before for Tohoku = 24, Iquique = 12
    Uparallel = df['U_parallel'].values.reshape(ncols,nrows).transpose().flatten()
    Uperp = df['U_perp'].values.reshape(ncols,nrows).transpose().flatten()
    slip = 
    dip = df['dip'].values.reshape(ncols,nrows).transpose().flatten()
    depth = df['depth'].values.reshape(ncols,nrows).transpose().flatten()
    strike = df['strike'].values.reshape(ncols,nrows).transpose().flatten()
    
    parameter_dirs = dir_creator(working_dir,name,model_type,parameter)
    

# relevant directories
def normal_shear_stress_map(data_folder_name,model_type,strike_arr,dip_arr,strike_val,nrows,ncols,patch_size,hypocenter):

    working_dir = os.getcwd()
    okada_files_dir = os.path.join(working_dir,'Data_from_Okada')
    case_files_dir = os.path.join(okada_files_dir,data_folder_name) # change input
    disp_files_dir = os.path.join(case_files_dir,'displacement')       
    stress_files_dir = os.path.join(case_files_dir,'stress')
    
    
    # list of files 
    disp_files = os.listdir(disp_files_dir)
    stress_files = os.listdir(stress_files_dir)
    
    # displacement: mean and std 
    disp_list = [np.loadtxt(os.path.join(disp_files_dir,d_f),delimiter=' ') for d_f in disp_files]
    disp_mean = np.mean(np.array(disp_list),axis=0)
    disp_std =  np.std(np.array(disp_list),axis=0)
    
    # stress: mean and std 
    stress_list = [np.loadtxt(os.path.join(stress_files_dir,sts_f),delimiter=' ') for sts_f in stress_files]
    stress_mean = np.mean(np.array(stress_list),axis=0)
    stress_std = np.std(np.array(stress_list),axis=0)
    
    # transform all to datafram for convenience 
    disp_pd = pd.DataFrame(disp_mean,columns=['x','y','z'])
    stress_pd = pd.DataFrame(stress_mean,columns=['xx','xy','xz','yy','yz','zz'])
    std_disp_pd = pd.DataFrame(disp_std,columns=['x','y','z'])
    std_stress_pd = pd.DataFrame(stress_std,columns=['xx','xy','xz','yy','yz','zz'])
    
    displacement_dict = {
                        'Mean':disp_pd,
                         'Std':std_disp_pd
                         }
    stress_dict = {
                    'Mean':stress_pd,
                    'Std':std_stress_pd
                         }
    
    # set fault geometry
    
    
    # create placeholder matrices
    
    # Create tnn, ts1, ts2 
    ns = len(stress_dict['Mean'].values)
    tnn = np.zeros(ns)
    ts1 = np.zeros(ns)
    ts2 = np.zeros(ns)
    # compute stress for each patch (total number: ns)
    for k in range(ns):
        
        Sxx = stress_dict['Mean']['xx'].loc[k]*1e-6
        Sxy = stress_dict['Mean']['xy'].loc[k]*1e-6
        Sxz = stress_dict['Mean']['xz'].loc[k]*1e-6
        Syy = stress_dict['Mean']['yy'].loc[k]*1e-6
        Syz = stress_dict['Mean']['yz'].loc[k]*1e-6
        Szz = stress_dict['Mean']['zz'].loc[k]*1e-6
        
        # assemble stress tensor
        
        tau = [[Sxx,Sxy,Sxz],[Sxy,Syy,Syz],[Sxz,Syz,Szz]]
        
        #tau = np.array([[S[k,0], S[k,1], S[k,2]],[S[k,1], S[k,3], S[k,4]], [S[k,2], S[k,4], S[k,5]]])
        
        strike = strike_arr[k]*(np.pi/180) # strike angle: 0 deg
        dip = dip_arr[k]*(np.pi/180); # dip angle: 20 deg
    
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
        
      
        
   
    filename = '{}_Stress_over_fault_{}_strike'.format(model_type,strike_val)
    fig, axes = plt.subplots(3,1,figsize=(12,10),dpi=600)
    traction_d = {'Normal':tnn,'Along-strike Shear':ts1,'Along-dip Shear':ts2}
    
    x = np.arange(patch_size/2,ncols*patch_size,patch_size)
    y = np.arange(patch_size/2,nrows*patch_size,patch_size)
    X, Y = np.meshgrid(x, y)
    for i,parameter_id in enumerate(traction_d.keys()):
        
        parameter = traction_d[parameter_id]
        title = f"{parameter_id} Stress across fault $\phi={strike_val}\degree$"
        # parameter Stress-related is already in correct order ie. row-wise as obtained from Okada
        im = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k',
                          cmap='bwr',
                          norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
        
        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        axes[i].margins(0)
                
        fig.colorbar(im, ax=axes[i],shrink=0.9,label='{}'.format('Stress (MPa)'))
        axes[i].plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markersize=18,markeredgecolor='black')
    
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
    fig.suptitle(f'Model:{model_type}',x=0.30,fontweight='bold')
    plt.tight_layout()
    plot_path = os.path.join(case_files_dir,f'{filename}.jpg')
    fig.savefig(plot_path)
    
 
#normal_shear_stress_map('output_Kinematic_Iquique_mean_model_over_fault_default_strike',strike_arr,dip_arr,'default',nrows,ncols,17,hypocenter)    
#normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_default_strike',strike_arr,dip_arr,'default',nrows,ncols,29,hypocenter)    
#normal_shear_stress_map('output_Kinematic_Tohoku_mean_model_over_fault_default_strike',strike_arr,dip_arr,'default',nrows,ncols,29,hypocenter)    


strike_arr = 90*np.ones(strike_arr.shape)
normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_90_strike','Static',strike_arr,dip_arr,'90',nrows,ncols,29,hypocenter)  

# normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_0_strike','Normal',strike_arr,dip_arr,'90',nrows,ncols,29,hypocenter)  
# normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_0_strike','Along_strike_Shear',strike_arr,dip_arr,'0',nrows,ncols,29,hypocenter)  
# normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_0_strike','Along_dip_Shear',strike_arr,dip_arr,'0',nrows,ncols,29,hypocenter)  


# normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_90_strike','Normal',strike_arr,dip_arr,'90',nrows,ncols,29,hypocenter)  
# normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_90_strike','Along_strike_Shear',strike_arr,dip_arr,'90',nrows,ncols,29,hypocenter)  
# normal_shear_stress_map('output_Static_Tohoku_mean_model_over_fault_90_strike','Along_dip_Shear',strike_arr,dip_arr,'90',nrows,ncols,29,hypocenter)  


#normal_shear_stress_map('output_Kinematic_Tohoku_mean_model_over_fault_0_strike',np.zeros(9*24),dip_arr,'0',9,24,29)    
#normal_shear_stress_map('output_Kinematic_Tohoku_mean_model_over_fault_90_strike',90*np.ones(9*24),dip_arr,'90',9,24,29)  