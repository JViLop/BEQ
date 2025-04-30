# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 13:56:52 2024

@author: joanv
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
import h5py
import os


observable = 'stress '
units = {'stress':'MPa','displacement':'m'} 
nobservables = 3

def observable_cov(name,observable,nsamples,model_type='kinematic'):
    
    
    
    
    working_dir = os.getcwd()
    
    output_dir = os.path.join(working_dir,'output')
    EQ_dir = os.path.join(output_dir,f'{name}')
    model_dir = os.path.join(EQ_dir,'model')
    model_type_dir = os.path.join(model_dir,model_type)
    mean_dir = os.path.join(model_type_dir,'mean')
    mean_csv_file_dir = os.path.join(mean_dir,f'{name}_mean_{model_type}_model.csv')
    df = pd.read_csv(mean_csv_file_dir)
    df= df.drop(df.columns[0],axis=1)    
    
    try:
        hypocenter = [df['Hypo_as'][0],-df['Hypo_dd'][0]]
    except:
        print('No Hypocenter available')
    
    okada_observables_dir = os.path.join(working_dir,'okada_observables_ensemble')
    samples_dir = os.path.join(okada_observables_dir,f'{nsamples}_samples')
    EQ_dir = os.path.join(samples_dir,name)
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        file_type_dir = os.path.join(model_type_dir,'binary_data')
        h5file_dir = os.path.join(file_type_dir,f'{name}_{observable}_nsamples_{nsamples}.h5')
    else:
        file_type_dir = os.path.join(EQ_dir,'binary_data')
        h5file_dir = os.path.join(file_type_dir,f'{name}_{observable}_nsamples_{nsamples}.h5')
    
    f = h5py.File(h5file_dir,'r')

    dset = np.array(f[observable])
    # covariance matrix
    mean0 = np.mean(dset,axis=0)

    cov = np.cov(dset.transpose())
    corr = np.corrcoef(dset.transpose())
    
    cov_okada_observables = os.path.join(working_dir,'cov_okada_observables')
    samples_dir = os.path.join(cov_okada_observables ,f'{nsamples}_samples')
    EQ_dir = os.path.join(samples_dir,name)
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        file_type_dir = os.path.join(model_type_dir,'cov_txt')
    else:
        file_type_dir = os.path.join(EQ_dir,'cov_txt')
    
    try:
        os.makedirs(file_type_dir)
        
    except FileExistsError:
        pass
    
    file_dir = os.path.join(file_type_dir ,f'cov_matrix_{name}_{observable}_nsamples_{nsamples}.txt')
    np.savetxt(file_dir,cov)
    
    
    # covariance of each observable (x,y,z) 
    nparameters = dset.shape[1]
    cov1 = cov[:nparameters//3,:nparameters//3]
    cov2 = cov[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    cov3 = cov[2*nparameters//3:,2*nparameters//3:]
    
    corr1 = corr[:nparameters//3,:nparameters//3]
    corr2 = corr[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    corr3 = corr[2*nparameters//3:,2*nparameters//3:]
    
    #cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
    # standard deviation (= square root of variance)
    std1 = np.sqrt(cov1.diagonal())
    std2 = np.sqrt(cov2.diagonal())
    std3 = np.sqrt(cov3.diagonal())
    

        


    # stress is defined on geometry for prescribed fault surface
    x0 = np.array(f['x0_geometry'][0])/1e3
    y0 = np.array(f['y0_geometry'][0])/1e3# modified
    
    n0cols,n0rows = len(x0),len(y0)
    
    slip = np.flip(df['Slip'].values.reshape(n0cols,n0rows).transpose(),axis=0).flatten()
    max_slip_patch = np.argmin(abs(slip-max(slip)))
    col0_max_slip,row0_max_slip = max_slip_patch%n0cols,max_slip_patch//n0cols
    
    
    # displacement is defined on customized geometry overlapping prescribed surface 
    
    x = np.array(f['x_geometry'][0])/1e3
    y = np.array(f['y_geometry'][0])/1e3# modified
    
    ncols,nrows = len(x),len(y)
    extra_cols = (ncols-n0cols)//2
    extra_rows = (nrows-n0rows)//2
    ncol_target = col0_max_slip + extra_cols 
    nrow_target = row0_max_slip + extra_rows
            
    max_observable_patch = nrow_target*ncols + ncol_target   
    print(max_slip_patch,col0_max_slip,row0_max_slip)
    X, Y = np.meshgrid(x, y)
    
    
    xflat = X.flatten()
    yflat = Y.flatten()
    x_target = xflat[max_observable_patch]
    y_target = yflat[max_observable_patch]
    
    r_target = np.array([x_target,y_target])
    r_all = np.column_stack((xflat,yflat))

    dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2)
    
    dr_sorted = np.sort(dr)
    ind_dr_sorted = np.argsort(dr)
    
    scaled_dr = (dr_sorted - np.min(dr_sorted))/(np.max(dr_sorted) - np.min(dr_sorted))
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        figure_type_dir = os.path.join(model_type_dir,'figures')
        file_type_dir = os.path.join(figure_type_dir,'direct_computation')
    else:
        figure_type_dir = os.path.join(EQ_dir,'figures')
        file_type_dir = os.path.join(figure_type_dir,'direct_computation')
    
    try:
        os.makedirs(file_type_dir)

    except FileExistsError:
        pass
    figure_dir = os.path.join(file_type_dir ,f'std_{name}_{observable}_nsamples_{nsamples}.jpg')

    if observable=='stress':
        std= {'normal':std1,'along-strike':std2,'along-dip':std3}
        corr_max = {'normal':corr1[:,max_observable_patch],
                    'along-strike':corr2[:,max_observable_patch],
                    'along-dip':corr3[:,max_observable_patch]}
    else:
        std= {'x':std1,'y':std2,'z':std3}
        corr_max = {'x':corr1[:,max_observable_patch],
                    'y':corr2[:,max_observable_patch],
                    'z':corr3[:,max_observable_patch]}
    
    print(max_observable_patch,ncol_target,nrow_target)
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    for i,parameter_id in enumerate(std.keys()):
            
      
        parameter = std[parameter_id]
        if observable=='stress':
            title = f"{name} uncertainty in {parameter_id} {observable} across fault"
        else:
            title = f"{name} uncertainty in {parameter_id} surface {observable} "

        im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow')
        parray = np.zeros((n0rows,n0cols))
        masked_arr = np.ma.masked_array(parray,parray<1)
        axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
        axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1)
        axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=3)
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        try:
            axes[i].plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
        except:
            pass 
        axes[i].margins(0)        
        fig.colorbar(im0, ax=axes[i],shrink=0.9,label='{} ({})'.format(observable,units[observable]))
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=12)
    plt.tight_layout()    
    fig.savefig(figure_dir)
    
    ### correlation  
    figure_dir = os.path.join(file_type_dir ,f'corr_curve_{name}_{observable}_nsamples_{nsamples}.jpg')

    if observable=='stress':
        std= {'normal':std1,'along-strike':std2,'along-dip':std3}
        corr_max = {'normal':corr1[:,max_observable_patch],
                    'along-strike':corr2[:,max_observable_patch],
                    'along-dip':corr3[:,max_observable_patch]}
    else:
        std= {'x':std1,'y':std2,'z':std3}
        corr_max = {'x':corr1[:,max_observable_patch],
                    'y':corr2[:,max_observable_patch],
                    'z':corr3[:,max_observable_patch]}
    
    print(max_observable_patch,ncol_target,nrow_target)
    fig, axes = plt.subplots(3,1,figsize=(5,8),dpi=600)
    for i,parameter_id in enumerate(corr_max.keys()):
        parameter = corr_max[parameter_id]    
        parameter = parameter[ind_dr_sorted]
    
        if observable=='stress':
            title = f"{name} Correlation as function of distance in {parameter_id} {observable} across fault"
        else:
            title = f"{name} Correlation as function of distance in {parameter_id} surface {observable} "

        axes[i].scatter(scaled_dr, parameter,s=3)
        axes[i].axhline(y = 0,linestyle='dashed',color='grey',alpha=0.5)
        # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        
        axes[i].margins(0)        
        axes[i].set_ylim(-1.2,1.2)
        axes[i].set_ylabel('Correlation',fontsize=8)
        axes[i].set_xlabel('Distance from patch (km)',fontsize=8)
        # axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold',fontsize=8)
        axes[i].tick_params(labelsize=8)
    plt.tight_layout()    
    fig.savefig(figure_dir)
    
    
 
    ###
    

    figure_dir = os.path.join(file_type_dir ,f'corr_max_patch_{name}_{observable}_nsamples_{nsamples}.jpg')
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    # correlation
    for i,parameter_id in enumerate(corr_max.keys()):
            
      
        parameter = corr_max[parameter_id]
        if observable=='stress':
            title = f"{name} correlation for max slip patch in {parameter_id} {observable} across fault"
        else:
            title = f"{name} correlation for max slip patch in {parameter_id} surface {observable} "

        
        im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),
                                 edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
        parray = np.zeros((n0rows,n0cols))
        masked_arr = np.ma.masked_array(parray,parray<1)
        axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
        axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1)
        axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        try:
            axes[i].plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
        except:
            pass 
        axes[i].margins(0)        
        fig.colorbar(im0, ax=axes[i],shrink=0.9,label='{} ({})'.format(observable,units[observable]))
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold',fontsize=11)
        axes[i].tick_params(labelsize=12)
    plt.tight_layout()   
    fig.savefig(figure_dir)
    
    
    
    ## validation figure dir ##
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        file_type_dir = os.path.join(model_type_dir,'validation_figures')
    else:
        file_type_dir = os.path.join(EQ_dir,'validation_figures')
    
    try:
        os.makedirs(file_type_dir)

    except FileExistsError:
        pass
    
    # if observable=='stress':
    #     observable0 ='along-dip '+ observable
    # else:
    #     observable0 ='z-'+observable
        
    figure_dir = os.path.join(file_type_dir ,f'{name}_mean_{observable}_nsamples_{nsamples}.jpg')

    
    # fig,ax = plt.subplots()
    # parameter = mean0[2*nparameters//3:].reshape(nrows,ncols)
    # im0 = ax.pcolormesh(x,y,parameter,edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
    # fig.colorbar(im0, ax=ax,shrink=0.6,label='{} ({})'.format(observable,units[observable]))
    # ax.set_ylabel('Down-dip distance (km)',fontsize=12)
    # ax.set_xlabel('Along-strike distance (km)',fontsize=12)
    # ax.set_aspect('equal', 'box')
    # title = f'{name} Mean {observable0} with {nsamples} samples'
    # ax.set_title(title,fontweight='bold')
    # ax.tick_params(labelsize=12)
    # plt.tight_layout()
    # plt.show() 
    # fig.savefig(figure_dir)
    if observable=='stress':
        
        meand = {'normal':mean0[:nparameters//3],'along-strike':mean0[nparameters//3:2*nparameters//3],'along-dip':mean0[2*nparameters//3:]}
    else:
        meand = {'x':mean0[:nparameters//3],'y':mean0[nparameters//3:2*nparameters//3],'z':mean0[2*nparameters//3:]}

        
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    for i,parameter_id in enumerate(meand.keys()):
        
        parameter = meand[parameter_id]
        if observable=='stress':
            title = f"{name} mean {parameter_id} {observable} across fault"
        else:
            title = f"{name} mean {parameter_id} {observable} "

     
 
        im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
        parray = np.zeros((n0rows,n0cols))
        masked_arr = np.ma.masked_array(parray,parray<1)
        axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
        axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1)
        axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        axes[i].margins(0)        
        try:
            axes[i].plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
        except:
            pass   
        
        fig.colorbar(im0, ax=axes[i],shrink=0.8,label='{} ({})'.format(observable,units[observable]))
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=12)
    plt.tight_layout()
    fig.savefig(figure_dir)    

    

    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        figure_type_dir = os.path.join(model_type_dir,'figures_std_spatial_col_variations')

    else:
        figure_type_dir = os.path.join(EQ_dir,'figures_std_spatial_col_variations')

    
    try:
        os.makedirs(figure_type_dir)

    except FileExistsError:
        pass
    figure_dir = os.path.join(figure_type_dir ,f'spatial_std_{name}_{observable}_nsamples_{nsamples}.jpg')
    
    fig,axes = plt.subplots(3,1,figsize=(6,8),dpi=600)
    for i,parameter_id in enumerate(meand.keys()):

        target_obs = meand[parameter_id].reshape(nrows,ncols)[:,ncol_target]
        std_up_target_obs = target_obs + std[parameter_id].reshape(nrows,ncols)[:,ncol_target]
        std_down_target_obs = target_obs - std[parameter_id].reshape(nrows,ncols)[:,ncol_target]
        axes[i].plot(y,target_obs,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='mean')

        axes[i].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        axes[i].legend()
        axes[i].set_xlabel('trench-perpendicular distance from trench')
        axes[i].set_ylabel(f'{observable}')
        axes[i].set_title(f'{name} trench-perpendicular mean {parameter_id} {observable} $\pm \sigma $')    
    plt.tight_layout()
    fig.savefig(figure_dir)    
    
    
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        figure_type_dir = os.path.join(model_type_dir,'figures_std_spatial_row_variations')

    else:
        figure_type_dir = os.path.join(EQ_dir,'figures_std_spatial_row_variations')

    
    try:
        os.makedirs(figure_type_dir)

    except FileExistsError:
        pass
    figure_dir = os.path.join(figure_type_dir ,f'spatial_std_{name}_{observable}_nsamples_{nsamples}.jpg')
    
    fig,axes = plt.subplots(3,1,figsize=(6,8),dpi=600)
    for i,parameter_id in enumerate(meand.keys()):

        target_obs = meand[parameter_id].reshape(nrows,ncols)[:,nrow_target]
        std_up_target_obs = target_obs + std[parameter_id].reshape(nrows,ncols)[:,nrow_target]
        std_down_target_obs = target_obs - std[parameter_id].reshape(nrows,ncols)[:,nrow_target]
        axes[i].plot(y,target_obs,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='mean')
        axes[i].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)

        # axes[i].plot(y,std_up_target_obs,'o', ls='-',color='red',label='mean $ + \sigma$')
        # axes[i].plot(y,std_down_target_obs,'o', ls='-',color='red',label='mean $ - \sigma$')
        axes[i].legend()
        axes[i].set_xlabel('trench-parallel distance from trench')
        axes[i].set_ylabel(f'{observable}')
        axes[i].set_title(f'{name} trench-parallel mean {parameter_id} {observable} $\pm \sigma $')
    
    plt.tight_layout()
    fig.savefig(figure_dir) 
    
    
    # corr_max_patch = corr[:,max_slip_patch]
    
    

def correlation_all(observable,nsamples,model_type='kinematic'):
        
        
        var = {}
        for name in ['Tohoku','Gorkha','Iquique','Illapel','Pedernales']:
            working_dir = os.getcwd()        
            output_dir = os.path.join(working_dir,'output')
            EQ_dir = os.path.join(output_dir,f'{name}')
            model_dir = os.path.join(EQ_dir,'model')
            model_type_dir = os.path.join(model_dir,model_type)
            mean_dir = os.path.join(model_type_dir,'mean')
            mean_csv_file_dir = os.path.join(mean_dir,f'{name}_mean_{model_type}_model.csv')
            df = pd.read_csv(mean_csv_file_dir)
            df= df.drop(df.columns[0],axis=1)    
            
            
            okada_observables_dir = os.path.join(working_dir,'okada_observables_ensemble')
            samples_dir = os.path.join(okada_observables_dir,f'{nsamples}_samples')
            EQ_dir = os.path.join(samples_dir,name)
            if name =='Tohoku':
                model_type_dir = os.path.join(EQ_dir,model_type)
                file_type_dir = os.path.join(model_type_dir,'binary_data')
                h5file_dir = os.path.join(file_type_dir,f'{name}_{observable}_nsamples_{nsamples}.h5')
            else:
                file_type_dir = os.path.join(EQ_dir,'binary_data')
                h5file_dir = os.path.join(file_type_dir,f'{name}_{observable}_nsamples_{nsamples}.h5')
            
            f = h5py.File(h5file_dir,'r')
    
            dset = np.array(f[observable])
            # covariance matrix
    

            corr = np.corrcoef(dset.transpose())
            
            cov_okada_observables = os.path.join(working_dir,'cov_okada_observables')
            samples_dir = os.path.join(cov_okada_observables ,f'{nsamples}_samples')
            EQ_dir = os.path.join(samples_dir,name)
            if name =='Tohoku':
                model_type_dir = os.path.join(EQ_dir,model_type)
                file_type_dir = os.path.join(model_type_dir,'cov_txt')
            else:
                file_type_dir = os.path.join(EQ_dir,'cov_txt')
            
            try:
                os.makedirs(file_type_dir)
                
            except FileExistsError:
                pass
            
            
            # covariance of each observable (x,y,z) 
            nparameters = dset.shape[1]
    
            corr1 = corr[:nparameters//3,:nparameters//3]
            corr2 = corr[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
            corr3 = corr[2*nparameters//3:,2*nparameters//3:]
            
            
        
    
            # stress is defined on geometry for prescribed fault surface
            x0 = np.array(f['x0_geometry'][0])/1e3
            y0 = np.array(f['y0_geometry'][0])/1e3# modified
            
            n0cols,n0rows = len(x0),len(y0)
            
            slip = np.flip(df['Slip'].values.reshape(n0cols,n0rows).transpose(),axis=0).flatten()
            max_slip_patch = np.argmin(abs(slip-max(slip)))
            col0_max_slip,row0_max_slip = max_slip_patch%n0cols,max_slip_patch//n0cols
            
            
            # displacement is defined on customized geometry overlapping prescribed surface 
            
            x = np.array(f['x_geometry'][0])/1e3
            y = np.array(f['y_geometry'][0])/1e3# modified
            
            ncols,nrows = len(x),len(y)
            extra_cols = (ncols-n0cols)//2
            extra_rows = (nrows-n0rows)//2
            ncol_target = col0_max_slip + extra_cols 
            nrow_target = row0_max_slip + extra_rows
                    
            max_observable_patch = nrow_target*ncols + ncol_target   
            X, Y = np.meshgrid(x, y)
            
            
            xflat = X.flatten()
            yflat = Y.flatten()
            x_target = xflat[max_observable_patch]
            y_target = yflat[max_observable_patch]
            
            r_target = np.array([x_target,y_target])
            r_all = np.column_stack((xflat,yflat))
    
            dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2)
            
            dr_sorted = np.sort(dr)
            ind_dr_sorted = np.argsort(dr)
            
            scaled_dr = (dr_sorted - np.min(dr_sorted))/(np.max(dr_sorted) - np.min(dr_sorted))
            if name =='Tohoku':
                model_type_dir = os.path.join(EQ_dir,model_type)
                figure_type_dir = os.path.join(model_type_dir,'figures')
                file_type_dir = os.path.join(figure_type_dir,'direct_computation')
            else:
                figure_type_dir = os.path.join(EQ_dir,'figures')
                file_type_dir = os.path.join(figure_type_dir,'direct_computation')
            
            try:
                os.makedirs(file_type_dir)
    
            except FileExistsError:
                pass

            if observable=='stress':
                corr_max = {'normal':corr1[:,max_observable_patch],
                            'along-strike':corr2[:,max_observable_patch],
                            'along-dip':corr3[:,max_observable_patch]}
            else:
                corr_max = {'along-strike':corr1[:,max_observable_patch],
                            'trench-normal':corr2[:,max_observable_patch],
                            'vertical':corr3[:,max_observable_patch]}
            
    
            
            ### correlation  
            figure_dir = os.path.join(file_type_dir ,f'corr_curve_{name}_{observable}_nsamples_{nsamples}.jpg')
    
            if observable=='stress':
                corr_max = {'normal':corr1[:,max_observable_patch],
                            'along-strike':corr2[:,max_observable_patch],
                            'along-dip':corr3[:,max_observable_patch]}
            else:
                corr_max = {'along-strike':corr1[:,max_observable_patch],
                            'trench-normal':corr2[:,max_observable_patch],
                            'vertical':corr3[:,max_observable_patch]}
            
            
            
            var[name]={'corr_max':corr_max,'dr_sorted':scaled_dr,'id_sorted':ind_dr_sorted} 
            
            
        path_all_curves = 'cov_okada_observables/1000_samples' 
        figure_dir = os.path.join(path_all_curves ,f'corr_curve_all_nsamples_{observable}_{nsamples}.jpg')
        colors = ['blue','red','orange','seagreen','brown']
        markers = ['.','^','s','p','D']
        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=800)
        for i,parameter_id in enumerate(corr_max.keys()):
            for j,(name,color,marker) in enumerate(zip(var.keys(),colors,markers)):
                parameter = var[name]['corr_max'][parameter_id]    
                parameter = parameter[var[name]['id_sorted']]
            
                if observable=='stress':
                    title = f"Correlation as function of distance in {parameter_id} fault {observable}"
                else:
                    title = f"Correlation as function of distance in {parameter_id} surface {observable} "
                scaled_dr =  var[name]['dr_sorted']
                axes[i].plot(scaled_dr, parameter,color = color,marker=marker,lw=0.20,ms=1.5,label=f'{name}')
                # axes[i].scatter(scaled_dr, parameter,color = color,marker=marker,s=1.5,label=f'{name}')
                axes[i].axhline(y = 0,linestyle='dashed',color='grey',alpha=0.25)
                # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
                #X0, Y0 = np.meshgrid(x0, y0)
                #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
                
                axes[i].margins(0)        
                axes[i].set_ylim(-1.2,1.2)
                axes[i].set_ylabel('Correlation',fontsize=8)
                axes[i].set_xlabel('Normalized Distance from patch',fontsize=8)
                # axes[i].set_aspect('equal', 'box')
                axes[i].set_title(title,fontweight='bold',fontsize=8)
                axes[i].tick_params(labelsize=8)
                axes[i].legend(loc="upper right",handlelength=1,ncol = 5,fontsize=6)
        plt.tight_layout()    
        fig.savefig(figure_dir)
    

