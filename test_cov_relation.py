# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 11:16:06 2024

@author: joanv
"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 


def observable_cov_GF(name,
                      observable,
                      geom,
                      nparameters,
                      patchsize,
                      nsamples,
                      RotAngle=None,
                      model_type='kinematic',
                      sample_type='mean',
                      rake=90):

    working_dir = os.getcwd()
    output_dir = os.path.join(working_dir,'output')
    EQ_dir = os.path.join(output_dir,name)
    model_dir =  os.path.join(EQ_dir,'model')
    model_type_dir =  os.path.join(model_dir,model_type)
    sample_type_dir = os.path.join(model_type_dir,sample_type)
    input_file_dir = os.path.join(sample_type_dir,f'{name}_{sample_type}_{model_type}_model.csv')
    
    nrows,ncols = geom[0],geom[1]
    # input data directory
    if name=='Pedernales':
        df = pd.read_csv(input_file_dir)
        df= df.drop(df.columns[0],axis=1)
        dip = np.flip(df['dip'].values.reshape(ncols,nrows).transpose(),axis=0).flatten()
        strike =np.flip(df['strike'].values.reshape(ncols,nrows).transpose(),axis=0).flatten()
    # binary data with slip model
    
    
    name_dir  = os.path.join(output_dir,f'{name}') 
    model_dir =  os.path.join(name_dir,'model')
    model_type_dir = os.path.join(model_dir, f'{model_type}')
    samples_dir = os.path.join(model_type_dir,f'{nsamples}_samples')
    bin_dir = os.path.join(samples_dir,'bin_data')
    bin_data_dir = os.path.join(bin_dir,f'{name}_{model_type}_n_100.dat')
    
    # txt file with static GF
    
    output_dir = os.path.join(working_dir,'GF')
    GF_type_dir =  os.path.join(output_dir,'static')
    parameter_dir = os.path.join(GF_type_dir,f'{observable}')
    name_dir  = os.path.join(parameter_dir,f'{name}')
    model_type_dir = os.path.join(name_dir, f'{model_type}')
    txt_data_dir = os.path.join(model_type_dir,f'{name}_{model_type}_Static_GF_displacement.txt')
    
    bin_data = np.fromfile(bin_data_dir,'float').reshape((nparameters,nsamples))
    
    Np = nrows*ncols
    units = {'stress':'MPa','displacement':'m'} 
    
    
    sm = np.zeros( bin_data[:2*Np,:].shape)
    # for i in range(slip_model.shape[1]):
    #     modelss = np.flip(slip_model[:Npatches,i].reshape(nrows,ncols,order='F'),axis=0).flatten()
    #     modelds = np.flip(slip_model[Npatches:,i].reshape(nrows,ncols,order='F'),axis=0).flatten()
    #     sm[:Npatches,i] = modelss
    #     sm[Npatches:2*Npatches,i] = modelds
    
    
    
    for i in range(bin_data.shape[1]):
        model = bin_data[:,i] 
        model0 = np.copy(model)
        if name =='Gorkha':
            model[:Np],model[Np:2*Np] = model0[Np:2*Np],model0[:Np]
        '''
        Pedernales: requires call to mean_model.csv because of strike,dip
        '''
        if name=='Pedernales':
            ss0 = np.zeros(Np)
            ds0 = np.zeros(Np)
            for p in range(Np):
                RotAngle2 = RotAngle*((np.pi) / 180.)
                
                rotation = np.arctan2(np.tan(strike[p]) - np.tan(RotAngle2), 
                                    np.cos(dip[p])*(1.+np.tan(RotAngle2)*np.tan(strike[p])))
    
                # If RotAngle within ]90, 270], change rotation
                if RotAngle > 90. and RotAngle<=270.:
                    rotation += np.pi
    
                rp = model[p]    # rake-perpendicular
                ar = model[p+Np] # rake-parallel
    
                ss0[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                ds0[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
            Uperp = ss0
            Uparallel = ds0    
        else:
            
            rake0 = (rake - 90)*(np.pi/180)
            Uperp =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0) 
            Uparallel = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)  
        
        sm[:Np,i] = np.flip(Uperp.reshape(ncols,nrows).transpose(),axis=0).flatten()
        sm[Np:2*Np,i] = np.flip(Uparallel.reshape(ncols,nrows).transpose(),axis=0).flatten()
    cov_sm = np.cov(sm)
    
    GF_d = np.loadtxt(txt_data_dir)

    cov_d = np.matmul(GF_d,np.matmul(cov_sm,GF_d.transpose()))
    cov = cov_d
    nparameters = cov_d.shape[1]
    cov1 = cov[:nparameters//3,:nparameters//3]
    cov2 = cov[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    cov3 = cov[2*nparameters//3:,2*nparameters//3:]
    
    #cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
    # standard deviation (= square root of variance)
    std1 = np.sqrt(cov1.diagonal())
    std2 = np.sqrt(cov2.diagonal())
    std3 = np.sqrt(cov3.diagonal())
    
    
    
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    if observable=='stress':
        std= {'normal':std1,'along-strike':std2,'along-dip':std3}
    else:
        std= {'x':std1,'y':std2,'z':std3}
    
    # stress is defined on geometry for prescribed fault surface
    x0 = np.arange((1/2)*patchsize,ncols*patchsize,patchsize)
    y0 = np.arange(-(nrows-1/2)*patchsize,0,patchsize)
    
    n0cols,n0rows = len(x0),len(y0)
    # displacement is defined on customized geometry overlapping prescribed surface 
    
    xs = np.arange(patchsize/2 - (ncols//4)*patchsize, ncols*patchsize+(ncols//4)*patchsize,patchsize)
    ys = np.arange(-(nrows//4)*patchsize - (nrows-1/2)*patchsize, (nrows//4)*patchsize,patchsize)
    
    ncols,nrows = len(xs),len(ys)                     
    X, Y = np.meshgrid(xs, ys)
    
    for i,parameter_id in enumerate(std.keys()):
        
        parameter = std[parameter_id]
        if observable=='stress':
            title = f"$\sigma$ in {parameter_id} {observable} across fault"
        else:
            title = f"$\sigma$ in {parameter_id} surface {observable} "
    
        im0 = axes[i].pcolormesh(xs,ys,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow')
        parray = np.zeros((n0rows,n0cols))
        masked_arr = np.ma.masked_array(parray,parray<1)
        axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white')
    
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        axes[i].margins(0)        
        fig.colorbar(im0, ax=axes[i],shrink=0.9,label='{} ({})'.format(observable,units[observable]))
        axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold')
        axes[i].tick_params(labelsize=12)
    plt.tight_layout()
    
    # output dir 
    okada_observables_dir = os.path.join(working_dir,'cov_okada_observables')
    samples_dir = os.path.join(okada_observables_dir,f'{nsamples}_samples')
    EQ_dir = os.path.join(samples_dir,name)
    
    
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        figure_dir = os.path.join(model_type_dir,'figures')
        figure_type_dir = os.path.join(figure_dir,'with_GF')
    else:
        figure_dir = os.path.join(EQ_dir,'figures')
        figure_type_dir = os.path.join(figure_dir,'with_GF')

    figure_file_dir  = os.path.join(figure_type_dir,f'std_from_GF_{name}_{observable}_nsamples_{nsamples}.jpg')
    try:
        os.makedirs(figure_type_dir)

    except FileExistsError:
        pass
    fig.savefig(figure_file_dir)
    
    
observable_cov_GF('Gorkha','displacement',(9,18),650,10,100)
observable_cov_GF('Iquique','displacement',(11,12),533,17,100)
observable_cov_GF('Illapel','displacement',(10,17),682,18,100)
observable_cov_GF('Pedernales','displacement',(8,10),331,15,100,RotAngle=360-99)
observable_cov_GF('Tohoku','displacement',(9,24),866,29,100)
observable_cov_GF('Tohoku','displacement',(9,24),432,29,100,model_type='static')

    