# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 21:25:44 2024

@author: joanv
"""
from utils.Okada_Covariance import observable_cov,correlation_all
import os
import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
#Gorkha
# observable_cov('Gorkha','stress',1000)
# observable_cov('Gorkha','displacement',1000)

#Pedernales
# observable_cov('Pedernales','stress',1000)
# observable_cov('Pedernales','displacement',1000)

# # Iquique
# observable_cov('Iquique','stress',1000)
# observable_cov('Iquique','displacement',1000)

# # Illapel
# observable_cov('Illapel','stress',1000)
# observable_cov('Illapel','displacement',1000)

# #Tohoku
# observable_cov('Tohoku','stress',1000,model_type='kinematic')
# observable_cov('Tohoku','displacement',1000,model_type='kinematic')

# # Tohoku
# observable_cov('Tohoku','stress',1000,model_type='static')
# observable_cov('Tohoku','displacement',1000,model_type='static')


# correlation_all('stress',1000)
# correlation_all('displacement',1000)


observable = 'displacement'
nsamples = 1000
model_type = 'kinematic'
var = {}
for name in ['Tohoku']:
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
    
    def order_corr(M,indeces):
        corr_ordered = M[indeces,:]
        corr_ordered = corr_ordered[:,indeces]
        return corr_ordered
    
    corr1_ordered = order_corr(corr1,ind_dr_sorted)
    corr2_ordered = order_corr(corr2,ind_dr_sorted)
    corr3_ordered= order_corr(corr3,ind_dr_sorted)
        
    # 2: trench-normal displacement
    
    fig,ax = plt.subplots(dpi = 900)
    im = ax.imshow(corr3_ordered,cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax = 1))
    ax.axes.get_xaxis().set_ticks([])
    ax.axes.get_yaxis().set_ticks([])
    ax.set_title('Tohoku Correlation Matrix \n of Trench-normal Displacement',fontsize=9.5,fontweight='bold')
    cbar = fig.colorbar(im,ax=ax)
    cbar.set_ticks([-1,-0.5,0,0.5,1])
    cbar.ax.tick_params(labelsize=8)
    
    plt.savefig('Tohoku_corr_matrix.png')
    
    
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

