# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 09:42:56 2025

@author: joanv
"""


import matplotlib.pyplot as plt 
import numpy as np
import matplotlib as mpl
import pandas as pd


import h5py

import os 
from matplotlib.colors import TwoSlopeNorm

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

observable = 'displacement'
nsamples = 100
model_type = 'kinematic'

var = {}



INDECES = {'Tohoku':[],'Illapel':[],'Iquique':[],'Gorkha':[],'Pedernales':[]}
GRIDS = {'Tohoku':[],'Illapel':[],'Iquique':[],'Gorkha':[],'Pedernales':[]}

names = ['Tohoku','Iquique','Pedernales']
nx = len(names)

remainder = np.arange(int(nx*3*2))
fig = plt.figure(figsize=(int(nx*4),4))

gs0 = fig.add_gridspec(3,int(nx*2),wspace=0.35)

for n,name in enumerate(names):
    n1 = 2*n
    n2 = 2*n + 1 
    x = remainder[remainder%int(nx*2)==n1]
    y = remainder[remainder%int(nx*2)==n2]
    z = np.concatenate((x,y))
    z = list(z)

    INDECES[name].append(z)

for n,name in enumerate(names):
    indeces = INDECES[name][0]
    for ind in indeces[:len(indeces)//2]:
        gs = gs0[ind].subgridspec(5,1,hspace=0.0)
        GRIDS[name].append(gs)
        
    for ind in indeces[len(indeces)//2:]:
         gs2 = gs0[ind].subgridspec(1,1,hspace=0.0)
         GRIDS[name].append(gs2)



'''
gs02 = gs0[2].subgridspec(1,1,hspace=0.0)
gs03 = gs0[3].subgridspec(1,1,hspace=0.0)
gs04 = gs0[4].subgridspec(1,1,hspace=0.0)
'''


for q,name in enumerate(names):
    print(q)
    working_dir = os.getcwd()        
    mean_csv_file_dir = os.path.join(working_dir,f'INPUT/{name}/model/kinematic/{nsamples}_samples/mean/{name}_mean_{model_type}_model.csv')
    df = pd.read_csv(mean_csv_file_dir)
    df= df.drop(df.columns[0],axis=1)    
    
    nrows, ncols = models[name]['geom'][0], models[name]['geom'][1]
    patch = models[name]['patch']
    nparam = models[name]['nparam']
    
       

    xstn,ystn = set_stn(2,4,patch,control=1)
    
    nrows_d,ncols_d = len(ystn),len(xstn)
    Xstn,Ystn = np.meshgrid(xstn,ystn)
    xstn_flat,ystn_flat =  Xstn.flatten(),Ystn.flatten()

    h5file_dir = os.path.join(working_dir,f'OUTPUT/{name}/model/kinematic/EDKS/{nsamples}_samples/EDKS_{name}_{observable}_nsamples_{nsamples}_parallel.h5')
    
    f = h5py.File(h5file_dir,'r')

    dset = np.array(f[observable])
    # covariance matrix

    mean = np.mean(dset,axis=0)
    corr = np.corrcoef(dset.transpose())
    nparameters = dset.shape[1]
    
    mean_dx = mean[:nparameters//3]
    mean_dy = mean[nparameters//3:2*nparameters//3]
    mean_dz = mean[2*nparameters//3:]
    
    mean_d = {'x':mean_dx,'y':mean_dy,'z':mean_dz}
    
   
    # covariance of each observable (x,y,z) 
    

    corr1 = corr[:nparameters//3,:nparameters//3]
    corr2 = corr[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
    corr3 = corr[2*nparameters//3:,2*nparameters//3:]
    
    corr_d = {'x':corr1, 'y':corr2,'z':corr3}



    d = np.sqrt(mean_dx**2 + mean_dy**2 + mean_dz**2)
    max_d_patch = np.argmin(abs(d-max(d)))
    d = np.reshape(d,(nrows_d,ncols_d))
    corr_d_ordered = {'x':np.zeros_like(corr1), 'y':np.zeros_like(corr2),'z':np.zeros_like(corr3)}
    dr_ordered = {'x':np.zeros_like(corr1), 'y':np.zeros_like(corr2),'z':np.zeros_like(corr3)}
    r_all = np.column_stack((xstn_flat,ystn_flat))
    for comp in 'xyz':
        for i in range(ncols_d*nrows_d):
            
            col0_max_d,row0_max_d = max_d_patch%ncols_d,max_d_patch//ncols_d
            
            
            # displacement is defined on customized geometry overlapping prescribed surface 
        
        
            
            x_target = xstn_flat[i]
            y_target = ystn_flat[i]
            
            
            r_target = np.array([x_target,y_target])
    
        
            dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2)
            
            dr_sorted = np.sort(dr)
            dr_ordered[comp][i,:] = dr_sorted/max(dr_sorted)
            ind_dr_sorted = np.argsort(dr)
            ind_d_sorted = np.argsort(mean_d[comp])
            corr_d_ordered[comp][i,:] = corr_d[comp][i,ind_dr_sorted]  
    
    
        for i in range(ncols_d*nrows_d):
            corr_d_ordered[comp][:,i] = corr_d_ordered[comp][ind_d_sorted,i]  
        # 2: trench-normal displacement
    
    drmax  = round(max(dr_sorted),0)
    gsleft = GRIDS[name][:3]
    gsright = GRIDS[name][3:]
    indleft = INDECES[name][:3]
    indright = INDECES[name][3:]
    for l,comp in enumerate(['x','y','z']): 

        for i in range(5):
            ax = fig.add_subplot(gsleft[l][i])
            size = ncols_d*nrows_d
            value = int(i*size//4)
            if i==4:
                value = -1
            arr = np.sort(mean_d[comp])[value]
            disp = round(arr,2)
            ax.text(0.785,0.74,'%s'%disp + 'm',transform=ax.transAxes,fontsize=4)
            ax.scatter(dr_ordered[comp][value,:],corr_d_ordered[comp][value,:],s=1e-4,marker='.')
            ax.set_ylim((-1.1,1.1)) 
            ax.tick_params(labelsize=5.5)
            ax.axhline(y=0,color='grey',ls='--',lw=0.1)
            if i == 0 and l==0:
                ax.text(1.0,1.5,f'{name}',transform=ax.transAxes,fontweight='bold',fontsize=9)

                ax.text(-0.025,1.5,'abcde'[q],transform=ax.transAxes,fontweight='bold',fontsize=9)
                
                ax.set_ylabel('Corr',fontsize=4)
            if i!=0 or l>0:
                ax.set_yticks([])
                
            if i!=0 and l>0:
                ax.set_yticks([])
                
            if i<4:
                ax.set_xticks([])
        
            if l <2:
                ax.set_xticks([])
            if i==4  and l==2:
                ax.set_xlabel('Normalized distance (' + '$r_{max}$ ~ '+ '%4.d'%(drmax) + 'km)' ,fontsize=4) 
            if name=='Tohoku' and i==2:
                ax.text(-0.3, 0.5,  r'$d_{}$'.format(comp) ,fontsize=10,fontweight = 'bold',ha='center',rotation = 90,transform=ax.transAxes)

                
                
                
        ax = fig.add_subplot(gsright[l][0])
        im = ax.imshow(corr_d_ordered[comp],cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax = 1))
        ax.axes.get_xaxis().set_ticks([])
        ax.axes.get_yaxis().set_ticks([])
        ax.text(-0.09, 0.7, r'$\longleftarrow$' + ' Increasing ' + r'$d_{}$'.format(comp)  , fontsize=4,va='center', rotation=90,transform=ax.transAxes) 
                #axes[k//2][k%2].annotate('', xy=(80, 0.), xytext=(0, 0.), arrowprops=dict(facecolor='black',arrowstyle='->'))
        ax.text(0.52, 1.03, 'Increasing point-to-point distance $r$' + r'$\longrightarrow$',  fontsize=4,ha='center',transform=ax.transAxes)
        #ax.set_title(f'{name} Correlation Matrix {comp} Displacement',fontsize=9.5,fontweight='bold')
        cbar = fig.colorbar(im,ax=ax,shrink=0.6,pad=0.07)
        cbar.ax.set_title('Corr.',fontsize=4,y=-0.2,pad=2)
        cbar.set_ticks([-1,-0.5,0,0.5,1])
        cbar.ax.tick_params(labelsize=5.5)
            
    
fig.savefig('all_spatial_corr_displacement.pdf')
                                                                                         # scaled_dr = (dr_sorted - np.min(dr_sorted))/(np.max(dr_sorted) - np.min(dr_sorted))
    # if name =='Tohoku':
    #     model_type_dir = os.path.join(EQ_dir,model_type)
    #     figure_type_dir = os.path.join(model_type_dir,'figures')
    #     file_type_dir = os.path.join(figure_type_dir,'direct_computation')
    # else:
    #     figure_type_dir = os.path.join(EQ_dir,'figures')
    #     file_type_dir = os.path.join(figure_type_dir,'direct_computation')
    
    # try:
    #     os.makedirs(file_type_dir)

    # except FileExistsError:
    #     pass
    
    # if observable=='stress':
    #     corr_max = {'normal':corr1[:,max_observable_patch],
    #                 'along-strike':corr2[:,max_observable_patch],
    #                 'along-dip':corr3[:,max_observable_patch]}
    # else:
    #     corr_max = {'along-strike':corr1[:,max_observable_patch],
    #                 'trench-normal':corr2[:,max_observable_patch],
    #                 'vertical':corr3[:,max_observable_patch]}
    

    
    # ### correlation  
    # figure_dir = os.path.join(file_type_dir ,f'corr_curve_{name}_{observable}_nsamples_{nsamples}.jpg')

    # if observable=='stress':
    #     corr_max = {'normal':corr1[:,max_observable_patch],
    #                 'along-strike':corr2[:,max_observable_patch],
    #                 'along-dip':corr3[:,max_observable_patch]}
    # else:
    #     corr_max = {'along-strike':corr1[:,max_observable_patch],
    #                 'trench-normal':corr2[:,max_observable_patch],
    #                 'vertical':corr3[:,max_observable_patch]}
    
    
    
    # var[name]={'corr_max':corr_max,'dr_sorted':scaled_dr,'id_sorted':ind_dr_sorted} 
