# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 09:18:37 2023

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 16:56:23 2023

@author: joanv
"""


import numpy as np
import pandas as pd 
import pygmt
import struct
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Ellipse
from scipy.stats import skew
from matplotlib.colors import TwoSlopeNorm

class Static_model():
    directory_folder = os.path.join(os.getcwd(),'static')
    directory_data =  'csv_data'
    directory_plots = 'plots'
    path_data = os.path.join(directory_folder, directory_data)  
    path_plots = os.path.join(directory_folder, directory_plots)
    
    def __init__(self,
                  name:str,
                  type_model:str,
                  model_file:str,
                  shape_model_file:tuple,
                  geom_file:str,
                  shape_geom:tuple,
                  patch_size:float,
                  precision='double',
                  sampling=False,
                  nsamples=20000):
        
        self.name = name
        self.type_model = type_model
        self.patch_size = patch_size
        self.nsamples = nsamples 
        
        ### Reading Source Model file ###
        
        data = np.fromfile(model_file,precision).reshape(shape_model_file)
        rng = np.random.default_rng()    # uniform sampler
        # uniform sampling `nsamples` from `data`
        
        if sampling==True:
            self.sampled_data = rng.choice(data,nsamples,axis=1) 
        else:
            self.sampled_data = data 
        
        # mean of each parameter at each node 
        self.mean =np.mean(self.sampled_data,axis=1) 
        
        # standard deviation of each parameter at each node 
        self.std = np.std(self.sampled_data,axis=1)
        
        #skewness
        self.skewness = skew(self.sampled_data,axis=1)
        ### Reading Geometry  ###
    
        geometry  = pd.read_csv(geom_file,skiprows=2,sep=' ') # subject to change
        
        geometry = geometry[['lat','lon','depth','strike','dip']] # taking only columns of interest
        
        # redifining column names
        geometry = geometry.rename(columns ={'lat':'latitude','lon':'longitude','depth':'depth_km','strike':'strike','dip':'dip'})
        self.geometry_d = geometry.to_dict(orient='list')  # save each column as list
        
        self.shape_geom =shape_geom
        
        self.resampled_models_()
        self.mean_()
        self.std_()
        self.skew_()
        self.cov_m()  # covariance_matrix
        self.full_model_()

    def resampled_models_(self):
        data = self.sampled_data.flatten()
        try: 
            os.mkdir(Static_model.directory_folder)  
            print(f'Directory {Static_model.directory_folder}  just created')
        except FileExistsError:
            print(f'Directory {Static_model.directory_folder} already created')
        filename  = f"n_{self.nsamples}_sampled_models_{self.type_model}_{self.name}.dat"    
        output_file = os.path.join(Static_model.directory_folder,filename)
        out_file = open(output_file,"wb")
        s = struct.pack('d'*len(data), *data)
        out_file.write(s)
        out_file.close()       
        
    def mean_(self):
         
        # subject to change
        mean_d = {
            'U_perp':self.mean[:216],           
            'U_parallel':self.mean[216:432],
            'Slip':  np.sqrt(self.mean[:216]**2 + self.mean[216:432]**2)  
                    }
        self.mean_d = mean_d
   
    def std_(self):
        
        # subject to change
        std_d = {
            'std_U_perp':self.std[:216],
            'std_U_parallel':self.std[216:]
                }
        self.std_d = std_d
        self.mean_d.update(self.std_d)
        
    def skew_(self):
        skew_d = {
            'skew_U_perp':self.skewness[:216],
            'skew_U_parallel':self.skewness[216:]    
                }
        self.skew_d = skew_d
        self.mean_d.update(self.skew_d)
        
    def cov_m(self):
        '''
        

        Parameters
        ----------
        patch_id : int
            runs from 0 - 215.

        Returns
        -------
        Covariance for all data 

        '''
        self.cov_m = np.cov(self.sampled_data)
        
        std_cov = {
            'std_U_perp_cov':np.sqrt(np.diag(self.cov_m)[:216]),
            'std_U_parallel_cov':np.sqrt(np.diag(self.cov_m)[216:432])
                }
        cov_pp = self.cov_m[:216,216:]
        perp_term = std_cov['std_U_perp_cov']*self.mean_d['U_perp']/self.mean_d['Slip']
        parallel_term = std_cov['std_U_parallel_cov']*self.mean_d['U_parallel']/self.mean_d['Slip']
        cross_term = self.mean_d['U_perp']*self.mean_d['U_parallel']*np.diag(cov_pp)/(self.mean_d['Slip']**2)
        std_U = np.sqrt(perp_term**2 + parallel_term**2 + 2*cross_term)
        std_cov['std_U_cov'] = std_U
        self.mean_d.update(std_cov)

        """
        U_perp = self.sampled_data[patch_id,:]
        U_parallel = self.sampled_data[patch_id + 216,:]
        m = np.vstack(U_perp,U_parallel)
        cov_d = np.cov(m)
        
        self.patch_cov = cov_d
        """
        
    def full_model_(self):
        '''
        

        Returns
        -------
        Modifies existing `self.mean_d` 
        Saves full_mean_model with both parameters and geometry

        '''
        self.mean_d.update(self.geometry_d)
        try: 
            os.makedirs(Static_model.path_data)  
            print(f'Directory {Static_model.path_data}  just created')
        except FileExistsError:
            print(f'Directory {Static_model.path_data} already created')
        
        data_path = os.path.join(Static_model.path_data,f'full_mean_{self.type_model}_model_{self.name}.csv')
        m_df = pd.DataFrame.from_dict(self.mean_d, orient='index') # handle `Value error exception` 
        m_df = m_df.transpose()
        m_df.to_csv(data_path)
       
 
    
    
    def distribution(self,parameter,color_bar_label,skew=False):
        nrows = self.shape_geom[0]
        ncols = self.shape_geom[1]
        #hypocenter = (self.mean_d['Hypo_as'],self.mean_d['Hypo_dd'])
        x = np.arange(self.patch_size/2,ncols*self.patch_size,self.patch_size)
        y = np.arange(self.patch_size/2,nrows*self.patch_size,self.patch_size)
        fig, ax = plt.subplots(figsize=(12,10),dpi=600)
        param = self.mean_d[parameter].reshape(ncols,nrows).transpose()
        if skew:
            
            minimum = min(param.min(),-param.max())
            maximum = max(-param.min(),param.max())
            im = ax.pcolormesh(x,y,param,edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=minimum, vmax=maximum))
        else:
            im = ax.pcolormesh(x,y,param,edgecolors='k',cmap='YlOrRd')
        X, Y = np.meshgrid(x, y)
        #ax.plot(hypocenter[0],hypocenter[1],marker='*',color='blue',markeredgecolor='black',markersize=25)
        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        ax.margins(0)
        if parameter=='Slip':
            U = self.mean_d['U_parallel'].reshape(ncols,nrows).transpose() # dip-slip displacement
            V = self.mean_d['U_perp'].reshape(ncols,nrows).transpose()     # strike-slip displacement (see p.4 Minson et al. (2013) Part II)
            q = ax.quiver(X,Y,V,U,scale=1.20,scale_units ='x', units='width',width=0.002,headwidth=4.5,headlength=6)
            ax.quiverkey(q, X=0.1, Y=1.1, U=100,
                          label='length = 100m', labelpos='E')
            factor = 1/q.scale    # scaling factor from UV unit to XY units
            offsetXY = np.stack((V,U),axis=2).flatten().reshape(nrows*ncols,2)   
            
            # minus sign in Y coordinate corrects for y-axis inversion
            
            offsetXY = np.array([[xy[0],-xy[1]] for xy in offsetXY])  
            coord = q.XY + offsetXY*factor
            # width and height had to be formatted in row-like order (initially in column-like order)
            ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
                            width=self.mean_d['std_U_perp'].reshape((nrows,ncols),order='F').flatten()[i]*factor, 
                            height=self.mean_d['std_U_parallel'].reshape((nrows,ncols),order='F').flatten()[i]*factor,
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
            os.makedirs(Static_model.path_plots)  
            print(f'Directory {Static_model.path_plots}  just created')
        except FileExistsError:
            print(f'Directory {Static_model.path_plots} already created')
        plot_path = os.path.join(Static_model.path_plots,f'{self.type_model}_{parameter}_{self.name}.jpg')
        fig.savefig(plot_path)



