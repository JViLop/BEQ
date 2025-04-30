# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:55:34 2023

@author: joanv
"""

import h5py
import numpy as np
import pandas as pd 

import struct
import matplotlib.pyplot as plt
import os
from matplotlib.patches import Ellipse


Np=132 # Number of subfaults
file_dir = './Data_Iquique_model/step_056.h5'
f = h5py.File(file_dir,'r')
samples = f['Sample Set']

strike_slip = samples[:,:Np]
dip_slip = samples[:,Np:2*Np]
nuisance_p = samples[:,2*Np:2*Np+3]
rise_time = samples[:,2*Np+3:3*Np+3]
rupture_velocity = samples[:,3*Np+3:4*Np+3] #must change see 
hypocenter = samples[:,4*Np+3:]


geom_file = './Data_Iquique_model/Iquique_FaultGeometry_v01-10km.3DsquareGrid'
geometry  = np.loadtxt(geom_file) # subject to change      

geometry_d =dict()
geometry_d['latitude'] =geometry[:,1]
geometry_d['longitude'] = geometry[:,0]
geometry_d['depth_km'] = geometry[:,4]
geometry_d['strike'] = geometry[:,5]
geometry_d['dip'] = geometry[:,6]


data = np.array(samples).transpose()



class Kinematic():
    
    directory_folder = os.path.join(os.getcwd(),'kinematic_class')
    directory_data =  'csv_data'
    directory_plots = 'plots'
    path_data = os.path.join(directory_folder, directory_data)  
    path_plots = os.path.join(directory_folder, directory_plots)
    
    def __init__(self,
                  name:str,
                  type_model:str,
                  model_file:str,
                  file_type:str,
                  shape_model_file:tuple,
                  geom_file:str,
                  shape_geom:tuple,
                  npatches:int,
                  patch_size:float,
                  precision='double',
                  sampling=False,
                  nsamples=20000):
        
        self.name = name
        self.type_model = type_model
        self.patch_size = patch_size
        self.npatches = npatches
        self.nsamples = nsamples 
        

        
        if file_type=='h5':
            f = h5py.File(model_file,'r')
            data = f['Sample Set']
            data = np.array(samples).transpose()

        ### Reading Source Model file ###
        elif file_type=='binary':  
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
        
        
        ### Reading Geometry  ###
        try: 
            
            geometry  = pd.read_csv(geom_file,skiprows=2,sep=' ') # subject to change            
            geometry = geometry[['lat','lon','depth','strike','dip']] # taking only columns of interest
            geometry = geometry.rename(columns ={'lat':'latitude','lon':'longitude','depth':'depth_km','strike':'strike','dip':'dip'})
            self.geometry_d = geometry.to_dict(orient='list')
        except:
            geometry  = np.loadtxt(geom_file) # subject to change      
            geometry_d =dict()
            geometry_d['latitude'] =geometry[:,1]
            geometry_d['longitude'] = geometry[:,0]
            geometry_d['depth_km'] = geometry[:,4]
            geometry_d['strike'] = geometry[:,5]
            geometry_d['dip'] = geometry[:,6]
            self.geometry_d = geometry_d
        # redifining column names
        
        # version before 
        #geometry = geometry.rename(columns ={'lat':'latitude','lon':'longitude','depth':'depth_km'})
        
        
        self.shape_geom =shape_geom
        self.resampled_models_()
        self.mean_()
        self.std_()
        self.cov_m()  # covariance_matrix
        self.full_model_()

    def resampled_models_(self):
        data = self.sampled_data.flatten()
        try: 
            os.mkdir(Kinematic.directory_folder)  
            print(f'Directory {Kinematic.directory_folder}  just created')
        except FileExistsError:
            print(f'Directory {Kinematic.directory_folder} already created')
        filename  = f"n_{self.nsamples}_sampled_models_{self.type_model}_{self.name}.dat"    
        output_file = os.path.join(Kinematic.directory_folder,filename)
        out_file = open(output_file,"wb")
        s = struct.pack('d'*len(data), *data)
        out_file.write(s)
        out_file.close()      
    

    def mean_(self):
        Np = self.npatches
        # subject to change
        mean_d = {
            'U_perp':self.mean[:Np],         # along-strike
            'U_parallel':self.mean[Np:2*Np], # along-dip
            'nuisance': self.mean[2*Np:2*Np+3],
            'Slip':  np.sqrt(self.mean[:Np]**2 + self.mean[Np:2*Np]**2),
            'Tr':self.mean[2*Np+3:3*Np+3],
            'Vr':self.mean[3*Np+3:4*Np+3],
            'Hypo_as': self.mean[4*Np+3:4*Np+4], # along-strike
            'Hypo_dd': self.mean[4*Np+4:] # down-dip
                    }
        self.mean_d = mean_d
   
    def std_(self):
        Np = self.npatches
        # subject to change
        std_d = {
            'std_U_perp':self.std[:Np],
            'std_U_parallel':self.std[Np:2*Np],
            'std_Tr':self.std[2*Np+3:3*Np+3],
            'std_Vr':self.std[3*Np+3:4*Np+3],
            'std_Hypo_as': self.std[4*Np+3:4*Np+4], # along-strike
            'std_Hypo_dd': self.std[4*Np+4:] # down-dip
                }
        self.std_d = std_d
        self.mean_d.update(self.std_d)
        
    def cov_m(self):
        '''
        

        Parameters
        ----------
    

        Returns
        -------
        Covariance for all data 

        '''
        self.cov_m = np.cov(self.sampled_data)
        
        std_cov = {
            'std_U_perp_cov':np.sqrt(np.diag(self.cov_m)[:Np]),
            'std_U_parallel_cov':np.sqrt(np.diag(self.cov_m)[Np:2*Np]),
            'std_Tr_cov':np.sqrt(np.diag(self.cov_m)[2*Np+3:3*Np+3]),
            'std_Vr_cov':np.sqrt(np.diag(self.cov_m)[3*Np+3:4*Np+3]),
            'std_Hyp_as_cov':np.sqrt(np.diag(self.cov_m)[4*Np+3:4*Np+4]),
            'std_Hyp_dd_cov':np.sqrt(np.diag(self.cov_m)[4*Np+4:]),
            }
        # computing std for U  = U(Uper,Uparallel)
        
        cov_pp = self.cov_m[:Np,Np:2*Np]
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
            os.makedirs(Kinematic.path_data)  
            print(f'Directory {Kinematic.path_data}  just created')
        except FileExistsError:
            print(f'Directory {Kinematic.path_data} already created')
        
        data_path = os.path.join(Kinematic.path_data,f'full_mean_{self.type_model}_model_{self.name}.csv')
       
        m_df = pd.DataFrame.from_dict(self.mean_d, orient='index') # handle `Value error exception` 
        m_df = m_df.transpose()
        m_df.to_csv(data_path)
       
 
    
    
    def distribution(self,parameter,color_bar_label,order_ok = True):
        hypocenter = (self.mean_d['Hypo_as'],self.mean_d['Hypo_dd'])
        nrows = self.shape_geom[0]
        ncols = self.shape_geom[1]
        x = np.arange(self.patch_size/2,ncols*self.patch_size,self.patch_size)
        y = np.arange(self.patch_size/2,nrows*self.patch_size,self.patch_size)
        fig, ax = plt.subplots(figsize=(12,10),dpi=600)
        if order_ok:
            
            im = ax.pcolormesh(x,y, self.mean_d[parameter].reshape(nrows,ncols),edgecolors='k',cmap='YlOrRd')
        else: 
            im = ax.pcolormesh(x,y, self.mean_d[parameter].reshape(ncols,nrows).transpose(),edgecolors='k',cmap='YlOrRd')
        X, Y = np.meshgrid(x, y)
        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        ax.margins(0)
        ax.plot(hypocenter[0],hypocenter[1],marker='*',color='blue',markeredgecolor='black',markersize=25)
        if parameter=='Slip':
            U = self.mean_d['U_parallel'].reshape(ncols,nrows).transpose() # dip-slip displacement
            V = self.mean_d['U_perp'].reshape(ncols,nrows).transpose()     # strike-slip displacement (see p.4 Minson et al. (2013) Part II)
            q = ax.quiver(X,Y,V,U,scale=0.75,scale_units ='x', units='width',width=0.0035,headwidth=2.5,headlength=2.5)
            ax.quiverkey(q, X=0.1, Y=1.1, U=15,
                          label='length = 15m', labelpos='E')
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
            os.makedirs(Kinematic.path_plots)  
            print(f'Directory {Kinematic.path_plots}  just created')
        except FileExistsError:
            print(f'Directory {Kinematic.path_plots} already created')
        plot_path = os.path.join(Kinematic.path_plots,f'{self.type_model}_{parameter}_{self.name}.jpg')
        fig.savefig(plot_path)

model_file = './Data_Iquique_model/step_056.h5'
geom_file = './Data_Iquique_model/Iquique_FaultGeometry_v01-10km.3DsquareGrid'

Iquique_kinematic = Kinematic('Iquique', 'kinematic',model_file,'h5',(1,1),geom_file,(11,12),132,17)
Iquique_kinematic.distribution('U_parallel',"$U_{||}$ (m)",order_ok = False)
Iquique_kinematic.distribution('U_perp',"$U_{\perp}$ (m)",order_ok = False)
Iquique_kinematic.distribution('Slip',"U (m)",order_ok = False)
Iquique_kinematic.distribution('Vr', '$V_{r}$ (km/s)',order_ok = False)
Iquique_kinematic.distribution('Tr', '$T_{r}$ (s)',order_ok = False)

Iquique_kinematic.distribution('std_U_parallel_cov',"$\sigma(U_{||})$ (m)",order_ok = False)
Iquique_kinematic.distribution('std_U_perp_cov',"$\sigma(U_{\perp})$ (m)",order_ok = False)
Iquique_kinematic.distribution('std_U_cov',"$\sigma(U)$ (m)",order_ok = False)
Iquique_kinematic.distribution('std_Vr_cov', '$\sigma(V_{r})$ (km/s)',order_ok = False)
Iquique_kinematic.distribution('std_Tr_cov', '$\sigma(T_{r})$ (s)' ,order_ok = False)



