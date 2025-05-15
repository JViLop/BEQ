# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 14:29:32 2024

@author: joanv
"""

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
from matplotlib.colors import TwoSlopeNorm
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

class EQ_model:
    def __init__ (self,name:str,
                   model_type:str,
                   step:int,
                   file_type:str,
                   geom_file_shape:tuple,
                   patchsize:int,
                   header = {'lon':0,'lat':1,'depth':4,'strike':5,'dip':6},
                   rake = 90,
                   RotAngle=None,
                   nrows_skip=1,
                   model_file_shape=None,
                   nramp = 0,
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
            lon, lat, depth, strike, dip    default works for Iquique and Illapel
        Returns
        -------
        None.
)
    '''
        self.name = name
        self.model_type = model_type
        self.step = step
        self.file_type = file_type
       
        self.model_file_shape = model_file_shape
        self.sampling = sampling
        if self.sampling:
            self.nsamples = nsamples
        else:
            self.nsamples = 'all'

        self.header = header
        self.nrows_skip = nrows_skip
        self.nramp = nramp
        self.precision = precision
        
        self.geom_file_shape = geom_file_shape
        self.nrows = geom_file_shape[0]
        self.ncols = geom_file_shape[1]
        self.patchsize = patchsize
        self.Np= int(self.nrows*self.ncols) # number of patches
        
        
        self.x = np.arange(self.patchsize/2,self.ncols*self.patchsize,self.patchsize)
        self.y = np.arange(self.patchsize/2,self.nrows*self.patchsize,self.patchsize)
        self.ydown = np.arange(-(self.nrows-1/2)*self.patchsize,0,self.patchsize)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
        
        self.rake = rake          # used for Gorkha
        self.RotAngle = RotAngle  # used for Pedernales
        # parent directories
        self.EQs_dir = os.path.join(os.getcwd(),'EQ')
        self.EQ_dir =  os.path.join(self.EQs_dir,self.name)
        
        # INPUT model and geometry directories 
        
        self.geometry_dir = os.path.join(self.EQ_dir,'geometry')
        self.geometry_formatted_dir = os.path.join(self.geometry_dir,'geometry_formatted')
        self.geometry_formatted_dir_plot = os.path.join(self.EQ_dir,'geometry_txt_for_plotting')
        self.geometry_files = os.listdir(self.geometry_dir)
        self.geometry_file_name = [ i for i in self.geometry_files if not i.endswith("formatted")][0]
        self.geometry_file_dir = os.path.join(self.geometry_dir, self.geometry_file_name) # must contain one file only
        

        self.model_dir = os.path.join(self.EQ_dir,'model')
        self.model_type_dir = os.path.join(self.model_dir,self.model_type)
        self.model_step_dir = os.path.join(self.model_type_dir,f'step_{self.step}')
        self.model_file_name = os.listdir(self.model_step_dir)[0]
        self.model_file_dir = os.path.join(self.model_step_dir,self.model_file_name)
    
        
        # OUTPUT directories
        self.output_EQs_dir = os.path.join(os.getcwd(),'INPUT')
        self.output_EQ_dir =  os.path.join(self.output_EQs_dir,self.name)
        self.output_model_dir = os.path.join(self.output_EQ_dir,'model')
        self.output_model_type_dir = os.path.join(self.output_model_dir,self.model_type)
        self.output_model_step_dir = os.path.join(self.output_model_type_dir,f'step_{self.step}')
        if sampling:
            self.output_sampling_type = os.path.join(self.output_model_step_dir,f'{self.nsamples}_samples')
            self.output_mean_model_dir = os.path.join(self.output_sampling_type,'mean')
            self.output_plots_model_dir = os.path.join(self.output_sampling_type ,'figures')
            self.output_bin_data_dir = os.path.join(self.output_sampling_type,'bin_data')
        else:
            self.output_sampling_type = os.path.join(self.output_model_step_dir,'all_samples')    
            self.output_mean_model_dir = os.path.join(self.output_sampling_type,'mean') 
            self.output_plots_model_dir = os.path.join(self.output_sampling_type ,'figures')
            self.output_bin_data_dir = os.path.join(self.output_sampling_type,'bin_data')
    ### Reading Geometry  ###
        
        try: 
            
            geometry  = pd.read_csv(self.geometry_file_dir,skiprows=nrows_skip,sep=' ') # subject to change            
            geometry = geometry[['lat','lon','depth','strike','dip']] # taking only columns of interest
            self.geometry_d = geometry.to_dict(orient='list')
        
        except: 
            try:
                geometry  = np.loadtxt(self.geometry_file_dir,skiprows=nrows_skip) # subject to change      
                i_lat = header['lat']
                i_lon = header['lon']
                i_depth = header['depth']
                i_strike = header['strike']
                i_dip = header['dip']
                self.geometry_d =dict()
                self.geometry_d['lat'] =geometry[:,i_lat]
                self.geometry_d['lon'] = geometry[:,i_lon]
                self.geometry_d['depth'] = geometry[:,i_depth]
                self.geometry_d['strike'] = geometry[:,i_strike]
                self.geometry_d['dip'] = geometry[:,i_dip]
            except:
                print(f'No geometry for {self.name} so cannot be saved')
                pass
            
        
        self.geometry_save()
        # self.geometry_checker()
        self.model_data_access()
        self.mean()
        self.standard_deviation()
        # self.skewness()
        self.model_save()
        
    
    def geometry_txt(self):
        gd = self.geometry_d
        col = np.ones(self.Np)
        areas = (self.patchsize**2)*col
        ids = np.arange(self.Np,dtype='int')
        D = np.array([gd['lon'],gd['lat'],col,col,
                      gd['depth'],gd['strike'],gd['dip'],
                      areas,ids])
        D = D.transpose()
        hd = 'lon lat E[km] N[km] Dep[km] strike dip Area ID' 
        fmt =  '%1.3f','%1.3f','%1.3f','%1.3f','%1.3f','%1.3f','%1.3f', '%1.3f','%d'
        try:
            os.makedirs(self.geometry_formatted_dir_plot)
        except FileExistsError:
            pass
        file_dir = os.path.join(self.geometry_formatted_dir_plot,f'{self.name}_{self.model_type}_geometry.txt')
        np.savetxt(file_dir,D,header=hd,fmt=fmt)
        return 

    def geometry_save(self):
        try: 
            self.csv_file_dir = os.path.join(self.geometry_formatted_dir,f'{self.name}_geometry.csv')
            geom_df = pd.DataFrame.from_dict(self.geometry_d, orient='index') # handle `Value error exception` 
            geom_df0 = geom_df.transpose()
            
            try:
                os.makedirs(self.geometry_formatted_dir)
            except FileExistsError:
                pass
            
            geom_df0.to_csv(self.csv_file_dir)
        
        except:
            print('geometry ok')
            pass
           
    
    def geometry_checker(self):
        df = pd.read_csv(self.csv_file_dir)
        lat = df['lat'].values
        
        lon = df['lon'].values
        
        ids =df.iloc[:,0].values
        plt.scatter(lon,lat,color='red')
        for i,id0 in enumerate(ids):
            # try:
                
            #     dlon = lon[i+1] - lon[i]
            #     dlat = lat[i+1] - lat[i]
            #     # print(f'dlon{i}:',dlon)
            #     # print(f'dlat{i}:',dlat)
            #     #print(f'd{i}',np.sqrt(dlon**2+dlat**2))
                
            # except:
            #  pass   
            plt.annotate(str(int(id0)),(lon[i],lat[i]),(lon[i],lat[i]),fontsize=8)
        plt.title(f'{self.name} {self.model_type}')
        plt.axis('equal')
        plt.show()
            
    
    
    def model_data_access(self):
        if self.file_type=='h5':
            f = h5py.File(self.model_file_dir,'r')
            try:
                
                data = f['Sample Set']
                data = np.array(data).transpose()
            
            except: 
                data = hstacking(f['ParameterSets']).transpose()
                
    
        ### Reading Source Model file ###
        elif self.file_type=='binary':  
            data = np.fromfile(self.model_file_dir,self.precision).reshape(self.model_file_shape)
        
    
        
        rng = np.random.default_rng()    # uniform sampler
        # uniform sampling `nsamples` from `data`
        
        if self.sampling==True:
            self.sampled_data = rng.choice(data,self.nsamples,axis=1) 
            dat = self.sampled_data.flatten()
            try: 
                os.makedirs(self.output_bin_data_dir)    
            except FileExistsError:
               pass   
            
            filename  = f"{self.name}_{self.model_type}_n_{self.nsamples}.dat"    
            output_file = os.path.join(self.output_bin_data_dir,filename)
            out_file = open(output_file,"wb")
            s = struct.pack('d'*len(dat), *dat)
            out_file.write(s)
            out_file.close()     
        else:
            self.sampled_data = data
            dat = self.sampled_data.flatten()
            try: 
                os.makedirs(self.output_bin_data_dir)    
            except FileExistsError:
               pass   
            
            filename  = f"{self.name}_{self.model_type}_beta_{self.step}_n_{self.nsamples}_.dat"    
            output_file = os.path.join(self.output_bin_data_dir,filename)
            out_file = open(output_file,"wb")
            s = struct.pack('d'*len(dat), *dat)
            out_file.write(s)
            out_file.close()
        
        
    # mean of each parameter at each node 
    # mean = np.mean(sampled_data,axis=1) 

    # # standard deviation of each parameter at each node 
    # std = np.std(sampled_data,axis=1)
    
    # # skewness of each parameter at each node
    
    # skewness = skew(sampled_data,axis=1)
    
    # # covariance matrix 
    # cov_m = np.cov(sampled_data)
    def mean(self):
        Np = self.Np
        nramp = self.nramp
        mean = np.mean(self.sampled_data,axis=1)
        mean0 = np.copy(mean)
        if self.name =='Gorkha':
              mean[:Np],mean[Np:2*Np] = mean0[Np:2*Np],mean0[:Np]
        if self.name == 'Pedernales':
            ss = np.zeros(Np)
            ds = np.zeros(Np)
            strike = self.geometry_d['strike']*((np.pi)/180)
            dip = self.geometry_d['dip']*((np.pi)/180)
            for p in range(Np):
                RotAngle2 = self.RotAngle*((np.pi) / 180.)
                
                rotation = np.arctan2(np.tan(strike[p]) - np.tan(RotAngle2), 
                                    np.cos(dip[p])*(1.+np.tan(RotAngle2)*np.tan(strike[p])))
    
                # If RotAngle within ]90, 270], change rotation
                if self.RotAngle > 90. and self.RotAngle<=270.:
                    rotation += np.pi
    
                rp = mean[p]    # rake-perpendicular
                ar = mean[p+Np] # rake-parallel
    
                ss[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                ds[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
            Uperp = ss
            Uparallel = ds    
        else:
            rake0 = (self.rake - 90)*(np.pi/180)
            Uperp =  mean[:Np]*np.cos(rake0)  - mean[Np:2*Np]*np.sin(rake0) 
            Uparallel = mean[:Np]*np.sin(rake0)  + mean[Np:2*Np]*np.cos(rake0)      
        if self.model_type=='static':
                self.mean_d = {
                    'U_perp':Uperp,         # along-strike
                    'U_parallel':Uparallel, # along-dip
                    'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2)
                                }
        else:
            self.mean_d = {
                    'U_perp':Uperp,         # along-strike
                    'U_parallel':Uparallel, # along-dip
                    'Slip': np.sqrt(mean[:Np]**2 + mean[Np:2*Np]**2),
                    'Tr':mean[2*Np+nramp:3*Np+nramp],
                    'Vr':mean[3*Np+nramp:4*Np+nramp],
                    'Hypo_as': mean[4*Np+nramp:4*Np+nramp+1], # along-strike
                    'Hypo_dd': mean[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                            }
    
    
    def max_patch(self):
        #return 
        # slip = np.flip(self.mean_d['Slip'].reshape(self.ncols,self.nrows).transpose(),axis=0).flatten()

        slip = self.mean_d['Slip']
        max_slip_patch = np.argmin(abs(slip-max(slip)))
        col0_max_slip,row0_max_slip = max_slip_patch%self.ncols,max_slip_patch//self.ncols
        self.ncol_target = col0_max_slip 
        self.nrow_target = row0_max_slip 
        max_observable_patch = self.nrow_target*self.ncols + self.ncol_target 
        self.max_slip = max_observable_patch # index of max_slip_patch if flatten
        
        slip_flipped = np.flip(self.mean_d['Slip'].reshape(self.ncols,self.nrows).transpose(),axis=0).flatten()
        max_flipped_slip_patch = np.argmin(abs(slip_flipped-max(slip_flipped)))
        col_max_flipped_slip,row_max_flipped_slip = max_flipped_slip_patch%self.ncols,max_flipped_slip_patch//self.ncols
        self.ncol_target_flipped = col_max_flipped_slip 
        self.nrow_target_flipped = row_max_flipped_slip
        
    
        return 
        
    def slip_corr_f(self):
            
        U1 = self.sampled_data[:self.Np,:]
        U2 = self.sampled_data[self.Np:2*self.Np,:] 
        slip_ensemble = np.sqrt(U1**2+U2**2)
        
        corr = np.corrcoef(slip_ensemble)
        self.slip_corr = corr
        
        return 
    
        
    def plot_corr_matrix(self):
        self.slip_corr_f()
        title = f'{self.name} Correlation Matrix'
        fig, ax = plt.subplots(dpi=600)
        ax.imshow(self.slip_corr)
        ax.set_title(title,fontweight='bold',fontsize=11)
        ax.set_aspect('equal', 'box')
        ax.tick_params(labelsize=12)
        plt.tight_layout()
        
        parameter_dir = os.path.join(self.output_plots_model_dir,'correlation_matrix')
        try: 
            os.makedirs(parameter_dir)  
        except FileExistsError:
            pass
        
        plot_path = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_beta_{self.step}_slip_correlation_matrix.jpg')
          
        fig.savefig(plot_path)
        
    
    
    def plot_corr(self):
        self.max_patch()
        self.slip_corr_f()
        
        fig, ax = plt.subplots(dpi=600)
        # correlation
    
        parameter = self.slip_corr[:,self.max_slip]

        title = f"{self.name} correlation for max slip patch in slip"

        # im0 = ax.pcolormesh(self.x,self.ydown,parameter.reshape(self.nrows,self.ncols),
        #                           edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
        im0 = ax.pcolormesh(self.x,self.ydown,np.flip(parameter.reshape(self.nrows,self.ncols,order='F'),axis=0),
                                  edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
        # im0 = ax.pcolormesh(self.x,self.ydown,np.flip(parameter.reshape(self.nrows,self.ncols,order='F'),axis=0),
        #                           edgecolors='k', cmap='bwr')
        ax.plot(self.x[self.ncol_target_flipped],self.ydown[self.nrow_target_flipped],marker='o',color='black',ms=3)

        # ax.plot(self.x[self.ncol_target],self.ydown[self.nrow_target],marker='o',color='black',ms=3)
        #X0, Y0 = np.meshgrid(x0, y0)
        #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        
        try:
            hypocenter  = (self.mean_d['Hypo_as'],-self.mean_d['Hypo_dd'])
            ax.plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
        except:
            pass 
        ax.margins(0)        
        fig.colorbar(im0, ax=ax,shrink=0.9,label='correlation')
        ax.set_ylabel('Down-dip distance (km)',fontsize=12)
        ax.set_xlabel('Along-strike distance (km)',fontsize=12)
        ax.set_aspect('equal', 'box')
        ax.set_title(title,fontweight='bold',fontsize=11)
        ax.tick_params(labelsize=12)
        plt.tight_layout()
        
        parameter_dir = os.path.join(self.output_plots_model_dir,'correlation_max')
        try: 
            os.makedirs(parameter_dir)  
        except FileExistsError:
            pass
        
        plot_path = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_beta_{self.step}_max_slip_correlation.jpg')
          
        fig.savefig(plot_path)
        
    
    
    
    
    

    

                
    def standard_deviation(self):
        
        Np = self.Np
        nramp = self.nramp
        # standard deviation of each parameter at each node 
        std = np.std(self.sampled_data,axis=1)
        U = np.sqrt(self.sampled_data[:Np,:]**2+self.sampled_data[Np:2*Np,:]**2)
        stdU = np.std(U,axis=1)
        std0 = np.copy(std)
        if self.name =='Gorkha':
              std[:Np],std[Np:2*Np] = std0[Np:2*Np],std0[:Np]
        
        if self.name == 'Pedernales':
            std_ss = np.zeros(Np)
            std_ds = np.zeros(Np)
            strike = self.geometry_d['strike']*((np.pi)/180)
            dip = self.geometry_d['dip']*((np.pi)/180)
            for p in range(Np):
                RotAngle2 = self.RotAngle*((np.pi) / 180.)
                
                rotation = np.arctan2(np.tan(strike[p]) - np.tan(RotAngle2), 
                                    np.cos(dip[p])*(1.+np.tan(RotAngle2)*np.tan(strike[p])))
    
                # If RotAngle within ]90, 270], change rotation
                if self.RotAngle > 90. and self.RotAngle<=270.:
                    rotation += np.pi
    
                rp = std[p]    # rake-perpendicular
                ar = std[p+Np] # rake-parallel
    
                std_ss[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                std_ds[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
            std_Uperp = std_ss
            std_Uparallel = std_ds
        else:
            rake0 = (self.rake - 90)*(np.pi/180)
            std_Uperp =  std[:Np]*np.cos(rake0)  - std[Np:2*Np]*np.sin(rake0) 
            std_Uparallel = std[:Np]*np.sin(rake0)  + std[Np:2*Np]*np.cos(rake0)     
        # skewness of each parameter at each node
    
        std_Uperp = np.abs(std_Uperp)
        std_Uparallel = np.abs(std_Uparallel)
        if self.model_type=='static':

                self.std_d = {
                    'std_U_perp':std_Uperp ,# std[:Np],         # along-strike
                    'std_U_parallel':std_Uparallel, # std[Np:2*Np], # along-dip
                    'std_U':stdU
                       }
        else:

            self.std_d = {
                   'std_U_perp':std_Uperp, #std[:Np],         # along-strike
                   'std_U_parallel':std_Uparallel, #std[Np:2*Np], # along-dip
                   'std_U':stdU,
                   'std_Tr':std[2*Np+nramp:3*Np+nramp],
                   'std_Vr':std[3*Np+nramp:4*Np+nramp],
                   'std_Hypo_as': std[4*Np+nramp:4*Np+nramp+1], # along-strike
                   'std_Hypo_dd': std[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                           }

                    
    def skewness(self):
          Np = self.Np
          nramp = self.nramp
          # skewness of each parameter at each node
          
          skewness = skew(self.sampled_data,axis=1)
          skewness0 = np.copy(skewness)
          if self.name =='Gorkha':
                skewness[:Np],skewness[Np:2*Np] = skewness0[Np:2*Np],skewness0[:Np]
          if self.model_type=='static':

                  self.skew_d = {
                      'skew_U_perp':skewness[:Np],
                      'skew_U_parallel':skewness[Np:2*Np]
                      }

          else:

              self.skew_d = {
                      'skew_U_perp':skewness[:Np],
                      'skew_U_parallel':skewness[Np:2*Np],
                      'skew_Tr':skewness[2*Np+nramp:3*Np+nramp],
                      'skew_Vr':skewness[3*Np+nramp:4*Np+nramp],
                      'skew_Hypo_as':skewness[4*Np+nramp:4*Np+nramp+1], # along-strike
                      'skew_Hypo_dd':skewness[4*Np+nramp+1:4*Np+nramp+2] # down-dip
                      }
                                             
    
    def model_save(self):
        full_model = self.mean_d | self.std_d 
        full_model = full_model | self.geometry_d 
        if skew==True:
            full_model = full_model | self.skew_d 
        # Geometry file for Pedernales is not definite; requires further reformatting
            
        
    #  PEDERNALES nramp = 9, GORKHA nramp = 0?
    # if name=='Pedernales':
    #     nramp = 9
        # If RotAngle within ]90, 270], change rotation
        # if RotAngle > 90. and RotAngle<=270.:
        #     rotation += np.pi

        # rp = alt.parameter('strike',samples[:,p]   ,bins=100,dtype=np.float32,color='red')
        # ar = alt.parameter('dip'   ,samples[:,p+Np],bins=100,dtype=np.float32,color='red')

        # ss = copy.deepcopy(rp)
        # ds = copy.deepcopy(ar)
        # ss.values = ar.values*np.cos(rotation) - rp.values*np.sin(rotation)
        # ds.values = ar.values*np.sin(rotation) + rp.values*np.cos(rotation)
        
    
        try: 
            os.makedirs(self.output_mean_model_dir)    
        except FileExistsError:
           pass   
        
        full_model_csv_file_dir = os.path.join(self.output_mean_model_dir,f'{self.name}_mean_{self.model_type}_model.csv')
        full_model_df = pd.DataFrame.from_dict(full_model, orient='index') # handle `Value error exception` 
        full_model_df = full_model_df.transpose()
        full_model_df.to_csv(full_model_csv_file_dir)
        
    



    def plot_distribution(self,parameter_type:str,parameter,figure_s,color_bar_label,padding):
        nrows,ncols = self.nrows, self.ncols
        X,Y = self.X, self.Y
        if parameter_type == 'mean':
            dictionary = self.mean_d
        elif parameter_type == 'std':
            dictionary = self.std_d
        elif parameter_type == 'skew':
            dictionary = self.skew_d
        
        fig, ax = plt.subplots(figsize=figure_s,dpi=600)
        param = dictionary[parameter].reshape(ncols,nrows).transpose()
        
        try:
            hypocenter = (dictionary['Hypo_as'],dictionary['Hypo_dd'])
            ax.plot(hypocenter[0],hypocenter[1],marker='*',color='blue',markeredgecolor='black',markersize=25)
        except:
            pass
            #hypocenter = (145,35) # for Gorkha
            #ax.plot(hypocenter[0],hypocenter[1],marker='*',color='blue',markeredgecolor='black',markersize=25)
        if parameter_type=='skew':
            minimum = min(param.min(),-param.max())
            maximum = max(-param.min(),param.max())
            im = ax.pcolormesh(self.x,self.y,param,edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=minimum, vmax=maximum))
        
        else:
            im = ax.pcolormesh(self.x,self.y,param,edgecolors='k',cmap='YlOrRd')
    
        #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        ax.margins(0)
        if parameter=='Slip':
            U = dictionary['U_parallel'].reshape(ncols,nrows).transpose() # dip-slip displacement
            V = dictionary['U_perp'].reshape(ncols,nrows).transpose()     # strike-slip displacement (see p.4 Minson et al. (2013) Part II)
            q = ax.quiver(X,Y,V,U,scale=0.9,scale_units ='x', units='width',width=0.004,headwidth=4.2,headlength=6.0)
            ax.quiverkey(q, X=0.11, Y=1.075, U=round(param.max(),-1),
                          label=' {} m'.format(round(param.max(),-1)), labelpos='E')
            factor = 1/q.scale    # scaling factor from UV unit to XY units
            offsetXY = np.stack((V,U),axis=2).flatten().reshape(nrows*ncols,2)   
            
            # minus sign in Y coordinate corrects for y-axis inversion
            
            offsetXY = np.array([[xy[0],-xy[1]] for xy in offsetXY])  
            coord = q.XY + offsetXY*factor
            # width and height had to be formatted in row-like order (initially in column-like order)
            ells = [Ellipse(xy=(coord[i][0],coord[i][1]),
                            width=self.std_d['std_U_perp'].reshape((nrows,ncols),order='F').flatten()[i]*factor, 
                            height=self.std_d['std_U_parallel'].reshape((nrows,ncols),order='F').flatten()[i]*factor,
                            angle=0,alpha=0.5,fill=False,edgecolor='k')
                    for i in range(q.N)]
            for e in ells:
                ax.add_artist(e)  
                
        fig.colorbar(im, ax=ax,shrink=0.33,label='{}'.format(color_bar_label))
        im.figure.axes[1].tick_params(axis='both', labelsize=10) 
        
        ax.set_ylabel('Down-dip distance (km)',fontsize=10)
        ax.set_xlabel('Along-strike distance (km)',fontsize=10)
        ax.set_aspect('equal', 'box')
        ax.tick_params(labelsize=10)
        ax.set_title(f'{self.name}_{self.model_type}_model'.replace("_"," ") + r' $ \beta $=%d'%(self.step) ,pad=padding,fontweight='bold',fontsize=10)
        ax.invert_yaxis()
        
        parameter_dir = os.path.join(self.output_plots_model_dir,f'{parameter_type}')
        try: 
            os.makedirs(parameter_dir)  
        except FileExistsError:
            pass
        
        plot_path = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_beta_{self.step}_{parameter_type}_{parameter}.jpg')
        
        fig.savefig(plot_path) 
    
    