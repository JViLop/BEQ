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
                   nsamples=20000,step = 44):
        

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
        self.X, self.Y = np.meshgrid(self.x, self.ydown)
        
        
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
        if self.name =='Tohoku':
            self.model_file_step = os.path.join(self.model_type_dir,'step_'+ str(step))
            self.model_file_name = os.listdir(self.model_file_step)[0]
            self.model_file_dir = os.path.join(self.model_file_step,self.model_file_name)
        
        
        else:
            
            self.model_file_name = os.listdir(self.model_type_dir)[0]
            self.model_file_dir = os.path.join(self.model_type_dir,self.model_file_name)
        
        
        # OUTPUT directories
        self.output_EQs_dir = os.path.join(os.getcwd(),'INPUT')
        self.output_EQ_dir =  os.path.join(self.output_EQs_dir,self.name)
        self.output_model_dir = os.path.join(self.output_EQ_dir,'model')
        self.output_model_type_dir = os.path.join(self.output_model_dir,self.model_type)
        if sampling:
            self.output_sampling_type = os.path.join(self.output_model_type_dir,f'{self.nsamples}_samples')
            self.output_mean_model_dir = os.path.join(self.output_sampling_type,'mean')
            self.output_plots_model_dir = os.path.join(self.output_sampling_type ,'figures')
            self.output_bin_data_dir = os.path.join(self.output_sampling_type,'bin_data')
        else:
            self.output_sampling_type = os.path.join(self.output_model_type_dir,'all_samples')    
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
        Np = self.Np
        nramp = self.nramp
        if self.file_type=='h5':
            f = h5py.File(self.model_file_dir,'r')
            try:
                
                data = f['Sample Set']
                data = np.array(data).transpose()
                if self.name =='Pedernales':
                    zeros_id = np.where(~data.any(axis=0))[0]
                    data = np.delete(data,zeros_id,axis=1)
                    zero_id_tr = np.where(data[2*Np+nramp:3*Np+nramp,:]==0)[1]
                    data = np.delete(data,zero_id_tr,axis=1)
                    zero_id_vr = np.where(data[3*Np+nramp:4*Np+nramp,:]==0)[1]
                    data = np.delete(data,zero_id_vr,axis=1)
                    zero_id_hyp = np.where(data[4*Np+nramp+1:,:]==0)[1]
                    data = np.delete(data,zero_id_hyp,axis=1)
                    zero_id_hyp = np.where(data[4*Np+nramp:4*Np+nramp+1,:]==0)[1]
                    data = np.delete(data,zero_id_hyp,axis=1)
            
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
            print(f'{self.name} shape:',data.shape)
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
    
    
    def max_slip_patch(self):
        #return 
        # slip = np.flip(self.mean_d['Slip'].reshape(self.ncols,self.nrows).transpose(),axis=0).flatten()

        slip = self.mean_d['Slip']
        max_slip_patch = np.argmin(abs(slip-max(slip)))
        col0_max_slip,row0_max_slip = max_slip_patch%self.ncols,max_slip_patch//self.ncols
        self.ncol_target = col0_max_slip 
        self.nrow_target = row0_max_slip 
        max_observable_patch = self.nrow_target*self.ncols + self.ncol_target 
        self.max_slip = max_observable_patch # index of max_slip_patch if flatten
        
        slip_flipped = np.flip(self.mean_d['Slip'].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
        max_flipped_slip_patch = np.argmin(abs(slip_flipped-max(slip_flipped)))
        col_max_flipped_slip,row_max_flipped_slip = max_flipped_slip_patch%self.ncols,max_flipped_slip_patch//self.ncols
        self.max_slip = max_flipped_slip_patch # index of max_slip_patch if flatten
        self.ncol_target_flipped = col_max_flipped_slip 
        self.nrow_target_flipped = row_max_flipped_slip
        
    
        return 
       
    def order_id_slip(self):
        slip = self.mean_d['Slip']
        
        
    def fault_patch(self,patch_id):
        col_slip,row_slip = patch_id%self.ncols,patch_id//self.ncols

        return col_slip,row_slip
        
    def slip_corr_f(self):
            
        U1 = self.sampled_data[:self.Np,:]
        U2 = self.sampled_data[self.Np:2*self.Np,:] 
        slip_ensemble = np.sqrt(U1**2+U2**2)
        
        corr = np.corrcoef(slip_ensemble)
        self.slip_corr = corr
        
        
        return 
    
    def spatial_corr(self):
        U1 = self.sampled_data[:self.Np,:]
        U2 = self.sampled_data[self.Np:2*self.Np,:] 
        
        U_ensemble = np.sqrt(U1**2+U2**2)
        Tr_ensemble = self.sampled_data[2*self.Np + self.nramp:3*self.Np + self.nramp,:]
        Vr_ensemble = self.sampled_data[2*self.Np + self.nramp:3*self.Np + self.nramp,:]
        U_Tr_ensemble = U_ensemble/Tr_ensemble

        
        U_corr = np.corrcoef(U_ensemble)
        Tr_corr = np.corrcoef(Tr_ensemble)
        Vr_corr = np.corrcoef(Vr_ensemble)
        U_Tr_corr = np.corrcoef(U_Tr_ensemble)
        
        for i in range(self.Np):
            U_corr[i,:] = np.flip(U_corr[i,:].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            Tr_corr[i,:] = np.flip(Tr_corr[i,:].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            Vr_corr[i,:] = np.flip(Vr_corr[i,:].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            U_Tr_corr[i,:] = np.flip(U_Tr_corr[i,:].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            
        for j in range(self.Np):
            U_corr[:,j] = np.flip(U_corr[:,j].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            Tr_corr[:,j] = np.flip(Tr_corr[:,j].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            Vr_corr[:,j] = np.flip(Vr_corr[:,j].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
            U_Tr_corr[:,j] = np.flip(U_Tr_corr[:,j].reshape(self.nrows,self.ncols,order='F'),axis=0).flatten()
        
        
        self.U_corr = U_corr
        self.Tr_corr =Tr_corr
        self.Vr_corr = Vr_corr
        self.U_Tr_corr = U_Tr_corr
     
        return 
    @staticmethod
    def order_corr(M,indeces):
        corr_ordered = M[indeces,:]
        corr_ordered = corr_ordered[:,indeces]
        return corr_ordered
    
    def corr_matrix(self,order_by_distance=True):
        self.spatial_corr()
        self.max_slip_patch()
        patch_id = self.max_slip
        U_corr_ordered = np.zeros_like(self.U_corr)
        Tr_corr_ordered = np.zeros_like(self.Tr_corr)
        Vr_corr_ordered = np.zeros_like(self.Vr_corr)
        U_Tr_corr_ordered = np.zeros_like(self.U_Tr_corr)

        parameters_dir = os.path.join(self.output_plots_model_dir,'spatial_correlation')
        os.makedirs(parameters_dir,exist_ok=True)
        fig_dir = os.path.join(parameters_dir,f'{self.name}_{self.model_type}_{self.nsamples}_spatial_ordered_corr_matrix.pdf')
        Slip = self.mean_d['Slip']
        Vr = self.mean_d['Vr']
        Tr = self.mean_d['Tr']
        U_Tr = Slip/Tr
        Slip_id_sorted = np.argsort(np.flip(Slip.reshape(self.nrows,self.ncols,order='F'),axis=0).flatten())
        Tr_id_sorted = np.argsort(np.flip(Tr.reshape(self.nrows,self.ncols,order='F'),axis=0).flatten())
        Vr_id_sorted = np.argsort(np.flip(Vr.reshape(self.nrows,self.ncols,order='F'),axis=0).flatten())
        U_Tr_id_sorted = np.argsort(np.flip(U_Tr.reshape(self.nrows,self.ncols,order='F'),axis=0).flatten())
        if order_by_distance:
            xflat = self.X.flatten()
            yflat = self.Y.flatten()
            for n in range(self.Np):
                x_target = xflat[n]
                y_target = yflat[n]
                r_target = np.array([x_target,y_target])
                r_all = np.column_stack((xflat,yflat))
                dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2)
                dr_sorted = np.sort(dr)
                ind_dr_sorted = np.argsort(dr)
                
                U_corr_ordered[n,:] = self.U_corr[n,ind_dr_sorted]
                Tr_corr_ordered[n,:] = self.Tr_corr[n,ind_dr_sorted]
                Vr_corr_ordered[n,:] = self.Vr_corr[n,ind_dr_sorted]
                U_Tr_corr_ordered[n,:] = self.U_Tr_corr[n,ind_dr_sorted]
            
            for n in range(self.Np):
                U_corr_ordered[:,n] = U_corr_ordered[Slip_id_sorted,n] 
                Tr_corr_ordered[:,n] = Tr_corr_ordered[Tr_id_sorted,n] 
                Vr_corr_ordered[:,n] = Vr_corr_ordered[Vr_id_sorted,n] 
                U_Tr_corr_ordered[:,n] = U_Tr_corr_ordered[U_Tr_id_sorted,n] 
            
            CORR = {'U':U_corr_ordered,
                    'Tr': Tr_corr_ordered,
                    'Vr': Vr_corr_ordered,
                    'U_Tr':U_Tr_corr_ordered}
            
            labels = {'U':'$U$ (m)',
                          'Tr':'$T_r$ (s)',
                          'Vr':'$V_r$ (km/s)',
                          'U_Tr':'$U/T_r$ (m/s)'}
            
            fig,axes = plt.subplots(2,2,figsize=(10,10))
            for k,param in enumerate(list(CORR.keys())):
                im = axes[k//2][k%2].imshow(CORR[param],cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax = 1))
                axes[k//2][k%2].axes.get_xaxis().set_ticks([])
                axes[k//2][k%2].axes.get_yaxis().set_ticks([])
                props = dict(boxstyle='round', facecolor='white',alpha=0.5)
                axes[k//2][k%2].text(0.78,0.91,f'{labels[param]}',fontsize=10,transform=axes[k//2][k%2].transAxes,bbox=props)
                cbar = fig.colorbar(im,ax=axes[k//2][k%2],shrink=0.6,label='correlation')
                cbar.set_ticks([-1,-0.5,0,0.5,1])
                cbar.ax.tick_params(labelsize=8)
                axes[k//2][k%2].set_aspect('equal', 'box')
                # Add arrow annotations for x-axis
                #axes[k//2][k%2].annotate('', xy=(0, 80), xytext=(0, 0), arrowprops=dict(facecolor='black',arrowstyle='->'))
                axes[k//2][k%2].text(-0.085, 0.7, r'$\longleftarrow$' + ' Increasing ' + r'{}'.format(labels[param]) , va='center', rotation=90,transform=axes[k//2][k%2].transAxes) 
                #axes[k//2][k%2].annotate('', xy=(80, 0.), xytext=(0, 0.), arrowprops=dict(facecolor='black',arrowstyle='->'))
                axes[k//2][k%2].text(0.5, 1.025, 'Increasing patch-to-patch distance ' + r'$\longrightarrow$', ha='center',transform=axes[k//2][k%2].transAxes)
            plt.tight_layout()
            fig.savefig(fig_dir)
            
            
            
    def spatial_corr_max(self,patch_id = None):
        self.spatial_corr()
        
        if patch_id == None:
            self.max_slip_patch()
            patch_id = self.max_slip
            ncol_target = self.ncol_target_flipped
            nrow_target =  self.nrow_target_flipped
        else:
            ncol_target,nrow_target = self.fault_patch(patch_id)
        
        parameters = {'U':self.U_corr[:,patch_id],
                      'Tr':self.Tr_corr[:,patch_id],
                      'Vr':self.Vr_corr[:,patch_id],
                      'U_Tr':self.U_Tr_corr[:,patch_id]}
        
        labels = {'U':'$U$ (m)',
                      'Tr':'$T_r$ (s)',
                      'Vr':'$V_r$ (km/s)',
                      'U_Tr':'$U/T_r$ (m/s)'}
        
        xflat = self.X.flatten()
        yflat = self.Y.flatten()
        x_target = xflat[patch_id]
        y_target = yflat[patch_id]
        r_target = np.array([x_target,y_target])
        r_all = np.column_stack((xflat,yflat))
        dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2)
        dr_sorted = np.sort(dr)
        ind_dr_sorted = np.argsort(dr)
        self.scaled_dr = (dr_sorted - np.min(dr_sorted))/(np.max(dr_sorted) - np.min(dr_sorted))
        self.corr_wrt_max = {'U':parameters['U'][ind_dr_sorted],
                      'Tr':parameters['Tr'][ind_dr_sorted],
                      'Vr':parameters['Vr'][ind_dr_sorted],
                      'U_Tr':parameters['U_Tr'][ind_dr_sorted]}
        
        
        
        
    def plot_spatial_corr(self,spacing,patch_id=None):

        self.spatial_corr()
        if patch_id == None:
            self.max_slip_patch()
            patch_id = self.max_slip
            ncol_target = self.ncol_target_flipped
            nrow_target =  self.nrow_target_flipped
        else:
            ncol_target,nrow_target = self.fault_patch(patch_id)
            
        fig,axes = plt.subplots(2,2,figsize = (6,6))
        parameters = {'U':self.U_corr[:,patch_id],
                      'Tr':self.Tr_corr[:,patch_id],
                      'Vr':self.Vr_corr[:,patch_id],
                      'U_Tr':self.U_Tr_corr[:,patch_id]}
        
        labels = {'U':'$U$ (m)',
                      'Tr':'$T_r$ (s)',
                      'Vr':'$V_r$ (km/s)',
                      'U_Tr':'$U/T_r$ (m/s)'}
        
        
        for k,param in enumerate(list(parameters.keys())):
            
            parameter = parameters[param]
            parameters_dir = os.path.join(self.output_plots_model_dir,'spatial_correlation')
            os.makedirs(parameters_dir,exist_ok=True)
            parameter_dir = os.path.join(parameters_dir,f'{param}')
            os.makedirs(parameter_dir,exist_ok=True)
            file_dir = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_{self.nsamples}_spatial_corr_{param}_wrt_max_slip.txt')
            np.savetxt(file_dir,parameter)
            title = f"{labels[param]}"
           
                
            # im0 = ax.pcolormesh(self.x,self.ydown,parameter.reshape(self.nrows,self.ncols),
            #                           edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
            #im0 = axes[k//2][k%2].pcolormesh(self.x,self.ydown,np.flip(parameter.reshape(self.nrows,self.ncols,order='F'),axis=0),
            # edgecolors='k', lw=0.2,cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
            im0 = axes[k//2][k%2].pcolormesh(self.x,self.ydown,parameter.reshape(self.nrows,self.ncols),
                                      edgecolors='k', lw=0.2,cmap='bwr',norm=TwoSlopeNorm(0,vmin=-1,vmax=1))
            # im0 = ax.pcolormesh(self.x,self.ydown,np.flip(parameter.reshape(self.nrows,self.ncols,order='F'),axis=0),
            #                           edgecolors='k', cmap='bwr')
            axes[k//2][k%2].plot(self.x[ncol_target],self.ydown[nrow_target],marker='o',color='black',ms=2)
    
            # ax.plot(self.x[self.ncol_target],self.ydown[self.nrow_target],marker='o',color='black',ms=3)
            #X0, Y0 = np.meshgrid(x0, y0)
            #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
            
            try:
                hypocenter  = (self.mean_d['Hypo_as'],-self.mean_d['Hypo_dd'])
                axes[k//2][k%2].plot(hypocenter[0],hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=10,markeredgewidth=0.55)
            except:
                pass 
            axes[k//2][k%2].margins(0)        
            fig.colorbar(im0, ax=axes[k//2][k%2],shrink=0.3,label='Correlation')
            if k==0:
                axes[k//2][k%2].set_ylabel('Down-dip distance (km)',fontsize=9)
                axes[k//2][k%2].set_xlabel('Along-strike distance (km)',fontsize=9)
            axes[k//2][k%2].set_aspect('equal', 'box')
            axes[k//2][k%2].set_title(title,fontsize=12)
            axes[k//2][k%2].tick_params(labelsize=9)
        plt.subplots_adjust(wspace=0.375,hspace=spacing)
            
        fig_dir = os.path.join(parameters_dir,f'{self.name}_{self.model_type}_{self.nsamples}_spatial_corr_wrt_max_slip.pdf')
          
        fig.savefig(fig_dir)
        
        
            
        
        
        
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
        
        plot_path = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_{self.nsamples}_slip_correlation_matrix.jpg')
          
        fig.savefig(plot_path)
        
    
    
    def plot_corr(self):
        self.max_slip_patch()
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
        
        plot_path = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_max_slip_correlation.jpg')
          
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
        full_model = {**self.mean_d, **self.std_d} 
        full_model = {**full_model,**self.geometry_d}
        if skew==True:
            full_model = {**full_model, **self.skew_d}
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
        X,Y = np.meshgrid(self.x,self.y)
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
        ax.set_title(f'{self.name}_{self.model_type}_model'.replace("_"," "),pad=padding,fontweight='bold')
        ax.invert_yaxis()
        
        parameter_dir = os.path.join(self.output_plots_model_dir,f'{parameter_type}')
        try: 
            os.makedirs(parameter_dir)  
        except FileExistsError:
            pass
        
        plot_path = os.path.join(parameter_dir,f'{self.name}_{self.model_type}_{parameter_type}_{parameter}.jpg')
        
        fig.savefig(plot_path) 
    
    
