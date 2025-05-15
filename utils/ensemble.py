import sys
import copy 

from csi.okadafull import displacement, stress, strain
from .observable import Observable
import pandas as pd
import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


class Ensemble_Displacement(Observable):
    def __init__(self,
                        parent_dir,
                        name,
                        model_type,
                        step,
                        want_step,
                        patchsize,
                        shape_geom,
                        observable = 'displacement',
                        samples=100,
                        new_station=False,
                        xstation = None,
                        ystation = None,
                        strikeOkada = 90,
                        factor = 4,
                        nparameters = 533,
                        RotAngle = None,
                        rake=90,
                        mu=30.0*1e9,
                        nu=0.25):
            
            self.nparameters = nparameters
            self.new_station = new_station
            self.RotAngle = RotAngle
            self.rake = rake
            self.mu = mu 
            self.nu = nu
            
            super().__init__(
                     parent_dir,
                     name,
                     model_type,
                     step,
                     want_step,
                     patchsize,
                     shape_geom,
                     observable = observable,
                     samples=samples,
                     new_station=new_station,
                     xstation = xstation ,
                     ystation = ystation,
                     strikeOkada = strikeOkada,
                     factor = factor)
            
            self.Multiple_Okada_displacement()
            self.save_statistics()
            self.plot_statistics()

    def make_dirs(self,file='data'):
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir=file)
        os.makedirs(folder_dir,exist_ok=True)
        return folder_dir

    def bin_data_reader(self):
        file_folder = self.dir_finder(len(self.list_in)-1,tail_dir='bin_data')
        file_dir = os.path.join(file_folder,f'{self.name}_{self.model_type}_n_{self.samples}.dat')
        self.bin_data = np.fromfile(file_dir,'float').reshape((self.nparameters,self.samples))
    
    def modify_stn(self):
        if self.new_station:
            self.zstn_flat = np.zeros(len(self.xstn)*len(self.ystn))

    def Multiple_Okada_displacement(self):
        self.bin_data_reader()
        self.mean_model_reader()
        self.source()
        self.station()
        self.modify_stn()
        
        
        h5f_name =  f'{self.name}_{self.observable}_nsamples_{self.samples}.h5'
        h5f_folder = self.make_dirs(file='data_ensemble')
        self.h5file_dir = os.path.join(h5f_folder,h5f_name)   
        h5file = h5py.File(self.h5file_dir,'w')
        
        displacement_dset = h5file.create_dataset(f'{self.observable}',shape=(self.samples,len(self.xstn)*len(self.ystn)*3))
        
        xstation = h5file.create_dataset('xstn',shape = (1,len(self.xstn)))
        ystation = h5file.create_dataset('ystn',shape = (1,len(self.ystn))) 
        xstation[:],ystation[:] = self.xstn, self.ystn 

        xsource = h5file.create_dataset('xsrc',shape = (1,len(self.xsrc)))
        ysource = h5file.create_dataset('ysrc',shape = (1,len(self.ysrc))) 
        xsource[:],ysource[:] = self.xsrc, self.ysrc

        for i in range(self.bin_data.shape[1]):

            Np = self.ncols*self.nrows
            model = self.bin_data[:,i] 
            model0 = np.copy(model)
            if self.name =='Gorkha':
                model[:Np],model[Np:2*Np] = model0[Np:2*Np],model0[:Np]
            
            if self.name=='Pedernales':
                ss0 = np.zeros(Np)
                ds0 = np.zeros(Np)
                for p in range(Np):
                    RotAngle2 = self.RotAngle*(np.pi/180.)
                    
                    rotation = np.arctan2(np.tan(self.strike[p]*(np.pi/180)) - np.tan(RotAngle2), 
                                        np.cos(self.dip_rad[p])*(1.+np.tan(RotAngle2)*np.tan(self.strike[p]*(np.pi/180))))

                    # If RotAngle within ]90, 270], change rotation
                    if self.RotAngle > 90. and self.RotAngle<=270.:
                        rotation += np.pi

                    rp = model[p]    # rake-perpendicular
                    ar = model[p+Np] # rake-parallel

                    ss0[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                    ds0[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
                Uperp = ss0
                Uparallel = ds0       
            else:
                rake0 = (self.rake - 90)*(np.pi/180)
                Uperp =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0) 
                Uparallel = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)      

            ss = self.array_formatter(Uperp)
            ds = self.array_formatter(Uparallel)
            ts = np.zeros(self.xsrc_flat.shape)
            Displacement = displacement(
                                self.xstn_flat,
                                self.ystn_flat,
                                self.zstn_flat, 
                                self.xsrc_flat, 
                                self.ysrc_flat, 
                                self.zsrc_flat, 
                                self.width, 
                                self.length, 
                                self.strike_rad, 
                                self.dip_rad, 
                                ss, 
                                ds, 
                                ts, nu=self.nu)
            
            
            displacement_dset[i,:] = Displacement.flatten(order='F')
         
        h5file.close()
    
    def read_h5file(self):
        # self.Multiple_Okada_displacement()

        f = h5py.File(self.h5file_dir,'r')
        dset = np.array(f[f'{self.observable}'])
        return dset
    
    def cov(self):
        dset  = self.read_h5file()
        self.covariance = np.cov(dset.transpose())
        
        nparameters = dset.shape[1]
        
        self.cov1 = self.covariance[:nparameters//3,:nparameters//3]
        self.cov2 = self.covariance[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
        self.cov3 = self.covariance[2*nparameters//3:,2*nparameters//3:]

        # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
        # standard deviation (= square root of variance)

        self.std1 = np.sqrt(self.cov1.diagonal())
        self.std2 = np.sqrt(self.cov2.diagonal())
        self.std3 = np.sqrt(self.cov3.diagonal())

    def corr(self):
        dset  = self.read_h5file()
        self.correlation = np.corrcoef(dset.transpose())
        
        nparameters = dset.shape[1]
        self.corr1 = self.correlation[:nparameters//3,:nparameters//3]
        self.corr2 = self.correlation[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
        self.corr3 = self.correlation[2*nparameters//3:,2*nparameters//3:]

    
    def plot_cov(self):
        self.cov()
        self.Okada_input()
        ncol_target,nrow_target,patchid = self.find_slip() 

        cov = {'Along-strike':self.std1,
                'Trench-normal':self.std2,
                'Vertical':self.std3}

        x,y = self.xstn*1e-3,self.ystn*1e-3   
        x0,y0 = self.xsrc*1e-3,self.ysrc*1e-3
        
        nrows = len(y)
        ncols = len(x)

        n0rows = len(y0)
        n0cols = len(x0) 

        # fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=900)
        for i,parameter_id in enumerate(cov.keys()):
            parameter = cov[parameter_id]

            title = f"{self.name} Uncertainty in {parameter_id} Surface Displacement "

            im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow',linewidth=0.25)
            parray = np.zeros((n0rows,n0cols))
            masked_arr = np.ma.masked_array(parray,parray<1)
            # axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=0.5)
            axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1.0)
            axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=3.0)
            # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
            #X0, Y0 = np.meshgrid(x0, y0)
            #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
            try:
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
            except:
                pass 
            axes[i].margins(0)        
            fig.colorbar(im0, ax=axes[i],shrink=0.65,label='Uncertainty (m)')
            axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
            axes[i].set_aspect('equal', 'box')
            axes[i].set_title(title,fontweight='bold')
            axes[i].tick_params(labelsize=12)
        
        
        # fig, axes = plt.subplots(1,3,figsize=(13,4.5),dpi=900)

        #### remember to shift work as trench is located far away ###
        # for i,parameter_id in enumerate(cov.keys()):
        #     parameter = cov[parameter_id]

        #     title = f"{self.name} Uncertainty in \n {parameter_id} surface {self.observable} "

        #     im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow')
        #     parray = np.zeros((n0rows,n0cols))
        #     masked_arr = np.ma.masked_array(parray,parray<1)
        #     axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
        #     axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1.5)
        #     axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=4)
        #     # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #     #X0, Y0 = np.meshgrid(x0, y0)
        #     #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        #     try:
        #         axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=10)
        #     except:
        #         pass 
        #     axes[i].margins(0)        
        #     cbar = fig.colorbar(im0, ax=axes[i],shrink=0.70,orientation='horizontal',label='Uncertainty (m)')
        #     axes[i].set_ylabel('Trench-normal distance (km)',fontsize=8)
        #     axes[i].set_xlabel('Along-strike distance (km)',fontsize=8)
        #     axes[i].set_aspect('equal', 'box')
        #     axes[i].set_title(title,fontsize=9.5,fontweight='bold')
        #     axes[i].tick_params(labelsize=8)
        #     cbar.ax.tick_params(labelsize=8)
        plt.tight_layout()
    
        fname =  f'std_{self.name}_{self.observable}_nsamples_{self.samples}.jpg'
        file_folder = self.make_dirs(file=f'figures_std')
        file_dir = os.path.join(file_folder,fname)   
        fig.savefig(file_dir)
        plt.close()

    def plot_corr(self):
        self.corr()
        self.Okada_input()
        # if patchID:
        #     ncol_target,nrow_target,patchid = self.find_slip(maxpatch = False,patchId=patchID)
        #     extra = str(patchid) + ' patch '
        # else:
        #     ncol_target,nrow_target,patchid = self.find_slip() 
        ncol_target,nrow_target,patchid = self.find_slip() 
        extra = 'max. slip patch'
        
        
        corr = {'along-strike':self.corr1[:,patchid],
                'along-dip':self.corr2[:,patchid],
                'vertical':self.corr3[:,patchid]}

        x,y = self.xstn*1e-3,self.ystn*1e-3   
        x0,y0 = self.xsrc*1e-3,self.ysrc*1e-3

        nrows = len(y)
        ncols = len(x)
        n0rows = len(y0)
        n0cols = len(x0) 

        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        for i,parameter_id in enumerate(corr.keys()):
            parameter = corr[parameter_id]
            # title = f"{self.name}  Correlation in {extra} {parameter_id} surface {self.observable} "
            title = f"{self.name}  Correlation in {parameter_id} surface {self.observable} "

            im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),
                                 edgecolors='k', linewidth=0.25, cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
            parray = np.zeros((n0rows,n0cols))
            masked_arr = np.ma.masked_array(parray,parray<1)
            # axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=0.5)
            axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1.0)
            axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=3.0)
            # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
            #X0, Y0 = np.meshgrid(x0, y0)
            #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
            try:
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
            except:
                pass 
            axes[i].margins(0)        
            fig.colorbar(im0, ax=axes[i],shrink=0.65,label='Correlation')
            axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
            axes[i].set_aspect('equal', 'box')
            axes[i].set_title(title,fontweight='bold')
            axes[i].tick_params(labelsize=12)
        plt.tight_layout()
    
        fname =  f'corr_{extra}_{self.name}_{self.observable}_nsamples_{self.samples}.jpg'
        file_folder = self.make_dirs(file=f'figures_corr')
        file_dir = os.path.join(file_folder,fname)   
        fig.savefig(file_dir)
        plt.close()

    def save_cov(self):
        self.cov()
        fname =  f'covariance_{self.name}_{self.observable}_nsamples_{self.samples}.txt'
        file_folder = self.make_dirs(file='covariance')
        file_dir = os.path.join(file_folder,fname) 
        np.savetxt(file_dir,self.covariance)
        
    def save_corr(self):
        self.corr()
        fname =  f'correlation_{self.name}_{self.observable}_nsamples_{self.samples}.txt'
        file_folder = self.make_dirs(file='covariance')
        file_dir = os.path.join(file_folder,fname)
        np.savetxt(file_dir,self.correlation) 
    
    def plot_statistics(self):
        self.plot_corr()
        self.plot_cov()
    
    def save_statistics(self):
        self.save_corr()
        self.save_cov()

    




class Ensemble_Stress(Observable):
    def __init__(self,
                        parent_dir,
                        name,
                        model_type,
                        step,
                        want_step,
                        patchsize,
                        shape_geom,
                        observable = 'Stress change',
                        samples=100,
                        new_station=False,
                        xstation = None,
                        ystation = None,
                        strikeOkada = 90,
                        factor = 4,
                        nparameters = 533,
                        RotAngle = None,
                        rake=90,
                        mu=30.0*1e9,
                        nu=0.25):
            
            self.nparameters = nparameters
            self.new_station = new_station
            self.RotAngle = RotAngle
            self.rake = rake
            self.mu = mu 
            self.nu = nu
        
            
            super().__init__(
                     parent_dir,
                     name,
                     model_type,
                     step,
                     want_step,
                     patchsize,
                     shape_geom,
                     observable = observable,
                     samples=samples,
                     new_station=new_station,
                     xstation = xstation ,
                     ystation = ystation,
                     strikeOkada = strikeOkada,
                     factor = factor)
            
            self.Multiple_Okada_stress()
            self.save_statistics()
            self.plot_statistics()

    def make_dirs(self,file='data'):
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir=file)
        os.makedirs(folder_dir,exist_ok=True)
        return folder_dir

    def bin_data_reader(self):
        file_folder = self.dir_finder(len(self.list_in)-1,tail_dir='bin_data')
        file_dir = os.path.join(file_folder,f'{self.name}_{self.model_type}_n_{self.samples}.dat')
        self.bin_data = np.fromfile(file_dir,'float').reshape((self.nparameters,int(self.samples)))
    
    def modify_stn(self):
        if self.new_station:
                  
            self.zstn_flat = -self.depth.flatten()*1e3
            self.strike_stn_rad = self.strike_rad.flatten()
            self.dip_stn_rad = self.dip_rad.flatten()

        else:
            self.strike_stn_rad = self.strike_rad
            self.dip_stn_rad = self.dip_rad
            self.zstn_flat = -self.depth.flatten()*1e3



    def access_h5file(self):
        h5f_name =  f'{self.name}_{self.observable}_nsamples_{self.samples}.h5'
        h5f_folder = self.make_dirs(file='data_ensemble')
        self.h5file_dir = os.path.join(h5f_folder,h5f_name) 


    def stress_drop(self,stress_change,slip,weighted = True):
        stress_drop = -stress_change
        if weighted:
            average_stress = np.average(stress_drop,weights = slip)
        else:
            average_stress = np.mean(stress_drop)
        
        return average_stress

    def stress_drop_Noda(self,as_shear,ad_shear,as_slip,ad_slip):
        
        # need minus sign for stress drop

        delta_stress = -np.array([as_shear,ad_shear]).T
        delta_u = np.array([as_slip,ad_slip]).T
    
        prod_stress_u = np.sum(delta_stress*delta_u)
        sum_du = np.linalg.norm(np.sum(delta_u,axis=0))  
        
        try:
            average_stress = prod_stress_u/sum_du 
        except:
             average_stress = None
        return average_stress

    def set_h5file(self):
        h5f_name =  f'{self.name}_{self.observable}_nsamples_{self.samples}.h5'
        h5f_folder = self.make_dirs(file='data_ensemble')
        self.h5file_dir = os.path.join(h5f_folder,h5f_name) 

    def Multiple_Okada_stress(self,stress_drop=False,comp=2):
        self.bin_data_reader()
        self.mean_model_reader()
        self.set_h5file()
        self.source()
        self.station()
        self.modify_stn()
        
        if stress_drop:
            self.stress_drops = []

          
        h5file = h5py.File(self.h5file_dir,'w')
        
        stress_dset = h5file.create_dataset(f'{self.observable}',shape=(self.samples,len(self.xstn)*len(self.ystn)*3))
        
        xstation = h5file.create_dataset('xstn',shape = (1,len(self.xstn)))
        ystation = h5file.create_dataset('ystn',shape = (1,len(self.ystn))) 
        xstation[:],ystation[:] = self.xstn, self.ystn 
        
        xsource = h5file.create_dataset('xsrc',shape = (1,len(self.xsrc)))
        ysource = h5file.create_dataset('ysrc',shape = (1,len(self.ysrc))) 
        xsource[:],ysource[:] = self.xsrc, self.ysrc

        for i in range(self.bin_data.shape[1]):

            Np = self.ncols*self.nrows
            model = self.bin_data[:,i] 
            model0 = np.copy(model)
            if self.name =='Gorkha':
                model[:Np],model[Np:2*Np] = model0[Np:2*Np],model0[:Np]
            
            if self.name=='Pedernales':
                ss0 = np.zeros(Np)
                ds0 = np.zeros(Np)
                for p in range(Np):
                    RotAngle2 = self.RotAngle*(np.pi/180.)
                    
                    rotation = np.arctan2(np.tan(self.strike[p]*(np.pi/180)) - np.tan(RotAngle2), 
                                        np.cos(self.dip_rad[p])*(1.+np.tan(RotAngle2)*np.tan(self.strike[p]*(np.pi/180))))

                    # If RotAngle within ]90, 270], change rotation
                    if self.RotAngle > 90. and self.RotAngle<=270.:
                        rotation += np.pi

                    rp = model[p]    # rake-perpendicular
                    ar = model[p+Np] # rake-parallel

                    ss0[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                    ds0[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
                Uperp = ss0
                Uparallel = ds0       
            else:
                rake0 = (self.rake - 90)*(np.pi/180)
                Uperp =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0) 
                Uparallel = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)      

            ss = self.array_formatter(Uperp)
            ds = self.array_formatter(Uparallel)
            ts = np.zeros(self.xsrc_flat.shape)
            Stress = stress(
                                self.xstn_flat,
                                self.ystn_flat,
                                self.zstn_flat, 
                                self.xsrc_flat, 
                                self.ysrc_flat, 
                                self.zsrc_flat, 
                                self.width, 
                                self.length, 
                                self.strike_rad, 
                                self.dip_rad, 
                                ss, 
                                ds, 
                                ts, nu=self.nu)[0]


            # new
            npatches = int(len(self.xstn)*len(self.ystn))  
                         
            tnn = np.zeros(npatches)
            ts1 = np.zeros(npatches)
            ts2 = np.zeros(npatches)

            Sxx = Stress[:,0] 
            Sxy = Stress[:,1] 
            Sxz = Stress[:,2] 
            Syy = Stress[:,3] 
            Syz = Stress[:,4]
            Szz = Stress[:,5]

            for k in range(npatches):

                sxx = Sxx[k]*1e-6
                sxy = Sxy[k]*1e-6
                sxz = Sxz[k]*1e-6
                syy = Syy[k]*1e-6
                syz = Syz[k]*1e-6
                szz = Szz[k]*1e-6

                tau = [[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]]
                
                CS, SS = np.cos(self.strike_stn_rad[k]), np.sin(self.strike_stn_rad[k])
                CD, SD = np.cos(self.dip_stn_rad[k]), np.sin(self.dip_stn_rad[k])
            
                # set direction vectors
                NVec = np.array([SD*CS, -SS*SD, CD]) # fault-normal
                Ss = np.array([SS, CS, 0]) # along-strike before it was  -CS
                Sd = np.cross(NVec,Ss) # see Jolivet stressfield.py strikedip2normal() ; manually np.array([-CD*CS, CD*SS, SD]) # along-dip
                
                # compute traction vector
                trn = np.matmul(tau, NVec) 
                tnn[k] = np.dot(trn, NVec) # fault-normal
                ts1[k] = np.dot(trn, Ss) # along-strike 
                ts2[k] = np.dot(trn, Sd) # along-dip
            
            
            Stress  = np.array([tnn,ts1,ts2])
            if stress_drop:
                ## work well; commented on 23-03
                # average_stress_drop = self.stress_drop(ts2,ds,weighted=False)
                # self.stress_drops.append(average_stress_drop)
                
                # # Noda formulation
                try:
                    average_stress_drop = self.stress_drop_Noda(ts1,ts2,ss,ds)
                    self.stress_drops.append(average_stress_drop)
                except:
                    continue
            
            Stress = Stress.flatten()

            stress_dset[i,:] = Stress 
        
        h5file.close()
        if stress_drop:
            self.plot_histogram(self.stress_drops)


    



    def read_h5file(self):
        
        f = h5py.File(self.h5file_dir,'r')
        dset = np.array(f[f'{self.observable}'])
        return dset
    
    def plot_histogram(self,array,nbins = 15):
        fig, ax = plt.subplots()

        # the histogram of the data
        n, bins, patches = ax.hist(array, nbins, density=True)

        # add a 'best fit' line
        ax.set_xlabel('stress drop (Mpa)')
        ax.set_ylabel('Probability Density')
        ax.set_title(f'{self.name} Histogram of stress drop for {self.samples} samples' )
        fig_name =  f'{self.name}_stress_drop_histogram_nsamples_{self.samples}.jpg'
        fig_folder= self.make_dirs(file='figures_histogram')
        fig_file_dir = os.path.join(fig_folder,fig_name)   
        fig.savefig(fig_file_dir)


    def cov(self):
        dset  = self.read_h5file()
        self.covariance = np.cov(dset.transpose())
        
        nparameters = dset.shape[1]
        
        self.cov1 = self.covariance[:nparameters//3,:nparameters//3]
        self.cov2 = self.covariance[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
        self.cov3 = self.covariance[2*nparameters//3:,2*nparameters//3:]

        # cov12 = cov[:nparameters//3,nparameters//3:2*nparameters//3]
        # standard deviation (= square root of variance)

        self.std1 = np.sqrt(self.cov1.diagonal())
        self.std2 = np.sqrt(self.cov2.diagonal())
        self.std3 = np.sqrt(self.cov3.diagonal())

    def corr(self):
        dset  = self.read_h5file()
        self.correlation = np.corrcoef(dset.transpose())
        
        nparameters = dset.shape[1]
        self.corr1 = self.correlation[:nparameters//3,:nparameters//3]
        self.corr2 = self.correlation[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
        self.corr3 = self.correlation[2*nparameters//3:,2*nparameters//3:]

    
    def plot_cov(self):
        self.cov()
        self.Okada_input()
        

        cov = {'Normal':self.std1,
                'Along-strike':self.std2,
                'Along-dip':self.std3}

        x,y = self.xstn*1e-3,self.ystn*1e-3   
        x0,y0 = self.xsrc*1e-3,self.ysrc*1e-3
        ncol_target,nrow_target,patchid = self.find_slip() 
        nrows = len(y)
        ncols = len(x)

        n0rows = len(y0)
        n0cols = len(x0) 

        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=900)
        for i,parameter_id in enumerate(cov.keys()):
            parameter = cov[parameter_id]

            title = f"{self.name} Uncertainty in {parameter_id} Stress Change "

            im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow',linewidth=0.5)
            parray = np.zeros((n0rows,n0cols))
            masked_arr = np.ma.masked_array(parray,parray<1)
            #  commented below just for visualization purposes
            # axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
            axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1.5)
            axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=4.5)
            # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
            #X0, Y0 = np.meshgrid(x0, y0)
            #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
            try:
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
            except:
                pass 
            axes[i].margins(0)        
            fig.colorbar(im0, ax=axes[i],shrink=0.65,label='Uncertainty (MPa)')
            axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
            axes[i].set_aspect('equal', 'box')
            axes[i].set_title(title,fontweight='bold')
            axes[i].tick_params(labelsize=12)
        plt.tight_layout()
        # fig, axes = plt.subplots(1,3,figsize=(13,4.5),dpi=900)
    
        # for i,parameter_id in enumerate(cov.keys()):
        #     parameter = cov[parameter_id]

        #     title = f"{self.name} Uncertatinty \n  in {parameter_id} {self.observable} "

        #     im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),edgecolors='k', cmap='rainbow')
        #     parray = np.zeros((n0rows,n0cols))
        #     masked_arr = np.ma.masked_array(parray,parray<1)
        #     axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
        #     axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1.5)
        #     axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=4)
        #     # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
        #     #X0, Y0 = np.meshgrid(x0, y0)
        #     #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
        #     try:
        #         axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=10)
        #     except:
        #         pass 
        #     axes[i].margins(0)        
        #     cbar = fig.colorbar(im0, ax=axes[i],shrink=0.70,orientation='horizontal',label='Uncertainty (MPa)')
        #     axes[i].set_ylabel('Down-dip distance (km)',fontsize=8)
        #     axes[i].set_xlabel('Along-strike distance (km)',fontsize=8)
        #     axes[i].set_aspect('equal', 'box')
        #     axes[i].set_title(title,fontweight='bold',fontsize=9.5)
        #     axes[i].tick_params(labelsize=8)
        #     cbar.ax.tick_params(labelsize=8)
        # plt.tight_layout()
    
        fname =  f'std_{self.name}_{self.observable}_nsamples_{self.samples}.jpg'
        file_folder = self.make_dirs(file=f'figures_std')
        file_dir = os.path.join(file_folder,fname)   
        fig.savefig(file_dir)
        plt.close()

    def plot_corr(self):
        self.corr()
        self.Okada_input()
        # if patchID:
        #     ncol_target,nrow_target,patchid = self.find_slip(maxpatch = False,patchId=patchID)
        #     extra = str(patchid) + ' patch '
        # else:
        #     ncol_target,nrow_target,patchid = self.find_slip() 
        ncol_target,nrow_target,patchid = self.find_slip() 
        extra = 'max. slip patch'
        
        
        corr = {'normal':self.corr1[:,patchid],
                'along-strike':self.corr2[:,patchid],
                'along-dip':self.corr3[:,patchid]}

        x,y = self.xstn*1e-3,self.ystn*1e-3   
        x0,y0 = self.xsrc*1e-3,self.ysrc*1e-3

        nrows = len(y)
        ncols = len(x)
        n0rows = len(y0)
        n0cols = len(x0) 

        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        for i,parameter_id in enumerate(corr.keys()):
            parameter = corr[parameter_id]

            title = f"{self.name} Correlation in {parameter_id} {self.observable} "

            im0 = axes[i].pcolormesh(x,y,parameter.reshape(nrows,ncols),
                                 edgecolors='k', linewidth = 0.5,cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
            parray = np.zeros((n0rows,n0cols))
            masked_arr = np.ma.masked_array(parray,parray<1)
            # commeted for vidualixation purposes
            # axes[i].pcolormesh(x0,y0,masked_arr,edgecolors='white',linewidth=1)
            axes[i].plot(x[ncol_target]*np.ones(len(y)),y,ls='dashed',color='black',linewidth=1.5)
            axes[i].plot(x[ncol_target],y[nrow_target],marker='o',color='black',ms=4.5)
            # axes[i].scatter(x[ncol_target]*np.ones(len(y)),y,marker='_',color='black',s=3)
            #X0, Y0 = np.meshgrid(x0, y0)
            #axes[i].plot(X0.flat, Y0.flat, 'o', color='black',markersize=2)
            try:
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markeredgecolor='black',markersize=15)
            except:
                pass 
            axes[i].margins(0)        
            fig.colorbar(im0, ax=axes[i],shrink=0.65,label='Correlation')
            axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
            axes[i].set_aspect('equal', 'box')
            axes[i].set_title(title,fontweight='bold')
            axes[i].tick_params(labelsize=12)
        plt.tight_layout()
    
        fname =  f'corr_{extra}_{self.name}_{self.observable}_nsamples_{self.samples}.jpg'
        file_folder = self.make_dirs(file=f'figures_corr')
        file_dir = os.path.join(file_folder,fname)   
        fig.savefig(file_dir)
        plt.close()
        
    def save_cov(self):
        self.cov()
        fname =  f'covariance_{self.name}_{self.observable}_nsamples_{self.samples}.txt'
        file_folder = self.make_dirs(file='covariance')
        file_dir = os.path.join(file_folder,fname) 
        np.savetxt(file_dir,self.covariance)
        
    def save_corr(self):
        self.corr()
        fname =  f'correlation_{self.name}_{self.observable}_nsamples_{self.samples}.txt'
        file_folder = self.make_dirs(file='covariance')
        file_dir = os.path.join(file_folder,fname)
        np.savetxt(file_dir,self.correlation) 
    
    def plot_statistics(self):
        self.plot_corr()
        self.plot_cov()
    
    def save_statistics(self):
        self.save_corr()
        self.save_cov()
   