import sys
 
# adding csi functions to the system path

from csi.okadafull import displacement, stress, strain
from .observable import Observable
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm



# added step
class Stress(Observable):
    def __init__(self,
                    parent_dir,
                    name,
                    model_type,
                    step,
                    want_step,
                    patchsize,
                    shape_geom,
                    observable = 'Stress change',
                    samples='all',
                    new_station=False,
                    xstation = None,
                    ystation = None,
                    strikeOkada = 90,
                    factor = 4,
                    mu=30.0*1e9,
                    nu=0.25):
        
        
        self.mu = mu 
        self.nu = nu
        self.new_station = new_station
        
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
            
        
        self.save_data()
        self.save_plot()
        
    def modify_stn(self):
        if self.new_station:
                  
            self.zstn_flat = -self.depth.flatten()*1e3
            self.strike_stn_rad = self.strike_rad.flatten()
            self.dip_stn_rad = self.dip_rad.flatten()
        else:
            self.strike_stn_rad = self.strike_rad
            self.dip_stn_rad = self.strike_dip



    def Okada_stress(self):
        self.Okada_input()
        self.modify_stn()
        self.Stress = stress(
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
                                    self.ss, 
                                    self.ds, 
                                    self.ts, nu=self.nu)[0]
        
        self.Stress_df = pd.DataFrame(self.Stress,columns=['xx','xy','xz','yy','yz','zz'])



    def stress_change(self):
        self.Okada_stress()
        self.npatches = int(len(self.xstn)*len(self.ystn))
        tnn = np.zeros(self.npatches)
        ts1 = np.zeros(self.npatches)
        ts2 = np.zeros(self.npatches)

        Sxx = self.Stress_df['xx'].values #.reshape(ncols,nrows).transpose().flatten()
        Sxy = self.Stress_df['xy'].values #.reshape(ncols,nrows).transpose().flatten()
        Sxz = self.Stress_df['xz'].values #.reshape(ncols,nrows).transpose().flatten()
        Syy = self.Stress_df['yy'].values #.reshape(ncols,nrows).transpose().flatten()
        Syz = self.Stress_df['yz'].values #.reshape(ncols,nrows).transpose().flatten()
        Szz = self.Stress_df['zz'].values #.reshape(ncols,nrows).transpose().flatten()
        
        for k in range(self.npatches):
            
            sxx = Sxx[k]*1e-6
            sxy = Sxy[k]*1e-6
            sxz = Sxz[k]*1e-6
            syy = Syy[k]*1e-6
            syz = Syz[k]*1e-6
            szz = Szz[k]*1e-6
            
            # assemble stress tensor
            
            tau = [[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]]
            
            #tau = np.array([[S[k,0], S[k,1], S[k,2]],[S[k,1], S[k,3], S[k,4]], [S[k,2], S[k,4], S[k,5]]])
            
            # set cosine/sine constants
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

        stressChange_d = {'Normal':tnn,'Along-strike Shear':ts1,'Along-dip Shear':ts2}
        self.stressChange_df  = pd.DataFrame(stressChange_d)
    


    def stress_drop(self,comp='ts2',weighted = True):
        self.stress_change()
        # need minus sign for stress drop
        stress_drop  = - self.stressChange_df[comp].values
        
        if weighted:
            average_stress = np.average(stress_drop,weights = dip_slip)
        
        else:
            average_stress = np.mean(stress_drop)
        
        self.average_stress = average_stress


    def stress_drop_Noda(self):
        self.stress_change()
        # need minus sign for stress drop
        as_shear = -self.stressChange_df['Along-strike Shear'].values
        ad_shear = -self.stressChange_df['Along-dip Shear'].values 

        delta_stress = np.array([as_shear,ad_shear]).T
        delta_u = np.array([self.ss,self.ds]).T
        prod_stress_u = np.sum(delta_stress*delta_u,axis=1)
        average_stress = np.sum(prod_stress_u)/np.linalg.norm(np.sum(delta_u,axis=1))        
        self.average_stress = average_stress

    def make_dirs(self,file='data'):
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir=file)
        os.makedirs(folder_dir,exist_ok=True)
        return folder_dir
    
    def save_data(self):
        self.stress_change()
        fname = f'{self.name}_{self.samples}_{self.model_type}_model_{self.observable}.txt'
        file_dir = os.path.join(self.make_dirs(),fname)
        np.savetxt(file_dir,self.stressChange_df)
        
    def plot_dir(self,extra=' ',suffix='png'):
        fname = f'{self.name}_{self.samples}_{self.model_type}_model_{self.observable}_{extra}.{suffix}'
        file_dir = os.path.join(self.make_dirs(file='figures'),fname)
        return file_dir
    
    def save_plot(self,suffix = 'png'):
        self.stress_change()
        xstn,ystn = self.xstn*1e-3,self.ystn*1e-3 # convert to km again
        Xstn, Ystn = np.meshgrid(xstn,ystn)  
        shape_stn = (len(ystn),len(xstn))
        fig,ax = plt.subplots(dpi=600)
        supertitle = f"{self.name} {self.observable} from {self.samples} samples"
        
        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=900)
        
        for i,parameter_id in enumerate(self.stressChange_df.keys()):
            
            parameter = self.stressChange_df[parameter_id].values
            title = f"{self.name} {parameter_id} {self.observable}"

            # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
            im = axes[i].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
            
            #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
            axes[i].margins(0)
                    
            fig.colorbar(im, ax=axes[i],shrink=0.96,label='{}'.format('Stress (MPa)'))
            try:
                
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markersize=18,markeredgecolor='black')
            except:
                pass

             
            axes[i].set_ylabel('Down-dip distance (km)',fontsize=11) # before 12
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=11)
            axes[i].set_aspect('equal', 'box') 
            axes[i].set_title(title,fontweight='bold',fontsize = 13) # before not specfied
            axes[i].tick_params(labelsize=12)

        # fig.suptitle(supertitle,x=0.5,y=1,fontweight='bold') # uncomment after poster
        plt.tight_layout()
        figure_dir = self.plot_dir()
        fig.savefig(figure_dir)
       