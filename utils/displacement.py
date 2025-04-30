

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

class Displacement(Observable):
    def __init__(self,
                    parent_dir,
                    name,
                    model_type,
                    step,
                    want_step,
                    patchsize,
                    shape_geom,
                    observable = 'displacement',
                    samples='all',
                    new_station=False,
                    xstation = None,
                    ystation = None,
                    strikeOkada = 90,
                    GF_on = True,
                    scale = 0.25,
                    larrow = 10,
                    factor = 4,
                    mu=30.0*1e9,
                    nu=0.25):
        
        self.scale =scale
        self.new_station = new_station
        self.larrow  =larrow
        self.mu = mu 
        self.nu = nu
        self.GF_on = GF_on
        
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
            
        self.Okada_execution()
        
    def modify_stn(self):
        if self.new_station:
            self.zstn_flat = np.zeros(len(self.xstn)*len(self.ystn))          

    def Okada_displacement(self):
        self.Okada_input()
        self.modify_stn()
        self.Displacement = displacement(
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
                                    self.ts, nu=self.nu)
        
        self.Displacement_df = pd.DataFrame(self.Displacement,columns=['x','y','z'])
        
    def make_dirs(self,file='data'):
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir=file)
        os.makedirs(folder_dir,exist_ok=True)
        return folder_dir
    
    def save_data(self):
        self.Okada_displacement()
        fname = f'{self.name}_{self.samples}_{self.model_type}_model_{self.observable}.txt'
        file_dir = os.path.join(self.make_dirs(),fname)
        np.savetxt(file_dir,self.Displacement)
        
    def save_GF(self):
        self.GF()
        fname = f'{self.name}_{self.samples}_{self.model_type}_model_{self.observable}_GF.txt'
        file_dir = os.path.join(self.make_dirs(file='GF'),fname)
        np.savetxt(file_dir,self.G)

        
    def plot_dir(self,extra=' ',suffix='png'):
        fname = f'{self.name}_{self.samples}_{self.model_type}_model_{self.observable}_{extra}.{suffix}'
        file_dir = os.path.join(self.make_dirs(file='figures'),fname)
        return file_dir
    

    def GF(self):
        
        Nobservables = len(self.xstn)*len(self.ystn)
        Npatches = self.nrows*self.ncols
        ts = np.zeros(self.xsrc_flat.shape)
        self.G = np.zeros((3*Nobservables,2*Npatches))
        for slip in ['ss','ds']:
            for j in range(Npatches):

                m = {'ss':np.zeros(Npatches),'ds':np.zeros(Npatches)}

                if slip == 'ds':
                    Np = Npatches
                else:
                    Np = 0

                m[slip][j] = 1
                Displacement = displacement(self.xstn_flat,
                                    self.ystn_flat,
                                    self.zstn_flat, 
                                    self.xsrc_flat, 
                                    self.ysrc_flat, 
                                    self.zsrc_flat, 
                                    self.width, 
                                    self.length, 
                                    self.strike_rad, 
                                    self.dip_rad, 
                                    m['ss'],
                                    m['ds'],  
                                    ts, nu=self.nu)
                dj = Displacement.flatten(order='F')
                self.G[:,j + Np] = dj

        return None

    def validate_GF(self):
        self.GF()

        mean_m = np.array([self.Uperp,self.Uparallel]).flatten()
        displacement_from_GF = np.matmul(self.G,mean_m).reshape(self.G.shape[0]//3,3,order='F')
        displacement_from_GF_df = pd.DataFrame(displacement_from_GF,columns=['x','y','z']) 
        self.save_plot(displacement_from_GF_df,suptitle_extra ='from GF',fname_extra= 'from GF')
    
    def Okada_execution(self):
        self.Okada_displacement()
        self.save_plot(self.Displacement_df)
        self.save_data()
        if self.GF_on:
            self.validate_GF()
            print("Validating ---------")
            self.save_GF()

            

    def save_plot(self,df,suptitle_extra ='',fname_extra= '',suffix = 'png'):
        
        xstn,ystn = self.xstn*1e-3,self.ystn*1e-3 # convert to km again
        if self.name == 'Gorkha':
            ystn = ystn - 76
            self.hypocenter[1] = self.hypocenter[1] - 76
        Xstn, Ystn = np.meshgrid(xstn,ystn)  
        shape_stn = (len(ystn),len(xstn))
        fig,ax = plt.subplots(dpi=900)
        # supertitle = f"{self.name} Surface Displacement from {self.samples} samples {suptitle_extra}"
        supertitle = f"{self.name} Surface Displacement"
        Ux = df['x'].values
        Uy = df['y'].values
        Uz = df['z'].values
    
        im = ax.pcolormesh(xstn,ystn,Uz.reshape(shape_stn),edgecolors='k',linewidths=0.1,cmap='bwr',norm = TwoSlopeNorm(0,vmin = min(Uz.min(),-Uz.max()),vmax=max(-Uz.min(),Uz.max()))) 
        fig.colorbar(im,ax = ax,shrink = 0.4,label= 'Uplift (m)')
        try:
            ax.plot(self.hypocenter[0],self.hypocenter[1],marker = '*',color='yellow',markersize=14,markeredgecolor='black')
        except:
            pass
        q = ax.quiver(Xstn,Ystn,Ux,Uy,scale=self.scale,scale_units ='x', units='width',width=0.002,headwidth=4.5,headlength=6)
        ax.quiverkey(q, X=0.04, Y=1.04, U=self.larrow,label=f'{self.larrow}m', labelpos='N',fontproperties={'size':7})
        ax.set_ylabel('Trench-normal distance (km)',fontsize=9)
        ax.set_xlabel('Along-strike distance (km)',fontsize=9)
        ax.set_aspect('equal', 'box')
        ax.tick_params(labelsize=9)
        ax.set_title(supertitle,x=0.5,y=1,fontsize=13, fontweight='bold') # size before 9
        plt.tight_layout()
        
        figure_dir = self.plot_dir(extra=f'allin1_{fname_extra}')
        fig.savefig(figure_dir)    
        plt.close()
        
        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        
        for i,parameter_id in enumerate(df.keys()):
            
            parameter = df[parameter_id].values
            title = f"{parameter_id} Surface displacement"

            # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
            im = axes[i].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
            
            #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
            axes[i].margins(0)
                    
            fig.colorbar(im, ax=axes[i],shrink=0.8,label='{}'.format('Displacement (m)'))
            try:
                
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markersize=18,markeredgecolor='black')
            except:
                pass

             
            axes[i].set_ylabel('Trench-normal distance (km)',fontsize=11)
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=11)
            axes[i].set_aspect('equal', 'box')
            axes[i].set_title(title,fontweight='bold',fontsize = 13)
            axes[i].tick_params(labelsize=12)
            #axes[i].invert_yaxis()
        # fig.suptitle(supertitle,x=0.5,y=1,fontweight='bold') # uncomment later
        plt.tight_layout()
        figure_dir = self.plot_dir(extra=f'{fname_extra}')
        fig.savefig(figure_dir)
        plt.close()
        
