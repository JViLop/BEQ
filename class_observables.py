# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 10:25:13 2024

@author: joanv
"""
import sys
 
# adding csi functions to the system path
sys.path.insert(0, '/home/josevilo/csi/csi')

from os import displacement, stress, strain
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm



class Observable():
    def __init__(self,
                 parent_dir,
                 name,
                 model_type,
                 patchsize,
                 shape_geom,
                 samples='all',
                 new_station='False',
                 xstation = None,
                 ystation = None,
                 strikeOkada = 90):
        
        
        self.parent_dir = parent_dir
        self.input_dir = os.path.join(parent_dir,'INPUT')
        self.name = name
        self.model_type = model_type
        self.samples = samples
        self.new_station = new_station
        self.xstation = xstation
        self.ystation = ystation
        
        self.list_out = [self.parent_dir,self.name,'model',self.model_type,self.observable,str(self.samples) + '_samples','data']
        self.list_in = [self.parent_dir,self.name,'model',self.model_type,str(self.samples) + '_samples','mean']
        
        self.patchsize = patchsize
        self.patchmeter = self.patchsize*1e3  # in meters
        self.shape_geom = shape_geom
        self.nrows = shape_geom[0] 
        self.ncols = shape_geom[1] 
        self.strikeOkada = strikeOkada
        self.save_data()
        self.save_plot()

    # def silent(func):
    #     def inner(self,*args,**kwargs):
    #         return func(self,*args,**kwargs)
    #     return inner
    
    
    def dir_manager(self,level,ldir,data_type,tail_dir='data'):
        ldir = ldir
        ldir[-1]=tail_dir
        
        if level ==0:
            path = os.path.join(self.parent_dir,data_type) 
            return path
        else:
            path = os.path.join(self.dir_manager(level-1,ldir,data_type,tail_dir=tail_dir),ldir[level])
            return path
        
    def dir_creator(self,level,tail_dir='data'):
        return self.dir_manager(level,self.list_out,'OUTPUT',tail_dir=tail_dir)
        
    def dir_finder(self,level,tail_dir='mean'):
        return self.dir_manager(level,self.list_in,'INPUT',tail_dir=tail_dir)

    def array_formatter(self,array):
        return np.flip(array.reshape(self.shape_geom,order='F'),axis=0).flatten()
        
    def mean_model_reader(self):
        file_folder = self.dir_finder(len(self.list_in)-1,tail_dir='mean')
        file_dir = os.path.join(file_folder,f'{self.name}_mean_{self.model_type}_model.csv')
        df = pd.read_csv(file_dir)
        df = df.drop(df.columns[0],axis=1)
        
        self.Uparallel = self.array_formatter(df['U_parallel'].values)
        self.Uperp = self.array_formatter(df['U_perp'].values)
        self.Slip = self.array_formatter(df['Slip'].values)
        self.dip = self.array_formatter(df['dip'].values)
        self.depth = self.array_formatter(df['depth'].values)
        self.strike = self.array_formatter(df['strike'].values)
        
        try:
            self.hypocenter = [df['Hypo_as'][0],-df['Hypo_dd'][0]] # negative down-dip component of hypocenter
        except:
            print('No Hypocenter available')
        
        return 
    
    def source(self):
        self.mean_model_reader()
     
        # DISLOCATION
        self.xsrc = np.arange((1/2)*self.patch,self.ncols* self.patch, self.patch)
        self.ysrc = np.arange(-(self.nrows-1/2)*self.patch,0, self.patch)
        Xsrc,Ysrc= np.meshgrid(self.xsrc,self.ysrc)
        
        self.xsrc_flat = Xsrc.flatten()
        self.ysrc_flat = Ysrc.flatten()
        self.zsrc_flat = self.depth*1e3
        
        self.length =  self.patch*np.ones(self.xsrc_flat.shape)
        self.width =  self.patch*np.ones(self.ysrc_flat.shape)
        self.dip_rad = self.dip*(np.pi)/180
        self.strike_rad = self.strikeOkada*np.ones(self.xsrc_flat.shape)*(np.pi/180)
    def station(self):
        if self.new_station:
            self.xstn = self.xstation
            self.ystn = self.ystation
        else:
            offsetx =  (self.ncols//self.factor)*self.patchmeter
            offsety =  (self.nrows//self.factor)*self.patchmeter
            self.xstn = np.arange(self.patch/2 - offsetx, self.ncols*self.patch + offsetx,self.patch)
            self.ystn = np.arange(-offsety - (self.nrows-1/2)*self.patch, offsety,self.patch)
        
        Xstn ,Ystn = np.meshgrid(self.xstn,self.ystn)
        
        self.xstn_flat = Xstn.flatten()
        self.ystn_flat = Ystn.flatten() 
        self.zstn_flat = np.zeros(self.xstn_flat.shape)   # for station, 0-depth array coincident with  seafloor
        
        
    def Okada_input(self):
        self.source()
        self.station()
        
        self.ts = np.zeros(self.xsrc_flat.shape)
        self.ss = self.Uperp
        self.ds = self.Uparallel
        
      



class Displacement(Observable):
    def __init__(self,parent_dir,
                     name,
                     model_type,
                     patchsize,
                     shape_geom,
                     new_station='False',
                     xstation = None,
                     ystation = None,
                     strikeOkada = 90,
                     scale = 0.25,
                     larrow = 10,
                     factor = 4,
                     observable='displacement',
                     samples='all',
                     mu=30.0*1e9,
                     nu=0.25):
        
        self.scale =scale
        self.larrow  =larrow
        self.factor = factor
        self.mu = mu 
        self.nu = nu
    
        
        super().__init__(self,parent_dir,
                         name,
                         model_type,
                         patchsize,
                         shape_geom,
                         samples='all',
                         new_station='False',
                         xstation = None,
                         ystation = None,
                         strikeOkada = 90)
            
        
        
                                        
    def Okada_displacement(self):
        self.Okada_input()
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
        
    def plot_dir(self,extra=' ',suffix='png'):
        fname = f'{self.name}_{self.samples}_{self.model_type}_model_{self.observable}_{extra}.{suffix}'
        file_dir = os.path.join(self.make_dirs(file='figures'),fname)
        return file_dir
    
    def save_plot(self,suffix = 'png'):
        self.Okada_displacement()
        xstn,ystn = self.xstn*1e-3,self.ystn*1e-3 # convert to km again
        Xstn, Ystn = np.meshgrid(xstn,ystn)  
        shape_stn = (len(ystn),len(xstn))
        fig,ax = plt.subplots(dpi=600)
        supertitle = f"{self.name} Surface Displacement from {self.samples} samples"
        Ux = self.Displacement_df['x'].values
        Uy = self.Displacement_df['y'].values
        Uz = self.Displacement_df['z'].values
    
        im = ax.pcolormesh(xstn,ystn,Uz.reshape(shape_stn),edgecolors='k',linewidths=0.25,cmap='bwr',norm = TwoSlopeNorm(0,vmin = min(Uz.min(),-Uz.max()),vmax=max(-Uz.min(),Uz.max()))) 
        fig.colorbar(im,ax = ax,shrink = 0.35,label= 'Uplift (m)')
        try:
            ax.plot(self.hypocenter[0],self.hypocenter[1],marker = '*',color='yellow',markersize=14,markeredgecolor='black')
        except:
            pass
        q = ax.quiver(Xstn,Ystn,Ux,Uy,scale=self.scale,scale_units ='x', units='width',width=0.002,headwidth=4.5,headlength=6)
        ax.quiverkey(q, X=0.04, Y=1.04, U=self.larrow,label=f'{self.larrow}m', labelpos='N',fontproperties={'size':7})
        ax.set_ylabel('Down-dip distance (km)',fontsize=9)
        ax.set_xlabel('Along-strike distance (km)',fontsize=9)
        ax.set_aspect('equal', 'box')
        ax.tick_params(labelsize=9)
        ax.set_title(supertitle,x=0.5,y=1,fontsize =9, fontweight='bold')
        plt.tight_layout()
        
        figure_dir = self.plot_dir(extra='allin1')
        fig.savefig(figure_dir)    
        plt.close()
        
        fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
        
        for i,parameter_id in enumerate(self.Displacement_df.keys()):
            
            parameter = self.Displacement_df[parameter_id].values
            title = f"{parameter_id} Surface displacement"

            # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
            im = axes[i].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
            
            #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
            axes[i].margins(0)
                    
            fig.colorbar(im, ax=axes[i],shrink=0.9,label='{}'.format('Displacement (m)'))
            try:
                
                axes[i].plot(self.hypocenter[0],self.hypocenter[1],marker='*',color='yellow',markersize=18,markeredgecolor='black')
            except:
                pass

             
            axes[i].set_ylabel('Down-dip distance (km)',fontsize=12)
            axes[i].set_xlabel('Along-strike distance (km)',fontsize=12)
            axes[i].set_aspect('equal', 'box')
            axes[i].set_title(title,fontweight='bold')
            axes[i].tick_params(labelsize=12)
            #axes[i].invert_yaxis()
        fig.suptitle(supertitle,x=0.5,y=1,fontweight='bold')
        plt.tight_layout()
        figure_dir = self.plot_dir()
        fig.savefig(figure_dir)
        plt.close()
        
