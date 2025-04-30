

import sys
 
# adding csi functions to the system path

from csi.okadafull import displacement, stress, strain
from .FastSweep import Hypocenter, Grid, FastSweep 
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


# added self.step 
class Observable():
    def __init__(self,
                 parent_dir,
                 name,
                 model_type,
                 step,
                 want_step,
                 patchsize,
                 shape_geom,
                 observable,
                 samples,
                 new_station,
                 xstation ,
                 ystation,
                 strikeOkada,
                 factor):
        
        
        self.parent_dir = parent_dir
        self.input_dir = os.path.join(parent_dir,'INPUT')
        self.name = name
        self.model_type = model_type
        #added step 
        self.step = step
        self.want_step = want_step
        self.samples = samples
        self.observable = observable
        self.new_station = new_station
        self.xstation = xstation
        self.ystation = ystation
        self.factor = factor
        

        # added self.step
        if want_step:
            self.list_out = [self.parent_dir,self.name,'model',self.model_type,'step_'+str(self.step),self.observable,str(self.samples) + '_samples','data']
            # added self.step 
            self.list_in = [self.parent_dir,self.name,'model',self.model_type,'step_'+str(self.step), str(self.samples) + '_samples','mean']
        else:
            self.list_out = [self.parent_dir,self.name,'model',self.model_type,self.observable,str(self.samples) + '_samples','data']
            # added self.step 
            self.list_in = [self.parent_dir,self.name,'model',self.model_type, str(self.samples) + '_samples','mean']
      

        self.patchsize = patchsize
        self.patch = self.patchsize*1e3  # in meters
        self.shape_geom = shape_geom
        self.nrows = shape_geom[0] 
        self.ncols = shape_geom[1] 
        self.strikeOkada = strikeOkada

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
        self.df = df.drop(df.columns[0],axis=1)
        
        self.Uparallel = self.array_formatter(self.df['U_parallel'].values)
        self.Uperp = self.array_formatter(self.df['U_perp'].values)
        self.Slip = self.array_formatter(self.df['Slip'].values)
        self.dip = self.array_formatter(self.df['dip'].values)
        self.depth = self.array_formatter(self.df['depth'].values)
        self.strike = self.array_formatter(self.df['strike'].values)
        
        
        
        
        try:
            self.Vr = self.array_formatter(self.df['Vr'].values).reshape(self.shape_geom) # need 2D array for Vr .flatten when looking for 1D array
            self.Tr = self.array_formatter(self.df['Tr'].values)
            self.hypocenter = [df['Hypo_as'][0],-df['Hypo_dd'][0]] # negative down-dip component of hypocenter
        except:
            print('Dealing with static model')
        
        return 
    
    def proj_ysrc_coords(self):
        self.mean_model_reader()
        dip = self.df['dip'].values[:self.nrows]
        proj_dysrc = -self.patch*np.cos(dip*np.pi/180) # in meters
        proj_ysrc = np.zeros_like(proj_dysrc)
        for i in range(len(proj_ysrc)):
            proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
        ysrc = np.flip(proj_ysrc)
        
        return ysrc
    
    def source(self):
        self.mean_model_reader()
     
        # DISLOCATION
        self.xsrc = np.arange((1/2)*self.patch,self.ncols*self.patch, self.patch)
        self.ysrc = np.arange(-(self.nrows-1/2)*self.patch,0, self.patch)
        # when ysrc must be correct; projected does not work with FastSweep
        #self.ysrc = self.proj_ysrc_coords()
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
            self.source()
            self.xstn = self.xstation
            self.ystn = self.ystation
            # below world for finer grid 
            # self.zstn_flat = np.zeros(len(self.xstn)*len(self.ystn))
            #  this was changed ; would not work for stress calculations
            # self.zstn_flat = -self.depth*1e3
            
        else:
            # works for displacement
            offsetx =  (self.ncols//self.factor)*self.patch
            offsety =  (self.nrows//self.factor)*self.patch
            self.xstn = np.arange(self.patch/2 - offsetx, self.ncols*self.patch + offsetx,self.patch)
            self.ystn = np.arange(-offsety - (self.nrows-1/2)*self.patch, offsety,self.patch)
            self.zstn_flat = np.zeros(len(self.xstn)*len(self.ystn))
        
        Xstn ,Ystn = np.meshgrid(self.xstn,self.ystn)
        
        self.xstn_flat = Xstn.flatten()
        self.ystn_flat = Ystn.flatten() 

            
    def find_slip(self):
        self.source()
        self.station()   

        n0cols,n0rows = len(self.xsrc),len(self.ysrc)
        
        max_slip_patch = np.argmin(abs(self.Slip-max(self.Slip)))

        
        col0_max_slip,row0_max_slip = max_slip_patch%n0cols, max_slip_patch//n0cols
        
        
        # displacement is defined on customized geometry overlapping prescribed surface 

        
        ncols,nrows = len(self.xstn),len(self.ystn)
        extra_cols = (ncols-n0cols)//2
        extra_rows = (nrows-n0rows)//2
        
        self.ncol_target = col0_max_slip + extra_cols 
        self.nrow_target = row0_max_slip + extra_rows
        self.max_slip_patch_in_stn = self.nrow_target*ncols + self.ncol_target 
        

        return self.ncol_target,self.nrow_target, self.max_slip_patch_in_stn
            
    def Okada_input(self):
        self.source()
        self.station()
        
        self.ts = np.zeros(self.xsrc_flat.shape)
        self.ss = self.Uperp
        self.ds = self.Uparallel

    def get_T0(self):
        self.source()
        X_src,Y_src = np.meshgrid(self.xsrc*1e-3,self.ysrc*1e-3)
        x_hyp = float(self.hypocenter[0])
        y_hyp = float(self.hypocenter[1])

        # hypoctr = Hypocenter(x=x_hyp,y=y_hyp)
        # grid0 = Grid(x=X_src,y =Y_src, vr= self.Vr) 
        solver = FastSweep()
        solver.setGrid(self.xsrc*1e-3,self.ysrc*1e-3,self.Vr)
        solver.setHypo(x_hyp,y_hyp)
        solver.fastSweep(verbose=False)
        self.T0 = solver.t0[1:-1,1:-1].flatten()

    def save_T0(self):
        self.get_T0()
        t0 = self.T0
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir='Rupture_data')
        os.makedirs(folder_dir,exist_ok=True)
        fname = f'{self.observable}_{self.name}_{self.samples}_{self.model_type}_model.txt'
        file_dir = os.path.join(folder_dir,fname)
        np.savetxt(file_dir,t0)
                
    def triangularSTF(self,t,slip,tr,t0):
        h =2*slip/tr
        P1 = (t0,0)
        P2 = (t0 + tr/2,h)
        P3 = (t0 + tr,0)
        m1  =  (P2[1] - P1[1])/(P2[0] - P1[0]) 
        m2  = (P3[1] - P2[1])/(P3[0] - P2[0])

        s = np.zeros_like(t)
        for i in range(len(t)):
            if t0<t[i]<t0 + tr/2:
                s[i] = P1[1] + (t[i] - P1[0])*m1
            elif t0 + tr/2<t[i]<t0 + tr:
                s[i] = P2[1] + (t[i] - P2[0])*m2
            else:
                s[i] = 0
        
        return s


    def cum_slip_triangularSTF(self,t,slip,tr,t0):
        """
        closed-form for triangular source time function 
        """

        h =2*slip/tr
        P1 = (t0,0)
        P2 = (t0 + tr/2,h)
        P3 = (t0 + tr,0)
        m1  =  (P2[1] - P1[1])/(P2[0] - P1[0]) 
        m2  = (P3[1] - P2[1])/(P3[0] - P2[0])

        s = np.zeros_like(t)
        for i in range(len(t)):
            if t[i]<=t0:
                s[i] = 0
            elif t0<t[i]<t0 + tr/2:
                s[i] = (P1[1] - P1[0]*m1)*(t[i]-t0) + (m1/2)*(t[i]**2 -t0**2) 
            elif t0 + tr/2<=t[i]<t0 + tr:
                s[i] = slip/2 + (P2[1] - P2[0]*m2)*(t[i] - (t0 + tr/2)) + (m2/2)*(t[i]**2 -(t0 + tr/2)**2) 
            elif t[i]>= t0 + tr:  
                s[i] = slip
        return s

    def propagate_rupture(self,time):
        self.mean_model_reader()
        self.save_T0()

        PropagateT = np.zeros((len(self.T0),len(time)))
        for i,(t0,tr,slip) in enumerate(zip(self.T0,self.Tr,self.Slip)):
            PropagateT[i,:] = self.cum_slip_triangularSTF(time,slip,tr,t0)-self.cum_slip_triangularSTF(time-10,slip,tr,t0)

        return PropagateT
    
    def plot_rupture(self,dt=10,yloc=0.88,smax=30,hspace=0.25):
        time  = np.arange(0,dt*15+dt,dt)
        rupture = self.propagate_rupture(time) # first dt at index 0 not meaningful 
        

        fig,axes =plt.subplots(5,3,figsize =(8,12),dpi = 500,sharey=True,sharex=True)
        plt.subplots_adjust(wspace=0.10, hspace=hspace)
        for i in range(15):
            im = axes[i//3][i%3].pcolormesh(self.xsrc*1e-3,self.ysrc*1e-3,rupture[:,i+1].reshape(self.nrows,self.ncols),
            edgecolors='none',cmap='hot_r',vmin =0, vmax =smax)
            textstr = f'{i*dt} s - {(i+1)*dt} s'

            # # these are matplotlib.patch.Patch properties
            # props = dict(boxstyle='round', facecolor='none')

            # # place a text box in upper left in axes coords
            # axes[i//3][i%3].text(0.85, 0.05, textstr, transform=axes[i//3][i%3].transAxes, fontsize=7,
            #     verticalalignment='top', bbox=props)
            axes[i//3][i%3].set_title(textstr,fontsize=8)
            axes[i//3][i%3].set_aspect('equal', 'box')
             
        fig.colorbar(im,ax= axes,shrink=0.4,pad=0.05,label= 'Slip (m)',orientation='horizontal')
        fig.suptitle(f'{self.name} Rupture evolution',y=yloc)
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir='Rupture_figure')
        os.makedirs(folder_dir,exist_ok=True)
        if self.want_step:

            fname = f'{self.observable}_{self.name}_{self.samples}_{self.model_type}_model_step_{self.step}.png'
        else:
            fname = f'{self.observable}_{self.name}_{self.samples}_{self.model_type}.png'

        file_dir = os.path.join(folder_dir,fname)
        fig.savefig(file_dir)
                
            
