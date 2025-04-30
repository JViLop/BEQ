# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 18:16:16 2024

@author: joanv
"""

from WaveInt_jobs import WaveInt
# from wavemod.sacpy import sac
import numpy as np
import os
import pandas as pd
import obspy as obs
import re
import multiprocessing as mp
import sys


class dynamicGF():
    def __init__(self,
                patchsize,
                patchID,
                dip,
                depth,
                Tr,
                Xsrc,
                model_file,
                xstn,
                ystn,
                output_dir):
    
        self.patchsize = patchsize
        self.patchID = patchID
        self.dip = dip
        self.depth = depth
        self.Tr = Tr 
        self.Xsrc = Xsrc
        self.model_file = model_file
        self.xstn = xstn
        self.ystn = ystn
        self.output_dir = output_dir
        self.build_model_dict()
        # self.set_stations()
    
    def set_stations(self):
        Xstn, Ystn = np.meshgrid(self.xstn,self.ystn)
        xstn_flat,ystn_flat = Xstn.flatten(), Ystn.flatten()
        IDs = np.arange(0,len(xstn_flat),1)
        d = np.array([IDs,xstn_flat,ystn_flat]).T
        stn_file_name = 'stations.txt'
        path = os.path.join(self.output_dir,stn_file_name)
        np.savetxt(path,d,fmt = ['%.4d','%.2f','%.2f'],header = "STNAME X_COORD Y_COORD")
        self.stn_file_path = path


    def get_mu(self):
        heights = self.model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(self.depth - depths))
        if self.depth > depths[i]:
             i +=1 
        rho = self.model_parameters['RHO'][i]*1e3 # convert to SI
        vs = self.model_parameters['VS'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        self.mu = mu
    
    def get_scalar_Moment(self,slip):
        self.get_mu()
        A = (self.patchsize*1e3)**2
        M0 = (A*slip*self.mu)*1e7  # dyne.cm 
        self.M0 = M0 

    def build_model_dict(self):
        path = os.path.join(self.output_dir,self.model_file)
        f = open(path, "r")
        content = f.readlines()
        model = content[12:]
        keys = re.findall('\w+',content[11],flags=re.DOTALL)
        fmt = re.compile('\d+.\d+')
        layers = [list(map(float,re.findall(fmt,m))) for m in model]
        layers = np.array(layers).T
        self.model_parameters = dict(zip(keys,layers))

    def calc_synthetics(self,npts,dt,rake,T0,observable,dir,stn_file_name):
        self.get_scalar_Moment(slip=1)
        # path_model = os.path.join(self.output_dir,self.model_file) 
        x = WaveInt(self.model_file,npts,dt,T0=T0,Xs=self.Xsrc)
        path_stn = os.path.join(self.output_dir,stn_file_name)
        x.readStatXY(path_stn)
        x.checkXs()
        x.calcDist(station_file=path_stn)
        x.writeDistFile(self.output_dir)
        # x.writeModelFile([5,5.7,6.0],[2.7,3.3,3.4],[2500,2700,2750],[2.8,5,8.2])
        x.calcKernel(self.output_dir)
        x.synthKernelSDR(observable,0,self.dip,rake,self.M0,'triangle',dir = dir, src_id = self.patchID,duration=self.Tr,output_dir=self.output_dir,rfile=None)

    def calc_GF(self,npts,dt,observable,sac_folder,stn_file_name,T0=-5):
        
        ss_sac_folder  = os.path.join(sac_folder,'along_strike') 
        dd_sac_folder  = os.path.join(sac_folder,'along_dip')
        os.makedirs(ss_sac_folder,exist_ok=True)
        os.makedirs(dd_sac_folder,exist_ok=True)

        self.calc_synthetics(npts,dt,0,T0,observable,ss_sac_folder,stn_file_name)
        self.calc_synthetics(npts,dt,90,T0,observable,dd_sac_folder,stn_file_name)







            



class DynamicGF():
    def __init__(self,
                parent_dir:str, 
                name:str,
                patchsize:float,
                geom:tuple,
                samples,
                variable:str,
                layered_model_file,
                xstns:np.array,
                ystns:np.array,
                job_id):
        
        self.parent_dir = parent_dir
        self.input_dir = os.path.join(parent_dir,'INPUT')
        
        self.name = name
        self.model_type = 'kinematic'
        self.patchsize = patchsize    # km
        self.shape_geom = geom
        self.nrows = geom[0] 
        self.ncols = geom[1] 
        
        self.variable = variable
        self.samples = samples
        self.list_in = [self.parent_dir,self.name,'model',self.model_type,str(self.samples) + '_samples','mean']
        self.list_out = [self.parent_dir,self.name,'model',self.model_type,self.variable,str(self.samples) + '_samples','data']

        self.layered_model_file = layered_model_file 
        self.npts = 10 
        self.dt = 0.5
        self.T0 = 0  # not to be confused with t0 (initial rise time array)
        self.xstns = xstns
        self.ystns = ystns
        self.job_id = job_id

        self.set_t0() # must be changed this is just for test purposes
        self.set_stations()
        self.get_responses_all_patches()
    
    def set_stations(self):
        self.folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir='job_'+self.job_id)
        os.makedirs(self.folder_dir,exist_ok=True)
        Xstn, Ystn = np.meshgrid(self.xstns,self.ystns)
        xstn_flat,ystn_flat = Xstn.flatten(), Ystn.flatten()
        IDs = np.arange(0,len(xstn_flat),1)
        d = np.array([IDs,xstn_flat,ystn_flat]).T
        stn_file_name = f'stations_{self.job_id}.txt'
        path_stn = os.path.join(self.folder_dir,stn_file_name)
        np.savetxt(path_stn,d,fmt = ['%.4d','%.2f','%.2f'],header = "STNAME X_COORD Y_COORD")
        self.stn_file_name = path_stn
    

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
    def set_t0(self):
        file_folder = self.dir_finder(len(self.list_in)-1,tail_dir='mean')
        file_dir = os.path.join(file_folder,f'RuptureTime_{self.name}_{self.samples}_kinematic_model.txt')
        self.t0 = np.loadtxt(file_dir)
    
    def proj_ysrc_coords(self):
        self.mean_model_reader()
        dip = self.df['dip'].values[:self.nrows]
        proj_dysrc = -self.patchsize*np.cos(dip*np.pi/180) # in meters
        proj_ysrc = np.zeros_like(proj_dysrc)
        for i in range(len(proj_ysrc)):
            proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
        ysrc = np.flip(proj_ysrc)
        
        return ysrc
    
    def source(self):
        self.mean_model_reader()
     
        # DISLOCATION
        self.xsrc = np.arange((1/2)*self.patchsize,self.ncols*self.patchsize, self.patchsize)
        # self.ysrc = np.arange(-(self.nrows-1/2)*self.patch,0, self.patch)
        self.ysrc = self.proj_ysrc_coords()
        Xsrc,Ysrc= np.meshgrid(self.xsrc,self.ysrc)
        
        self.xsrc_flat = Xsrc.flatten()
        self.ysrc_flat = Ysrc.flatten()
        self.zsrc_flat = self.depth  # it mush be in km (not meters)
        
        self.Xsrcs  = np.array([self.xsrc_flat,self.ysrc_flat,self.zsrc_flat]).T

    def get_response_single_patch(self,patchID,dip,depth,Tr,xsrc):
        
        self.folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir='job_'+self.job_id)
        # path_model = os.path.join(self.folder_dir,self.layered_model_file)
        os.makedirs(self.folder_dir,exist_ok=True)
        response = dynamicGF(self.patchsize,
                patchID,
                dip,
                depth,
                Tr,
                xsrc,
                self.layered_model_file,
                self.xstns,
                self.ystns,
                self.folder_dir)
        
        response.calc_GF(self.npts,self.dt,self.variable,self.folder_dir,self.stn_file_name,T0=self.T0)
    
    def get_responses_all_patches(self):
        self.source()
        
        for p in range(len(self.Xsrcs)):
            self.get_response_single_patch(p,self.dip[p],self.depth[p],round(self.Tr[p],1),self.Xsrcs[p])


    def time_shift(self):
        self.t0_min = min(self.t0)
        self.id_t0_min = np.argmin(abs(self.t0_min - self.t0))
        
        self.dt_offset = self.t0  - self.t0_min

    def time_shift_waveform(self,comp,receiverID,direction):
        self.mean_model_reader()
        self.time_shift()
        self.folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir='job_'+self.job_id)
        self.slip_dir_folder =  os.path.join(self.folder_dir,direction)  
        self.comp_folder = os.path.join(self.slip_dir_folder,comp)
        comp_folder_files = os.listdir(self.comp_folder)
        # initialize Stream object to store all response along comp at receiver
        strm = obs.Stream()

        fmt = re.compile('p\d{4}_%s_%s.\w+'%(receiverID,comp))

        # accessing files with minimum initial rise time to which all traces are shifted 

        file_t0_min ='p%.4d_%s_%s.SAC'%(self.id_t0_min,receiverID,comp) 
        file_t0_min_dir = os.path.join(self.comp_folder,file_t0_min)
        st_t0_min = obs.read(file_t0_min_dir)
        global_tshift_left = 50.0
        global_tshift_right = 500.0
        initial_start = st_t0_min[0].stats.starttime
        starttime_ref = initial_start - global_tshift_left  
        endtime_ref = st_t0_min[0].stats.endtime + global_tshift_right
        st_t0_min[0].trim(starttime=starttime_ref,endtime=endtime_ref,pad=True,fill_value=0)
        st_t0_min[0].stats.station = file_t0_min[:12].replace("_","")
        comp_folder_files.remove(file_t0_min)
        strm +=st_t0_min
        for file in comp_folder_files:
            check = re.findall(fmt,file)
            if check:
                id = int(file[1:5])
                file_dir = os.path.join(self.comp_folder,file)
                st = obs.read(file_dir)
                st[0].stats.starttime = initial_start + self.T0 + round(self.dt_offset[id],1)
                st[0].stats.station = file[:12].replace("_","")
                # fit size of reference trace with minimun initial time rise 
                st[0].trim(starttime= starttime_ref, endtime=endtime_ref,pad=True,fill_value=0)
                strm += st
        
        return strm

    def time_shift_both(self,comp,receiverID,extra=''):

        as_strm = self.time_shift_waveform(comp,receiverID,'along_strike')
        ad_strm = self.time_shift_waveform(comp,receiverID,'along_dip')

        full_as_dir = os.path.join(self.comp_folder,'full')
        full_ad_dir = os.path.join(self.comp_folder,'full')
        os.makedirs(full_as_dir,exist_ok=True)
        os.makedirs(full_ad_dir,exist_ok=True)

        as_dir = os.path.join(full_as_dir,f"{self.name}_{self.samples}_full_r_{receiverID}_{comp}_{extra}_as.mseed")
        ad_dir = os.path.join(full_ad_dir,f"{self.name}_{self.samples}_full_r_{receiverID}_{comp}_{extra}_ad.mseed")

        as_strm.write(as_dir)  
        ad_strm.write(ad_dir)  

    

    



        
        



# Works well; need to clean it up 
# gorkha = DynamicGF(pwd,'Gorkha',10,(9,18),100,'A',u,'model_test',np.array([1,2]),np.array([-1,-2]))
"""
import sys 
patch_index = int(sys.argv[1])
python test_wavemod.py 1

"""


pwd = os.getcwd()



patch= 29
ncols,nrows = 24,9
factor = 4
offsetx =  (ncols//factor)*patch
offsety =  (nrows//factor)*patch

# xstn = np.arange(patch/2 - offsetx, ncols*patch + offsetx,patch)
# ystn = np.arange(-offsety - (nrows-1/2)*patch, offsety,patch)

xstn = np.linspace(patch/2 - offsetx, ncols*patch,5)
ystn = np.linspace(-offsety - (nrows-1/2)*patch,0,5)
X,Y = np.meshgrid(xstn,ystn)

stations = np.stack((X.flatten(),Y.flatten()),axis=1)



def send_job(i):
    x_stn = np.array([stations[i][0]])
    y_stn = np.array([stations[i][1]])
    tohoku = DynamicGF(pwd,'Tohoku',29,(9,24),100,'D','model_test',x_stn,y_stn,str(i))
    tohoku.time_shift_both('Z','0000',extra=f'j_{str(i)}') 
    tohoku.time_shift_both('N','0000',extra=f'j_{str(i)}')
    tohoku.time_shift_both('E','0000',extra=f'j_{str(i)}')

def parallel_jobs():
    pool = mp.Pool()
    pool.map(send_job,range(len(stations)))

if __name__=='__main__':
    parallel_jobs()

