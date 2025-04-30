# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 21:25:32 2024

@author: joanv
"""

from mpi4py import MPI
import time
import math
import os
import numpy as np
from obspy import read
import obspy 

home=os.path.join('/home/josevilo','MudPy-1.3/examples')
project_name = 'Tohoku'
model_name='tohoku.mod'   #Velocity model
rupture_name='tohoku.rupt'   #Rupture model, not needed for inversion
fault_name='tohoku.fault'    #Fault geometry
station_file='tohoku.sta'   #Station distribution
GF_list='gps.gflist'
vord = 'disp'

G_name = 'Tohoku_mpi'
station_file= home +'/'+project_name+'/data/station_info/'+GF_list

staname= np.genfromtxt(station_file,dtype="U",usecols=0)
Nsta=len(staname)
#Load fault model
source=np.loadtxt(home+'/'+project_name+'/data/model_info/'+fault_name,ndmin=2)
Nfaults=source.shape[0] #Number of subfaults



if __name__=='__main__':
    t0 = time.time()
    
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nprocs = comm.Get_size()
    
    # number of integration steps
    nsteps = Nsta 
    
    
    if rank == 0:
        # determine the size of each sub-task
        ave, res = divmod(nsteps, nprocs)
        counts = [ave + 1 if p < res else ave for p in range(nprocs)]
    
        # determine the starting and ending indices of each sub-task
        starts = [sum(counts[:p]) for p in range(nprocs)]
        ends = [sum(counts[:p+1]) for p in range(nprocs)]
    
        # save the starting and ending indices in data  
        data = [(starts[p], ends[p]) for p in range(nprocs)]
    else:
        data = None
    
    data = comm.scatter(data, root=0)
    
    # compute partial contribution to pi on each process
    
    Ess=obspy.Stream()
    Nss=obspy.Stream()
    Zss=obspy.Stream()
    Eds=obspy.Stream()
    Nds=obspy.Stream()
    Zds=obspy.Stream()
    for ksta in range(data[0],data[1]):
            print('... ... reading green functions for station #'+str(ksta+1)+' of '+str(Nsta))
            
    
            
            for kfault in range(Nfaults):
                print("... reading for station:"+str(ksta) + " -> source:" + str(kfault))
                nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
                nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
                strdepth='%.4f' % source[kfault,3]
                syn_path=home+'/'+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
                #Get synthetics
                Ess+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
                Nss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
                Zss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
                Eds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
                Nds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
                Zds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
    	
    '''
    Ess = comm.reduce(Ess, root=0)
    Nss = comm.reduce(Nss, root=0)
    Zss = comm.reduce(Zss, root=0)
    Eds = comm.reduce(Eds, root=0)
    Nds = comm.reduce(Nds, root=0)
    Zds = comm.reduce(Zds, root=0)
    
    '''
    Ess = comm.reduce(Ess,op=MPI.SUM, root=0)
    Nss = comm.reduce(Nss,op=MPI.SUM, root=0)
    Zss = comm.reduce(Zss,op=MPI.SUM, root=0)
    Eds = comm.reduce(Eds,op=MPI.SUM, root=0)
    Nds = comm.reduce(Nds,op=MPI.SUM, root=0)
    Zds = comm.reduce(Zds,op=MPI.SUM, root=0)
    if rank==0:
        
        
        
        Ess.write(home+'/'+project_name+'/GFs/matrices/'+G_name+'.Ess.'+vord+'.mseed',format='MSEED')
        Nss.write(home+'/'+project_name+'/GFs/matrices/'+G_name+'.Nss.'+vord+'.mseed',format='MSEED')
        Zss.write(home+'/'+project_name+'/GFs/matrices/'+G_name+'.Zss.'+vord+'.mseed',format='MSEED')
        Eds.write(home+'/'+project_name+'/GFs/matrices/'+G_name+'.Eds.'+vord+'.mseed',format='MSEED')
        Nds.write(home+'/'+project_name+'/GFs/matrices/'+G_name+'.Nds.'+vord+'.mseed',format='MSEED')
        Zds.write(home+'/'+project_name+'/GFs/matrices/'+G_name+'.Zds.'+vord+'.mseed',format='MSEED')
       
    
        tf = time.time()
        print(str(tf-t0)+'s')
        