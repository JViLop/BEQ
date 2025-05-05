# -*- coding: utf-8 -*-
"""
Created on Sat Sep  7 14:37:46 2024

@author: joanv
"""

'''
D. Melgar 02/2014

Forward modeling routines
'''

#from numba import jit




from mudpy import runslip_ensemble,forward_ensemble,forward,runslip

from obspy.core import UTCDateTime
from numpy import array
import sys
import os
########                            GLOBALS                             ########
working_dir = os.getcwd()
name_eq = str(sys.argv[1])

nsamples = 100
n = int(sys.argv[2])
# home = f'/home/josevilo/dynamic/MudPy/{name_eq}/{nsamples}_samples/'
home= os.path.join(working_dir,f'Dynamic_Simulations/{name_eq}/{nsamples}_samples/')
gf_project_name = '1'
project_name=f'{n}'
run_name='fwd'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
# =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
ncpus=60
hot_start=0  #Start at a certain subfault number
static=0  #=1 computes static GFs only, =0 computes the complete waveform
tsunami=False
model_name=f'{name_eq}.mod'   #Velocity model
rupture_name=f'{name_eq}.rupt'   #Rupture model, not needed for inversion
fault_name=f'{name_eq}.fault'    #Fault geometry
station_file=f'{name_eq}.sta'   #Station distribution
GF_list=f'{name_eq}.gflist'#What GFs are to be computed for each station
G_from_file=False
G_name=f'{name_eq}'
NFFT=256 ; dt=1  #Time parameters
dk=1 ; pmin=0 ; pmax=1 ; kmax=12   #fk integration parameters
custom_stf=None
################################################################################

############                 Synthetics parameters               ###############
models = {'Tohoku':{'time_epi':UTCDateTime('2011-03-11T05:46:24'),'epicenter':array([38.322,142.369,29]),'Mw':9.0},
          'Iquique':{'time_epi':UTCDateTime('2014-04-01T23:46:45'),'epicenter':array([-19.572,-70.908,10.0]),'Mw':8.1},
          'Illapel':{'time_epi':UTCDateTime('2015-09-16T22:54:31'),'epicenter':array([-31.25,72,23.0]),'Mw':8.3},
          'Gorkha':{'time_epi':UTCDateTime('2015-04-15T06:11:26'),'epicenter':array([28.1473,84.7079,10.0]),'Mw':7.8},
          'Pedernales':{'time_epi':UTCDateTime('2016-04-16T23:58:36'),'epicenter':array([0.31,-80.15,19.6]),'Mw':7.8}}

############                 Synthetics parameters               ###############
time_epi=models[name_eq]['time_epi']
epicenter=models[name_eq]['epicenter'] 
Mw = models[name_eq]['Mw']
resample=None #Resample synthetics to this rate (in Hz)
integrate=1 #=0 produces velocities, =1 makes displacements
beta=0 #Rake offset, usually a good idea to keep at zero
num_windows=1
#rupture_speed=3.0  #Only necessary if onset times are not identified in rupt file
stf_type='triangle'
source_time_function =stf_type
################################################################################

    
    
zeta=0.2,
stf_falloff_rate=4.0
hot_start=0
   
from numpy import genfromtxt,array,loadtxt
import datetime
from obspy import Stream,read
from mpi4py import MPI

print('Solving for kinematic problem(s)')
#Time for log file
now=datetime.datetime.now()
now=now.strftime('%b-%d-%H%M')
vord = 'disp'
rupture_list=None
#load source names
if rupture_list==None:
    #all_sources=array([home+project_name+'/forward_models/'+rupture_name])   
    all_sources=array([rupture_name])  
else:
    all_sources=genfromtxt(home+gf_project_name+'/data/'+rupture_list,dtype='U')

print('... loading all synthetics into memory')

station_file=home+gf_project_name+'/data/station_info/'+GF_list
staname=genfromtxt(station_file,dtype="U",usecols=0)
Nsta=len(staname)
#Load fault model
source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
Nfaults=source.shape[0] #Number of subfaults
kindex=0


if __name__ =='__main__':
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
    
    Ess=Stream()
    Nss=Stream()
    Zss=Stream()
    Eds=Stream()
    Nds=Stream()
    Zds=Stream()
    for ksta in range(data[0],data[1]):
        print('Reading green functions for station #'+str(ksta+1)+' of '+str(Nsta))
        for kfault in range(Nfaults):
            print(f'MPI: processor {rank} -> reading green functions for station # '+str(ksta+1)+' of '+str(Nsta) + 'for fault # '+str(kfault))
            #Get subfault GF directory
            nsub='sub'+str(int(source[kfault,0])).rjust(4,'0')
            nfault='subfault'+str(int(source[kfault,0])).rjust(4,'0')
            strdepth='%.4f' % source[kfault,3]
            syn_path=home+project_name+'/GFs/dynamic/'+model_name+'_'+strdepth+'.'+nsub+'/'
            #Get synthetics
    
            Ess+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
            os.remove(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
            Nss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
            os.remove(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
            Zss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
            os.remove(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
            Eds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
            os.remove(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
            Nds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
            os.remove(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
            Zds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
            os.remove(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
            kindex+=1
        
    Ess = comm.reduce(Ess,op=MPI.SUM, root=0)
    Nss = comm.reduce(Nss,op=MPI.SUM, root=0)
    Zss = comm.reduce(Zss,op=MPI.SUM, root=0)
    Eds = comm.reduce(Eds,op=MPI.SUM, root=0)
    Nds = comm.reduce(Nds,op=MPI.SUM, root=0)
    Zds = comm.reduce(Zds,op=MPI.SUM, root=0)     
    
    
    if rank==0:
        
        
        print('Writting synthetics to miniSEED, hang on this might take a minute or two.')
        #Ess.write(home+project_name+'/GFs/matrices/'+G_name+'.Ess.'+vord+'.mseed',format='MSEED')
        #Nss.write(home+project_name+'/GFs/matrices/'+G_name+'.Nss.'+vord+'.mseed',format='MSEED')
        #Zss.write(home+project_name+'/GFs/matrices/'+G_name+'.Zss.'+vord+'.mseed',format='MSEED')
        #Eds.write(home+project_name+'/GFs/matrices/'+G_name+'.Eds.'+vord+'.mseed',format='MSEED')
        #Nds.write(home+project_name+'/GFs/matrices/'+G_name+'.Nds.'+vord+'.mseed',format='MSEED')
        #Zds.write(home+project_name+'/GFs/matrices/'+G_name+'.Zds.'+vord+'.mseed',format='MSEED')
       	
    
    
    
    #Now loop over rupture models
        for ksource in range(hot_start,len(all_sources)):
            print('... solving for source '+str(ksource)+' of '+str(len(all_sources)))
            rupture_name=all_sources[ksource]
            print(rupture_name)
            
            if rupture_list!=None:
                #Get epicentral time
                epicenter,time_epi=forward_ensemble.read_fakequakes_hypo_time(home,gf_project_name,rupture_name)
                forwardBool=False
            else:
                forwardBool=True #This controls where we look for the rupture file
            
            # Put in matrix
            m,G=forward_ensemble.get_fakequakes_G_and_m(Nss,Ess,Zss,Nds,Eds,Zds,home,project_name,gf_project_name,rupture_name,time_epi,GF_list,epicenter,NFFT,source_time_function,stf_falloff_rate,zeta,forward=forwardBool)
            # Solve
            waveforms=G.dot(m)
            #Write output
            forward_ensemble.write_fakequakes_waveforms(home,project_name,gf_project_name,rupture_name,waveforms,GF_list,NFFT,time_epi,dt,loc_epi = epicenter ,Mw = Mw)
            
    

            
                            
