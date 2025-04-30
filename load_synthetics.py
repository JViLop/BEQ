# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 21:31:02 2024

@author: joanv
"""


from obspy.core import UTCDateTime
from obspy import Stream,read
from numpy import genfromtxt,loadtxt
from mpi4py import MPI
import sys
import os
########                            GLOBALS                             ########
working_dir = os.getcwd()
n = int(sys.argv[1])

nsamples = 100
name_eq = 'Tohoku'
# home = f'/home/josevilo/dynamic/MudPy/{name_eq}/{nsamples}_samples/'
home= os.path.join(working_dir,f'{name_eq}/{nsamples}_samples/')
gf_project_name = '1'
project_name=f'{n}'
run_name='fwd'
################################################################################


#####              What-do-you-want-to-do flags, 1=do, 0=leave be          #####

init=0 #Initalize project
make_green=int(sys.argv[2])#Compute GFs
make_synthetics=int(sys.argv[3]) #Compute synthetics for a given model at given stations
solve=int(sys.argv[4]) # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
ncpus=4
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



vord = 'disp'
station_file=home+gf_project_name+'/data/station_info/'+GF_list
staname=genfromtxt(station_file,dtype="U",usecols=0)
Nsta=len(staname)
#Load fault model
source=loadtxt(home+project_name+'/data/model_info/'+fault_name,ndmin=2)
Nfaults=source.shape[0] #Number of subfaults
kindex=0

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
        Zds=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
        Ess+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.e')
        Nss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.n')
        Zss+=read(syn_path+staname[ksta]+'.'+nfault+'.SS.'+vord+'.z')
        Eds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.e')
        Nds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.n')
        Zds+=read(syn_path+staname[ksta]+'.'+nfault+'.DS.'+vord+'.z')
        kindex+=1
        


if rank==0:
    
    Ess = comm.reduce(Ess,op=MPI.SUM, root=0)
    Nss = comm.reduce(Nss,op=MPI.SUM, root=0)
    Zss = comm.reduce(Zss,op=MPI.SUM, root=0)
    Eds = comm.reduce(Eds,op=MPI.SUM, root=0)
    Nds = comm.reduce(Nds,op=MPI.SUM, root=0)
    Zds = comm.reduce(Zds,op=MPI.SUM, root=0)
    print('Writting synthetics to miniSEED, hang on this might take a minute or two.')
    Ess.write(home+project_name+'/GFs/matrices/'+G_name+'.Ess.'+vord+'.mseed',format='MSEED')
    Nss.write(home+project_name+'/GFs/matrices/'+G_name+'.Nss.'+vord+'.mseed',format='MSEED')
    Zss.write(home+project_name+'/GFs/matrices/'+G_name+'.Zss.'+vord+'.mseed',format='MSEED')
    Eds.write(home+project_name+'/GFs/matrices/'+G_name+'.Eds.'+vord+'.mseed',format='MSEED')
    Nds.write(home+project_name+'/GFs/matrices/'+G_name+'.Nds.'+vord+'.mseed',format='MSEED')
    Zds.write(home+project_name+'/GFs/matrices/'+G_name+'.Zds.'+vord+'.mseed',format='MSEED')

