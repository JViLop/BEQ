'''
Diego Melgar, 09/2015

Parameter file
Project: Coquimbo
Comment: 
'''

from mudpy import runslip, forward,runslip_ensemble,forward_ensemble
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
make_green=int(sys.argv[3])#Compute GFs
make_synthetics=int(sys.argv[4]) #Compute synthetics for a given model at given stations
solve=int(sys.argv[5]) # =1 solves forward problem or runs inverse calculation, =0 does nothing
###############################################################################

###############            Green function parameters               #############
ncpus=12  #
hot_start=0  #Start at a certain subfault number
static=0  #=1 computes static GFs only, =0 computes the complete waveform
tsunami=False
model_name=f'{name_eq}.mod'   #Velocity model
rupture_name=f'{name_eq}.rupt'   #Rupture model, not needed for inversion
fault_name=f'{name_eq}.fault'    #Fault geometry
station_file=f'{name_eq}.sta'   #Station distribution
GF_list=f'{name_eq}.gflist'#What GFs are to be computed for each station
G_from_file=True
G_name=f'{name_eq}'
NFFT=256 ; dt=1  #Time parameters
dk=1 ; pmin=0 ; pmax=1 ; kmax=12   #fk integration parameters
custom_stf=None
################################################################################

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
################################################################################


#Initalize project folders
if init==1:
    runslip.init(home,project_name)

# Run green functions          
if make_green==1: 
    if ncpus<2:
        runslip.make_green(home,project_name,station_file,fault_name,model_name,
            dt,NFFT,static,tsunami,hot_start,dk,pmin,pmax,kmax)
    else:
        runslip_ensemble.make_parallel_green(home,project_name,station_file,fault_name,
            model_name,dt,NFFT,static,tsunami,hot_start,dk,pmin,pmax,kmax,ncpus)


#Now make synthetics for source/station pairs
if make_synthetics==1:
    if ncpus<2: 
        runslip.make_synthetics(home,project_name,station_file,fault_name,model_name,integrate,
                static,tsunami,beta,hot_start,time_epi)
    else:

        runslip_ensemble.make_parallel_synthetics(home,project_name,gf_project_name,station_file,fault_name,model_name,integrate,
                static,tsunami,beta,hot_start,time_epi,ncpus,custom_stf=None)

            
    
#Run forward comptuation or solve for inverse problem
if solve==1:
    if static==0: #Forward problem (full waveforms)
        forward_ensemble.waveforms_fakequakes(home,project_name,gf_project_name,fault_name,None,GF_list,
                model_name,run_name,dt,NFFT,G_from_file,G_name,source_time_function=stf_type,
                stf_falloff_rate=4.0,rupture_name=rupture_name,epicenter=epicenter,time_epi=time_epi,Mw= Mw)
    if static==1: #Forward problem (coseismics)
        forward.coseismics_matrix(home,project_name,rupture_name,station_file,G_from_file,model_name,G_name)
        
        

    
    
    



