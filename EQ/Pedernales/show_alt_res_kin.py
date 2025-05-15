#!/usr/bin/env python

# Import externals
import sys
import copy
import shutil
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import h5py

# Import personal libraries
import csi.RectangularPatchesKin as rectFault # Use only for kinematic modeling
import csi.multifaultsolve as multiflt
import csi.geodeticplot as geoplt
import csi.insar as ir
import csi.gps   as gr
import csi.tsunami    as tsunami
import csi.seismic    as seis
import csi.faultpostproc as faultpostproc

import AltarExplore as alt

# Import Wavemod
import wavemod as wm
# Load Arguments
from Arguments import *
from filt import *

n_ramp_param = 9
utmzone = 17
#------------------------------------------------------------------
#------------------------------------------------------------------
# Load Arguments
#from Arguments import *
Compute_model = True     # Compute the mean/median/mode from posterior samples?
BETA_STEP     = 53        # Beta step id
res_dir       = 'results' # Path to results directory
comp_GFs      = False     # Re-compute Green's functions from EDKS kernels?
comp_KinGFs   = False
comp_bigG     = False
post_plot     = 'mean'    # Which model do you want to plot (mean/median/mode)

comp_HD_GFs = False
HD_GFdir    = './GFs_HD'

#----------------------------------------------------------------
#---------------------------- FAULT  ----------------------------
# Initialize faults
fault = rectFault('Megathrust', lon0=lon0, lat0=lat0)
fault.readPatchesFromFile(faultGeometry, readpatchindex=False, donotreadslip=True)
fault.setTrace(delta_depth=0.1)
fault.initializeslip()
fault.computeEquivRectangle()        
fault.computeArea()

fault.initializekinmodel()

# Set hypocenter
fault.setHypoXY(lon_hypo,lat_hypo,UTM=False) # Hypocenter (for kinematic modeling)

# Set mapping
fault.setFaultMap(Nstrike,Ndip,leading='dip',check_depth=True)

# Set Mu
fault.setMu(earth_model)

# Define Sub Grid
fault.buildSubGrid(Npt,Npt)    

print 'NUMPATCH : ', len(fault.patch)



#----------------------------------------------------------------
#------------------------ SEISMIC DATA --------------------------
# Init waveform engine and seismic data
print("Init waveform engine")
sm_sac_files = []
   
with open(sac_lst_sm) as f:
    for l in f:
        items = l.strip().split()        
        sm_sac_files.append(join(seis_dir_sm,items[0]))

gps_sac_files = []

with open(sac_lst_gps) as f:
    for l in f:
        items = l.strip().split()        
        gps_sac_files.append(join(seis_dir_gps,items[0]))

sm_data = seis('Strong_motion',lon0=lon0, lat0=lat0)
sm_data.readSac(sm_sac_files)
sm_data.readCdFromBinaryFile('Cd.Strong_motion.bin')
hrgps_data = seis('HRGPS_kin',lon0=lon0, lat0=lat0)
hrgps_data.readSac(gps_sac_files)
hrgps_data.readCdFromBinaryFile('Cd.HRGPS.bin')

seismic_data = [hrgps_data,sm_data]        ############# TO BE CORRECTED

# Waveform engine
waveDB = wm.WaveDB(gf_db_dir,gf_db_scale)


# Compute/Load Kinematic Green's functions                
for i in range(len(seismic_data)):
    data = seismic_data[i]    
    # Compute/Load/Save GFs
    if comp_KinGFs: 
        # Compute GFs
        print(data.name)
        fault.buildKinGFsFromDB(data,waveDB,1.,rake1, rake_key=rakes_key[0], filter_coef=sos, differentiate=False)
        fault.buildKinGFsFromDB(data,waveDB,1.,rake2, rake_key=rakes_key[1], filter_coef=sos, differentiate=False)
        # Save GFs
        fault.saveKinGF(data,outputDir=GFdir,rmdir=False)    
    elif comp_bigG:
        # Load GFs
        fault.loadKinGF(data,rakes_key,inputDir=GFdir)

# Build Big G matrix
print("Assemble bigG and bigD")
bigDfile = os.path.join(GFdir,'kinematicG.data.bin')
bigCdfile = os.path.join(GFdir,'kinematicG.Cd.bin')
bigGfile = os.path.join(GFdir,'kinematicG.gf.bin')
fault.setBigDmap(seismic_data)
if comp_bigG or comp_KinGFs:
    fastSweep = wm.FastSweep()
    fault.buildBigGD(fastSweep,seismic_data,rakes_key,-1.,Ntriangles,Dtriangles,dtype='float32')
    fault.saveBigGD(bigDfile,bigGfile,dtype='float32')
else:
    fault.loadBigGD(bigDfile,bigGfile,dtype='float32')



#--------------------------------------------------------------------
#--------------------------------------------------------------------
# Get the solution
Np = len(fault.patch)
if Compute_model:    
    post_mode   = np.zeros((4*Np+n_ramp_param+2,),dtype='double')
    post_mean   = np.zeros((4*Np+n_ramp_param+2,),dtype='double')
    post_median = np.zeros((4*Np+n_ramp_param+2,),dtype='double')
    print('n-param',2*Np+n_ramp_param)
    fault.Cm = np.zeros((2*Np,3),dtype='double')

    # Read hdf5 file
    f = h5py.File(join(res_dir,'step_%03d.h5'%(BETA_STEP)),'r') # hdf5 file object
    assert f['Sample Set'].shape[1] == 4*Np + n_ramp_param +2, 'Incorrect number of parameters'
    samples = f['Sample Set'].value # Samples
    f.close()
    
    # Compute Marginals mean/median/mode
    for p in range(Np): #Slip model

        RotAngle2 = RotAngle
        xc, yc, zc, width, length, strike, dip = fault.getpatchgeometry(p, center=True) 
                
        # Convert angle in radians
        RotAngle2 *= ((np.pi) / 180.)
        rotation = np.arctan2(np.tan(strike) - np.tan(RotAngle2), 
                            np.cos(dip)*(1.+np.tan(RotAngle2)*np.tan(strike)))

        # If RotAngle within ]90, 270], change rotation
        if RotAngle > 90. and RotAngle<=270.:
            rotation += np.pi

        rp = alt.parameter('strike',samples[:,p]   ,bins=100,dtype=np.float32,color='red')
        ar = alt.parameter('dip'   ,samples[:,p+Np],bins=100,dtype=np.float32,color='red')

        ss = copy.deepcopy(rp)
        ds = copy.deepcopy(ar)
        ss.values = ar.values*np.cos(rotation) - rp.values*np.sin(rotation)
        ds.values = ar.values*np.sin(rotation) + rp.values*np.cos(rotation)
        
        ssi = np.where(ss.marginal.h==ss.marginal.h.max())[0][0]
        dsi = np.where(ds.marginal.h==ds.marginal.h.max())[0][0]
        # Mode
        post_mode[p]    = ss.marginal.x[ssi]
        post_mode[p+Np] = ds.marginal.x[dsi]
        # Mean
        post_mean[p]    = ss.mean()
        post_mean[p+Np] = ds.mean()
        # Median
        post_median[p]    = ss.median()
        post_median[p+Np] = ds.median()
        
        M   = np.array([ss.values,ds.values])    
        cov = np.cov(M)
        fault.Cm[p,:] = np.array([cov[0,0],cov[1,1],cov[0,1]])    
            
    for i in range(2*Np,2*Np+n_ramp_param): #Orbital parameters
        p = alt.parameter('orbital',samples[:,i],bins=100,dtype=np.float32)
        pi = np.where(p.marginal.h==p.marginal.h.max())[0][0]
        post_mode[i]   = p.marginal.x[pi]
        post_mean[i]   = p.mean()
        post_median[i] = p.median()
        #p.marginal.show()

    for i in range(2*Np+n_ramp_param,4*Np+n_ramp_param+2): #Kinematic parameters
        p = alt.parameter('kinematic',samples[:,i],bins=100,dtype=np.float32)
        pi = np.where(p.marginal.h==p.marginal.h.max())[0][0]
        post_mode[i]   = p.marginal.x[pi]
        post_mean[i]   = p.mean()
        post_median[i] = p.median()        

    # Save model
    np.savetxt(join(res_dir,'kinematicG.final_mode'),post_mode)
    np.savetxt(join(res_dir,'kinematicG.final_mean'),post_mean)
    np.savetxt(join(res_dir,'kinematicG.final_median'),post_median)
    np.savetxt(join(res_dir,'kinematicG.Cm'),fault.Cm)
    
# Load model
post     = np.loadtxt(join(res_dir,'kinematicG.final_'+post_plot))

if os.path.exists(join(res_dir,'kinematicG.Cm')):
    fault.Cm = np.loadtxt(join(res_dir,'kinematicG.Cm'))
else:
    Np = len(fault.patch)
    fault.Cm = np.zeros((2*Np,3),dtype='double')        

#----------------------------------------------------
#----------------- Plot solution -------------------

# solve the least square pb
#slv.mpost=post.copy()
#slv.distributem()
#slv_HD.mpost=post.copy()
#slv_HD.distributem()

#---------------------------------------------------
#---------- Compute magnitude/mechanisms -----------
if Compute_model:
    fpp = faultpostproc(fault.name,fault,Mu=fault.mu,utmzone=utmzone,samplesh5=join(res_dir,'step_%03d.h5'%(BETA_STEP)))
    fpp.h5_init()
    fpp.h5_finalize()
    fpp.computeMomentTensor()
    Mo=fpp.computeScalarMoment()
    Mw=fpp.computeMagnitude(plotHist=res_dir)

    areath,Moth,slipth,StressDrop = fpp.stressdrop(threshold=0.35,threshold_rand=False,return_Area_Mo_Slip=True)
    MaxSlip=np.sqrt(samples[:,:Np]*samples[:,:Np]+samples[:,Np:2*Np]*samples[:,Np:2*Np]).max(axis=1)
    
    plt.figure()
    plt.hist(StressDrop/1e6,bins=40,color='g',normed=True)
    plt.xlabel('\n$\Delta\sigma,$ MPa',fontsize=14)
    plt.savefig(join(res_dir,'stressdrop.pdf'))

    plt.figure()
    plt.hist(slipth,bins=100,color='g',normed=True)
    plt.xlabel('\nMean Slip, m',fontsize=14)
    plt.savefig(join(res_dir,'meanslip.pdf'))
    
    plt.figure()
    plt.hist(MaxSlip,bins=100,color='g',normed=True)
    plt.xlabel('Maximum Slip, m',fontsize=14)
    plt.savefig(join(res_dir,'maxslip.pdf'))

    plt.figure()
    plt.hist(Mw,bins=100,color='g',normed=True)
    plt.xlabel('Mw',fontsize=14)
    plt.savefig(join(res_dir,'Mw.pdf'))

    plt.figure()
    plt.hist(Mo,bins=100,color='g',normed=True)
    plt.xlabel('Seismic moment, N.m',fontsize=14)
    plt.savefig(join(res_dir,'Mo.pdf') )   
    
    plt.show()
    print('Magnitude Mw=%.2f'%(Mw.mean()))
    print('Moment M0=%.3e N-m'%(Mo.mean()))
    
else:
    fpp = faultpostproc(fault.name,fault,Mu=fault.mu,utmzone=17)
    fpp.computeMomentTensor()
    Mo=fpp.computeScalarMoment()
    Mw=fpp.computeMagnitude()
    print('Magnitude Mw=%.2f'%(Mw))
    print('Moment M0=%.3e N-m'%(Mo))

#---------------------------------------------------------------
#-------------------- COMPUTE SYNTHETICS -----------------------

## Build synthetics for kinematic data ##
# Assign fault parameters
fault.slip[:,0] = post[:Np]
fault.slip[:,1] = post[Np:2*Np]
fault.tr = post[2*Np+n_ramp_param:3*Np+n_ramp_param]
fault.vr = post[3*Np+n_ramp_param:4*Np+n_ramp_param]
h_strike = post[4*Np+n_ramp_param]
h_dip    = post[4*Np+n_ramp_param+1]
fault.setHypoOnFault(h_strike,h_dip)
fpp2 = faultpostproc(fault.name,fault,Mu=fault.mu,utmzone=utmzone)
fpp2.computeMomentTensor()

# Fast sweeping
eik_solver = wm.FastSweep()
eik_solver.setGridFromFault(fault,pSize/Nmesh) # MEME QUE ALTAR !!!!!!!!
eik_solver.fastSweep()

# G x M
G = fault.bigG
D = fault.bigD
m = np.zeros((G.shape[1],))
STFMT = np.zeros((Ntriangles,3,3))
rup_time = []
for p in range(len(fault.patch)):
    # Location at the patch center
    p_x, p_y, p_z, p_width, p_length, p_strike, p_dip = fault.getpatchgeometry(p,center=True)
    dip_c, strike_c = fault.getHypoToCenter(p,True)
    # Grid location
    grid_size_dip = p_length/Npt
    grid_size_strike = p_length/Npt
    grid_strike = strike_c+np.arange(0.5*grid_size_strike,p_length,grid_size_strike) - p_length/2.
    grid_dip    = dip_c+np.arange(0.5*grid_size_dip   ,p_width ,grid_size_dip   ) - p_width/2.
    time = np.arange(Ntriangles)*Dtriangles#+Dtriangles
    T  = np.zeros(time.shape)
    Tr2 = fault.tr[p]/2.
    rup_time.append(eik_solver.getT0([dip_c],[strike_c])[0])
    for i in range(Npt):
        for j in range(Npt):
            t = eik_solver.getT0([grid_dip[i]],[grid_strike[j]])[0]
            tc = t+Tr2
            ti = np.where(np.abs(time-tc)<Tr2)[0]            
            T[ti] += (1/Tr2 - np.abs(time[ti]-tc)/(Tr2*Tr2))*Dtriangles
    for nt in range(Ntriangles):
        m[2*nt*Np+p]     += T[nt] * fault.slip[p,0]/float(Npt*Npt)
        m[(2*nt+1)*Np+p] += T[nt] * fault.slip[p,1]/float(Npt*Npt)
        STFMT[nt] += T[nt] * fpp2.MTs[p]/float(Npt*Npt)

STF = np.zeros((Ntriangles,))
for nt in range(Ntriangles):
    STF[nt] = np.sqrt(0.5 * np.sum(STFMT[nt]**2,axis=(0,1)))

rup_time = np.array(rup_time)
print(G.shape,D.shape,m.shape,2*Np*Ntriangles)
seismic_predictions = np.dot(G,m)

plt.figure()
plt.plot(D)
plt.plot(seismic_predictions)

#---------------------------------------------------------------
#----------------------PREPARE GMT FILES -----------------------
f = open(fault.name+'.hypo','wt')
f.write('%.3f, %.3f, %.3f\n'%(fault.hypo_lon,fault.hypo_lat,fault.hypo_z))
f.close()
    
fault.writePatches2File(fault.name+'.fault',add_slip='total')
fault.writeSlipDirection2File(fault.name+'.slipdir',scale='total',factor=3.0,ellipse=True,nsigma=2.4477)

sslip = fault.slip[:,0].copy()
fault.slip[:,0] = rup_time.copy()
fault.writePatches2File(fault.name+'.trup',add_slip='strikeslip')

fault.slip[:,0] = fault.vr.copy()
fault.writePatches2File(fault.name+'.vr',add_slip='strikeslip')

fault.slip[:,0] = fault.tr.copy()
fault.writePatches2File(fault.name+'.tr',add_slip='strikeslip')

fault.slip[:,0] = sslip.copy()

fid = open(fault.name+'.fault_xyz','w')
for p in range(len(fault.patch)):
    slp = np.sqrt(fault.slip[p,0]**2 + fault.slip[p,1]**2)
    p_x, p_y, p_z, width, length, strike_rad, dip_rad = fault.getpatchgeometry(p,center=True)
    lon,lat = fault.xy2ll(p_x,p_y)
    fid.write('{} {} {} \n'.format(lon,lat,slp))
fid.close()

fid = open(fault.name+'.trup_xyz','w')
for p in range(len(fault.patch)):
    slp = rup_time[p]
    p_x, p_y, p_z, width, length, strike_rad, dip_rad = fault.getpatchgeometry(p,center=True)
    lon,lat = fault.xy2ll(p_x,p_y)
    fid.write('{} {} {} \n'.format(lon,lat,slp))
fid.close()


fid = open(fault.name+'.tr_xyz','w')
for p in range(len(fault.patch)):
    slp = fault.tr[p]
    p_x, p_y, p_z, width, length, strike_rad, dip_rad = fault.getpatchgeometry(p,center=True)
    lon,lat = fault.xy2ll(p_x,p_y)
    fid.write('{} {} {} \n'.format(lon,lat,slp))
fid.close()


#---------------------------------------------------------------
#----------------------STRESS DROP I FILE ----------------------

#---------------------------------------------------------
#-------------------- SHOW RESULTS -----------------------
    

