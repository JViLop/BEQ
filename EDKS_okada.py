# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 17:27:17 2024

@author: joanv
"""

 
# adding csi functions to the system path
import numpy as np
import sys
import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# Okada
import okada4py as ok92

# EDKS
from layered_disloc import layered_disloc



# NOTE: In this convention, ss(+) means left-lateral, and ss(-) means right-lateral

#--------------------------------------------------
# Check inputs
def ArraySizes(*args):
    '''
    Only requirement is that each arguments has the same size and can be converted to a numpy array
    Returns : Numpy arrays
    '''
    
    # Create a list of sizes
    Sizes = []
    Arrays = []

    # Check class
    for arg in args:
        if arg.__class__ in (list, tuple):
            arg = np.array(arg)
        elif arg.__class__ in (float, np.float64, int):
            arg = np.array([arg])
        Arrays.append(arg)
        Sizes.append(arg.shape)
    
    # Assert sizes
    assert (len(np.unique(Sizes))==1), 'The {} provided arrays are not the same size'.format(len(args))

    # All done
    return Arrays

#--------------------------------------------------
# Displacements only
def displacement(xs, ys, zs, xc, yc, zc, width, length, strike, dip, ss, ds, ts, nu=0.25):
    '''
    Returns the displacements at the stations located on (xs, ys, zs) for patches
        with centers on (xc, yc, zc). All arguments can be float, list or array.
    '''

    # Here Mu can be anything. RJ tested it and the displacement is not-sensitive to Mu as it should be.
    # Although, it does not work with Mu = 0.0 GPa... So we take a random value of 30GPa
    mu = 30e9

    # Nu does matter here, and it is by default 0.25

    # Check 
    xs, ys, zs = ArraySizes(xs, ys, zs)
    xc, yc, zc, width, length, strike, dip, ss, ds, ts = ArraySizes(xc, yc, zc, width, length, strike, dip, ss, ds, ts)

    # Normally, StaticInv does angles in Radians
    dip = dip*180./np.pi
    strike = strike*180./np.pi

    # Run okada
    u, d, s, flag, flag2 = ok92.okada92(xs, ys, zs, xc, yc, zc, length, width, dip, strike, ss, ds, ts, mu, nu)

    # Check if things went well
    if not (flag==0).all():
        if not np.where(flag!=0)==[]:
            print(' Error: {}'.format(tuple(np.where(flag!=0))))
            print('Something went wrong in okada4py... You should check...')

    # Reshape the displacement
    u = u.reshape((len(xs), 3))

    # All Done
    return u
        


        


# source


def save_fig(test_id,xstn,ystn,dfs,zsrc,dip,w,l,ss,ds):
    shape_stn = (len(ystn),len(xstn))
    cols = ['Okada','EDKS','Residual']
    fig, axes = plt.subplots(3,3,figsize=(12,10),dpi=600,sharex=True)
    for j, df in enumerate(dfs):
        
        for i,parameter_id in enumerate(df.keys()):
            
            parameter = df[parameter_id].values
            title = f"Test{test_id}: "+ " $z_{source}=$"+f"{zsrc}km" + r" $\delta=$"+f"{dip}° "+r"$A =$"+f"{w*l}km2 "+ "$s$="+f"({ss},{ds})m"
            file_name = f"Test{test_id}"
            # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
            if j!=2:
                im = axes[i][j].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
                fig.colorbar(im, ax=axes[i][j],shrink=0.6,label='{}'.format(f'{parameter_id} (m)'))
            else:
                im = axes[i][j].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k',cmap = 'viridis')
                fig.colorbar(im, ax=axes[i][j],shrink=0.6,label='{}'.format('abs. residual (m)'))
 #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
            
                    
            #fig.colorbar(im, ax=axes[i][j],shrink=0.6,label='{}'.format(f'{parameter_id} (m)'))
        
        
                
            axes[i][j].set_ylabel('Trench-normal distance (km)',fontsize=7)
            axes[i][j].set_xlabel('Along-strike distance (km)',fontsize=7)
            axes[i][j].set_aspect('equal', 'box')
            axes[i][j].tick_params(labelsize=7)

            axes[i][j].set_title(f'{cols[j]}' ,fontweight='bold',fontsize=11)
            #axes[i].invert_yaxis()
    fig.suptitle(title,x=0.5,y=1,fontweight='bold') # uncomment later
    plt.tight_layout()
    main_dir = os.getcwd()
    fig_dir  = os.path.join(main_dir,'Tests3')
    os.makedirs(fig_dir,exist_ok =True)
    fig_name = os.path.join(fig_dir,file_name.replace(" ","_") + '.png')
    fig.savefig(fig_name)
    plt.close()
	

def rel_error(diff,a):
    return np.divide(diff,a,out=np.zeros_like(diff),where=a==0) 

    
    
def calc_static(test_id,dr = 0.5,km=1000.0,w = 10.0,l = 10.0,z_source =10.0,dip_angle = 90.0,Slip = np.array([0.0,50.0])):

    # station 
    xstn = np.arange(-30.0,30.0 + dr,dr)*km
    ystn = np.arange(-30.0,30.0 + dr,dr)*km
    
    Xstn,Ystn = np.meshgrid(xstn,ystn)
    xstn_flat,ystn_flat =  Xstn.flatten(),Ystn.flatten()
    zstn_flat = np.zeros(xstn_flat.shape)
    
    
    xsrc = np.array([0])
    ysrc = np.array([0])
    zsrc = np.array([z_source])*km
    
    width = np.array([w])*km
    length = np.array([l])*km
    strike = np.array([90.0])
    strike_rad = strike*(np.pi/180)
    dip = np.array([dip_angle])
    dip_rad= dip*(np.pi/180)
    
    
    ss = np.array([Slip[0]])
    ds = np.array([Slip[1]])
    ts = np.array([0.0])
    
    nu = 0.3400051746245543
    
    
    
    Displacement = displacement(xstn_flat,ystn_flat,zstn_flat, 
                                        xsrc, 
                                        ysrc, 
                                        zsrc, 
                                        width, 
                                        length, 
                                        strike_rad, 
                                        dip_rad, 
                                        ss, 
                                        ds, 
                                        ts, nu=nu)
    
    E_ok = Displacement[:,0]
    N_ok = Displacement[:,1]
    Z_ok = Displacement[:,2]
    
    df_ok = pd.DataFrame(Displacement,columns=['x','y','z'])
    

    # EDKS initialization
        
    edks = 'poisson_half_space.edks'
    prefix = 'poisson_half_space'
    BIN_EDKS = '${HOME}/edks/bin'
    
    # reshape to form vectors.
    # source coordinates
    xS = xsrc
    yS = ysrc
    
    zS = zsrc
    aS = np.ones_like(dip)*(width[0]*length[0])
    
    SSlip = ss
    DSlip = ds
    
    # calculate the GF's
    rake = 0 * np.ones(SSlip.shape)
    slip = 1 * np.ones(SSlip.shape)
    GFeSS, GFnSS, GFzSS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, xstn_flat, ystn_flat,\
                                  edks, prefix, BIN_EDKS)  
    
    rake = 90* np.ones(SSlip.shape)
    slip = 1 * np.ones(SSlip.shape)
    GFeDS, GFnDS, GFzDS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, xstn_flat, ystn_flat,\
                                  edks, prefix, BIN_EDKS)
    
    
     # compute forward displacement calculation.

    E_edks = np.dot(GFeSS,SSlip) + np.dot(GFeDS, DSlip)
    N_edks = np.dot(GFnSS,SSlip) + np.dot(GFnDS, DSlip)
    Z_edks = np.dot(GFzSS,SSlip) + np.dot(GFzDS, DSlip)
    

    
    dict_edks = {'x':E_edks,'y':N_edks,'z':Z_edks}
    df_edks = pd.DataFrame(dict_edks)

    delta_E =  np.abs(E_ok - E_edks)
    delta_N =   np.abs(N_ok - N_edks)
    delta_Z =   np.abs(Z_ok - Z_edks)
    dict_residual = {'x':delta_E ,'y':delta_N,'z':delta_Z} 
    df_residual = pd.DataFrame(dict_residual)
    
    dfs = [df_ok,df_edks,df_residual]
    save_fig(test_id,xstn/km,ystn/km,dfs,z_source,dip_angle,w,l,ss[0],ds[0])

    
    
    
    
    
    
    
    # supertitle = f"{self.name} Surface Displacement from {self.samples} samples {suptitle_extra}"
A =  np.array([1,4,25])
Zsrc = np.array([4,10])
Dip  = np.array([3])
DSlip = np.array([0,10])
SSlip = np.array([0,10])


A = A.astype(float)
#W = W.astype(float)

#L = L.astype(float)

Zsrc = Zsrc.astype(float)

Dip = Dip.astype(float)

DSlip = DSlip.astype(float)
SSlip = SSlip.astype(float)
count  = 0
for a in A:
    for z in Zsrc:
        for dip in Dip:
            for dslip in DSlip:
                
                for sslip in SSlip:
                    
                    count +=1
                    try:
                        calc_static(count,w=np.sqrt(a),l=np.sqrt(a),z_source=z,dip_angle=dip,Slip=np.array([sslip,dslip]))
                    except:
                        print('Invalid choice')
