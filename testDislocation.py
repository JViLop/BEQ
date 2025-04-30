import numpy as np
from layered_disloc import layered_disloc
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
def test1():

    km = 1000
    # station 
    xstn = np.arange(-100,100.0 + dr,dr)*km
    ystn = np.arange(-100.0,100.0 + dr,dr)*km

    Xstn,Ystn = np.meshgrid(xstn,ystn)
    xstn_flat,ystn_flat =  Xstn.flatten(),Ystn.flatten()

    zstn_flat = np.zeros(xstn_flat.shape)
    
    
    # source
    test_name = 'Dip Slip'
    xsrc = np.array([0.0])
    ysrc = np.array([0.0])
    zsrc = np.array([40.0])*km
    
    width = np.array([10.0])*km
    length = np.array([10.0])*km
    strike = np.array([90.0])
    dip = np.array([90.0])
    
    test_name = 'Dip Slip'
    
    ss = np.array([0.0])
    ds = np.array([50.0])
    ts = np.array([0.0])
   # definitions
    name = 'Gorkha'
    edks = 'poisson_half_space.edks'
    prefix = 'poisson_half_space'
    BIN_EDKS = '${HOME}/edks/bin'

    # reshape to form vectors.
    # source coordinates
    xS = xsrc
    yS = ysrc
    
    zS = zsrc
    aS = np.ones_like(dip)*(width[0]*length[0])
    strike = np.ones_like(dip)*90
    SSlip = ss
    DSlip = ds
    
    # calculate the GF's
    rake = 0 * np.ones(SSlip.shape)
    slip = 1 * np.ones(SSlip.shape)
    GFeSS, GFnSS, GFzSS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, xstn_flat, ystn_flat,\
                                  edks, prefix, BIN_EDKS)  
    
    rake = 90* np.ones(SSlip.shape)
    slip = 1 * np.ones(SSlip.shape)
    GFeDS, GFnDS, GFzDS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, xstn_flat, y_flat,\
                                  edks, prefix, BIN_EDKS)
    
    
     # compute forward displacement calculation.
    print(GFeDS.shape, GFeSS.shape)	
    dE = np.dot(GFeSS,SSlip) + np.dot(GFeDS, DSlip)
    dN = np.dot(GFnSS,SSlip) + np.dot(GFnDS, DSlip)
    dZ = np.dot(GFzSS,SSlip) + np.dot(GFzDS, DSlip)
    
    print(dE.shape)
    print(max(dE))
    print(max(dZ))
    print(max(dN))
    
    xstn ,ystn = xstn/km,ystn/km
    Xstn, Ystn = np.meshgrid(xstn,ystn) 
    shape_stn = (len(ystn),len(xstn))
    fig,ax = plt.subplots(dpi=900)
    
    # supertitle = f"{self.name} Surface Displacement from {self.samples} samples {suptitle_extra}"
    
    Ux = dE
    Uy = dN
    Uz = dZ
    

    fig, axes = plt.subplots(3,1,figsize=(5,10),dpi=600)
    dict = {'x':dE,'y':dN,'z':dZ} 
    df = pd.DataFrame(dict)    
    for i,parameter_id in enumerate(df.keys()):
       
       parameter = df[parameter_id].values
       supertitle = f"EDKS Synthetic Test {test_name}"
    
       # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
       im = axes[i].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
       
       #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
       axes[i].margins(0)
            
       fig.colorbar(im, ax=axes[i],shrink=0.8,label='{}'.format('Displacement (m)'))

    
      
       axes[i].set_ylabel('Trench-normal distance (km)',fontsize=11)
       axes[i].set_xlabel('Along-strike distance (km)',fontsize=11)
       axes[i].set_aspect('equal', 'box')
       axes[i].tick_params(labelsize=12)
       #axes[i].invert_yaxis()
    fig.suptitle(supertitle,x=0.5,y=0.98,fontsize=13,fontweight='bold') # uncomment later
    plt.tight_layout()
    
    fig.savefig(supertitle.replace(" ","_") + '.png')
    plt.close()
    '''
    dEE = dE[:,0].reshape((Ny, Nx))
    dNN = dN[:,0].reshape((Ny, Nx))
    dZZ = dZ[:,0].reshape((Ny, Nx))
   
   # plot the vector fields.
   import pylab as PL
   #PL.figure()
   norm = np.sqrt(dE * dE + dN * dN)
   #PL.quiver(x, y, dE/norm, dN/norm)
   #PL.axis('equal')
   #PL.figure()
   #PL.quiver(x, y, dZ * 0.0, dZ/np.abs(dZ))
   #PL.axis('equal')

   PL.figure()
   PL.subplot(1,3,1)
   PL.pcolor(XX, YY, dEE)
   PL.colorbar()
   PL.axis('equal')
   PL.subplot(1,3,2)
   PL.pcolor(XX, YY, dNN)
   PL.colorbar()
   PL.axis('equal')
   PL.subplot(1,3,3)
   PL.pcolor(XX, YY, dZZ)
   PL.axis('equal')
   PL.colorbar()


   PL.savefig('d.png')
   print(PL)
	'''
if __name__ == '__main__':
   test1()
