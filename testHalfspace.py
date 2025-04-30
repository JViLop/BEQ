import numpy as np
from layered_disloc import layered_disloc
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import pandas as pd
import sys
import os

file_name = os.path.basename(sys.argv[0])
test_name = filename[4:-2]
names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]

def model_dict(names,geoms,patches,arrow_sizes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
    return model

def set_stn(c,control=0):
      
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn
  
def proj_ysrc_coords(patch,dip):
      proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
      proj_ysrc = np.zeros_like(proj_dysrc)
      for i in range(len(proj_ysrc)):
          proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
      ysrc = np.flip(proj_ysrc)
        
      return ysrc
   


def run():
    
    models = model_dict(names,geoms,patches,arrow_sizes)
    name = sys.argv[0]
    
       # units
    km = 1000.0 # in meters
    # definitions
    edks = 'poisson_half_space.edks'
    prefix = 'poisson_half_space'
    BIN_EDKS = '${HOME}/edks/bin'
    input_dir  = f'{name}_mean_kinematic_model.csv'
    
    nrows, ncols = models[name]['geom'][0], models[name]['geom'][1]
    patch = models[name]['patch']*km	
    factor = int(sys.argv[1])
    control = int(sys.argv[2])
    c = int(sys.argv[3])
    
    
    xstn,ystn = set_stn(c,control=1)
    
       
    XX, YY  = np.meshgrid(xstn,ystn)
    x = XX.flatten()
    y = YY.flatten()
    df = pd.read_csv(input_dir)
    df= df.drop(df.columns[0],axis=1)
    
    	
    dip = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
    strike =np.flip(df['strike'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
    depth = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
    Uperp = np.flip(df['U_perp'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
      
    
    
    
    Uparallel = np.flip(df['U_parallel'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
    Slip =  np.flip(df['Slip'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
    Slip = np.sqrt(Uperp**2 + Uparallel**2)
    x_hyp = df['Hypo_as'].values[0]
    y_hyp = -df['Hypo_dd'].values[0]
    xS = np.arange(patchsize/2 , ncols*patch,patch)
    yS = np.arange(-(nrows-1/2)*patch,0,patch)
       # shift accordingly at surface
    yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
    XS, YS = np.meshgrid(xS,yS)
    
       # reshape to form vectors.
       # source coordinates
    xS = XS.flatten()
    yS = YS.flatten()
    
    zS = depth*km
    aS = np.ones_like(dip)*(patch*patch)
    strike = np.ones_like(dip)*90
    SSlip = Uperp
    DSlip = Uparallel
    
    # calculate the GF's
    rake = 0 * np.ones(SSlip.shape)
    slip = 1 * np.ones(SSlip.shape)
    GFeSS, GFnSS, GFzSS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, x, y,\
                                         edks, prefix, BIN_EDKS)  
    
    rake = 90* np.ones(SSlip.shape)
    slip = 1 * np.ones(SSlip.shape)
    GFeDS, GFnDS, GFzDS = layered_disloc(xS, yS, zS, strike, dip, rake, slip, aS, x, y,\
                                 edks, prefix, BIN_EDKS)
    
    
    # compute forward displacement calculation.
    print(GFeDS.shape, GFeSS.shape)	
    dE = np.dot(GFeSS,SSlip) + np.dot(GFeDS, DSlip)
    dN = np.dot(GFnSS,SSlip) + np.dot(GFnDS, DSlip)
    dZ = np.dot(GFzSS,SSlip) + np.dot(GFzDS, DSlip)
    
    print(dE.shape)
    print(max(dE))
    print(max(dN))
    print(max(dZ))
    
    xstn ,ystn = xstn/km,ystn/km
    Xstn, Ystn = np.meshgrid(xstn,ystn) 
    shape_stn = (len(ystn),len(xstn))
    fig,ax = plt.subplots(dpi=900)
    name  = 'Iquique'
    # supertitle = f"{self.name} Surface Displacement from {self.samples} samples {suptitle_extra}"
    supertitle = f"{name} Surface Displacement with EDKS ({test_name}) "
    Ux = dE
    Uy = dN
    Uz = dZ
    
    im = ax.pcolormesh(xstn,ystn,Uz.reshape(shape_stn),edgecolors='k',linewidths=0.1,cmap='bwr',norm = TwoSlopeNorm(0,vmin = min(Uz.min(),-Uz.max()),vmax=max(-Uz.min(),Uz.max()))) 
    fig.colorbar(im,ax = ax,shrink = 0.4,label= 'Uplift (m)')
    try:
    	ax.plot(x_hyp,y_hyp,marker = '*',color='yellow',markersize=14,markeredgecolor='black')
    except:
    	pass
    larrow = 10
    scale = 0.3
    q = ax.quiver(Xstn,Ystn,Ux,Uy,scale=scale,scale_units ='x', units='width',width=0.002,headwidth=4.5,headlength=6)
    ax.quiverkey(q, X=0.04, Y=1.04, U=larrow,label=f'{larrow}m', labelpos='N',fontproperties={'size':7})
    ax.set_ylabel('Trench-normal distance (km)',fontsize=9)
    ax.set_xlabel('Along-strike distance (km)',fontsize=9)
    ax.set_aspect('equal', 'box')
    ax.tick_params(labelsize=9)
    ax.set_title(supertitle,x=0.5,y=1.015,fontsize=10, fontweight='bold') # size before 9
    plt.tight_layout()
    
    
    fig.savefig(f'{name}_{test_name}_all_in_1')    
    plt.close()
    
    fig, axes = plt.subplots(3,1,figsize=(8,10),dpi=600)
    dict = {'x':dE,'y':dN,'z':dZ} 
    df = pd.DataFrame(dict)    
    for i,parameter_id in enumerate(df.keys()):
        parameter = df[parameter_id].values
        title = f"{parameter_id} Surface displacement"
    
    # parameter Displacements-related is already in correct order ie. row-wise as obtained from Okada
       	im = axes[i].pcolormesh(xstn,ystn,parameter.reshape(shape_stn),edgecolors='k', cmap='bwr',norm=TwoSlopeNorm(0,vmin=min(parameter.min(),-parameter.max()),vmax=max(-parameter.min(),parameter.max())))
    
    #ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
        axes[i].margins(0)
        fig.colorbar(im, ax=axes[i],shrink=0.8,label='{}'.format('Displacement (m)'))
        try:
            axes[i].plot(x_hyp,y_hyp,marker = '*',color='yellow',markersize=14,markeredgecolor='black')
        except:
            pass
        larrow = models[name]['arrow_size']
        scale = 0.3
        axes[i].set_ylabel('Trench-normal distance (km)',fontsize=11)
        axes[i].set_xlabel('Along-strike distance (km)',fontsize=11)
        axes[i].set_aspect('equal', 'box')
        axes[i].set_title(title,fontweight='bold',fontsize = 11)
        axes[i].tick_params(labelsize=12)
    #axes[i].invert_yaxis()
    fig.suptitle(supertitle,x=0.5,y=0.98,fontsize=13,fontweight='bold') # uncomment later
    plt.tight_layout()
    
    fig.savefig(f'{name}_{test_name}_all')
    plt.close()
    

if __name__ == '__main__':
   run()
