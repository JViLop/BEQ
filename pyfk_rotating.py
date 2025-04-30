# -*- coding: utf-8 -*-
"""
Created on Wed May  8 19:25:19 2024

@author: joanv
"""

import numpy as np
import obspy as obs
import matplotlib.pyplot as plt
import pandas as pd

def array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()

def proj_ysrc_coords(dip,nrows):
    dip = dip[:nrows]
    proj_dysrc = -patchsize*np.cos(dip*np.pi/180) # in meters
    proj_ysrc = np.zeros_like(proj_dysrc)
    for i in range(len(proj_ysrc)):
        proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
    ysrc = np.flip(proj_ysrc)    
    return ysrc

#df = pd.read_csv('Tohoku_mean_kinematic_model.csv')

#nrows,ncols = 9,24
#shape = (nrows,ncols)
#patchsize = 29

## Illapel 
df = pd.read_csv('Illapel_mean_kinematic_model.csv')

nrows,ncols = 10,17
shape = (nrows,ncols)
patchsize = 18



Uperp = array_formatter(df['U_perp'].values,shape)
Uparallel =  array_formatter(df['U_parallel'].values,shape)
SLIP = array_formatter(df['Slip'].values,shape)
RAKE = array_formatter(np.arctan2(Uparallel,Uperp)*(180/np.pi),shape)
DIP = array_formatter(df['dip'].values,shape)
DEPTH = array_formatter(df['depth'].values,shape)
DURATION = array_formatter(df['Tr'].values,shape)


xsrc = np.arange((1/2)*patchsize,ncols*patchsize, patchsize)
ysrc = proj_ysrc_coords(DIP,nrows)
Xsrc,Ysrc= np.meshgrid(xsrc,ysrc)
        
xsrc_flat = Xsrc.flatten()
ysrc_flat = Ysrc.flatten()
zsrc_flat = DEPTH  # it mush be in km (not meters)
Xsrcs  = np.array([ysrc_flat,xsrc_flat,zsrc_flat]).T

n = len(Xsrcs)
# Tohoku
#x,y = -50,400

# Illapel
x,y = -40,160

receiver = np.array([x,y])

STRIKE = 90 
NPTS = 300
DT = 0.1
A = (patchsize*1e3)**2
DR = np.linalg.norm(receiver  - Xsrcs[:,:2],axis=1)
AZ = np.arctan2(receiver[1]-Xsrcs[:,1],receiver[0]-Xsrcs[:,0])

# for comp in 'ZRT':
#     strm = obs.read(f'traces_{comp}.mseed')
#     traces = 0
#     for i in range(216):
#         traces += strm[i].data/100
        
#     plt.plot(traces,label=comp)
#     plt.legend()


# last to be used see test 21 ####
#
# folder = f'{x}_{y}/'
# path = './pyfk_test/test21/' + folder
# strmR = obs.read(path + 'traces_R.mseed')
# strmT = obs.read(path +'traces_T.mseed')
# strmR.sort()
# strmT.sort()
# strmN = obs.Stream()  
# strmE = obs.Stream() 

# strm_dict = {'Z':obs.read(path + 'traces_Z.mseed')}   
# dt = strm_dict['Z'][0].stats.delta
# npts =strm_dict['Z'][0].stats.npts
# t = np.arange(0,npts*dt,dt)
# for i in range(n):
    
#     caz = np.cos(AZ[i])
#     saz = np.sin(AZ[i])
#     N = (caz*strmR[i].data - saz*strmT[i].data)
#     E = (saz*strmR[i].data + caz*strmT[i].data)
#     traceN,traceE = obs.Trace(N),obs.Trace(E)
#     strmN +=traceN
#     strmE +=traceE
    
# strm_dict['N'] = strmN
# strm_dict['E'] = strmE
# keys = {'N':'y','E':'x','Z':'z'}
# fig,ax = plt.subplots(dpi=900)
# for key in list(strm_dict.keys()):
#     strm  = strm_dict[key]
#     traces = 0
#     for i in range(n):
#         traces += strm[i].data*1e-2

#     ax.plot(t,traces,lw = 1,label=keys[key])

#     ax.axhline(y = 0,color='k', ls='dashed',lw=0.75)
#     ax.set_title(f'Receiver x = {y}km, y ={x}km')
#     ax.set_ylabel('acceleration (m/s)')
#     ax.set_xlabel('Time (s)')
# ax.legend(ncols=1)

####



#From this line downward, creates the animation
comp = 'N'
# xstn = np.flip(np.array([50,0,-50,-100,-150,-200,-250,-300]))
# ystn = np.array([0,50,100,150,200,250,300,350,400,450,500,550,600])

xstn = np.flip(np.array([50,0,-50,-100,-150,-200,-250]))
ystn = np.array([0,50,100,150,200,250,300,350])
Xstn,Ystn = np.meshgrid(xstn,ystn)
xflat,yflat = Xstn.flatten(order='F'),Ystn.flatten(order='F')
receivers = np.column_stack((xflat,yflat))
ncols,nrows = len(ystn),len(xstn)
G = np.zeros((nrows*ncols,405))
for k,r in enumerate(receivers):
    
    x,y = r[0],r[1]
    receiver = np.array([x,y])
    
    STRIKE = 90 
    NPTS = 512
    DT = 0.1
    A = (patchsize*1e3)**2
    DR = np.linalg.norm(receiver  - Xsrcs[:,:2],axis=1)
    AZ = np.arctan2(receiver[1]-Xsrcs[:,1],receiver[0]-Xsrcs[:,0])
    
    # for comp in 'ZRT':
    #     strm = obs.read(f'traces_{comp}.mseed')
    #     traces = 0
    #     for i in range(216):
    #         traces += strm[i].data/100
            
    #     plt.plot(traces,label=comp)
    #     plt.legend()
    
    folder = f'{x}_{y}/'
    path = './pyfk_test/test22/' + folder
    strmR = obs.read(path + 'traces_R.mseed')
    strmT = obs.read(path +'traces_T.mseed')
    strmR.sort()
    strmT.sort()
    strmN = obs.Stream()  
    strmE = obs.Stream() 
    
    strm_dict = {'Z':obs.read(path + 'traces_Z.mseed')}  
    dt = strm_dict['Z'][0].stats.delta
    npts =strm_dict['Z'][0].stats.npts
    t = np.arange(0,npts*dt,dt)
    for i in range(n):
        
        caz = np.cos(AZ[i])
        saz = np.sin(AZ[i])
        N = (caz*strmR[i].data - saz*strmT[i].data)
        E = (saz*strmR[i].data + caz*strmT[i].data)
        traceN,traceE = obs.Trace(N),obs.Trace(E)
        strmN +=traceN
        strmE +=traceE
        
    strm_dict['N'] = strmN
    strm_dict['E'] = strmE
    keys = {'N':'y','E':'x','Z':'z'}
    # for key in list(strm_dict.keys()):
    for key in [comp]:
        strm  = strm_dict[key]
        traces = 0
        for i in range(n):
            traces += strm[i].data*1e-2
        
        d = traces
        G[k,:] = d
        
    
    



G = G.reshape(nrows,ncols,405)
import numpy as np
from matplotlib import pyplot as plt, animation
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True


G0 = G[:,:,::50]
nframes = G0.shape[2]
fig, ax = plt.subplots()
# xstn = np.flip(np.array([50,0,-50,-100,-150,-200,-250,-300]))
# ystn = np.array([0,50,100,150,200,250,300,350,400,450,500,550,600])


xstn = np.flip(np.array([50,0,-50,-100,-150,-200,-250]))
ystn = np.array([0,50,100,150,200,250,300,350])

cax = ax.pcolormesh(ystn, xstn, G0[:, :, 0], vmin = min(np.min(G0),-np.max(G0)),vmax =  max(-np.min(G0),np.max(G0)),cmap='bwr')
ax.set_title(f'{comp} displacement')
fig.colorbar(cax)

def animate(i):
    cax.set_array(G0[:, :, i].flatten())

anim = animation.FuncAnimation(fig, animate,frames =nframes,interval = 100)
anim.save(filename=f"{comp}_Illapel.mp4", writer="ffmpeg")
plt.show()
    
 