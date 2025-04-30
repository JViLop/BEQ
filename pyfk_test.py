# -*- coding: utf-8 -*-
"""
Created on Tue May  7 22:59:00 2024

@author: joanv
"""
import pandas as pd
import numpy as np
import re
from pyfk import SourceModel, SeisModel, Config
from pyfk import calculate_sync
import os
from pyfk.gf.gf import calculate_gf
from pyfk import generate_source_time_function



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

df = pd.read_csv('Tohoku_mean_kinematic_model.csv')

nrows,ncols = 9,24
shape = (nrows,ncols)
patchsize = 29


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





def Mt(M0,strike,dip,rake):
    strike = strike*(np.pi/180)
    rake = rake*(np.pi/180)
    dip = dip*(np.pi/180)
    
    Mxx = -(np.sin(dip)*np.cos(rake)*np.sin(2*strike) + np.sin(2*dip)*np.sin(rake)*(np.sin(strike))**2)
    Mxy = np.sin(dip)*np.cos(rake)*np.cos(2*strike) + (1/2)*np.sin(2*dip)*np.sin(rake)*np.sin(2*strike)
    Mxz = -(np.cos(dip)*np.cos(rake)*np.cos(strike) + np.cos(2*dip)*np.sin(rake)*np.sin(strike))
    Myy = np.sin(dip)*np.cos(rake)*np.cos(2*strike) - np.sin(2*dip)*np.sin(rake)*(np.cos(strike))**2
    Myz = -(np.cos(dip)*np.cos(rake)*np.sin(strike) - np.cos(2*dip)*np.sin(rake)*np.cos(strike))
    Mzz = np.sin(2*dip)*np.sin(rake)
    
    return [Mxx,Mxy,Mxz,Myy,Myz,Mzz]


def build_model(path_model):
        path = os.path.join(path_model)
        f = open(path, "r")
        content = f.readlines()
        model = content[12:]
        keys = re.findall('\w+',content[11],flags=re.DOTALL)
        fmt = re.compile('\d+.\d+')
        layers = [list(map(float,re.findall(fmt,m))) for m in model]
        layers = np.array(layers)
        model_parameters = dict(zip(keys,layers))
        l = np.zeros((5,6))
        vs = layers[:,2]
        qs = layers[:,5]
        vp = layers[:,1]
        qp = layers[:,4]
        l[:,0] = layers[:,0]
        l[:,3] = layers[:,3]
        l[:,1] = vs
        l[:,2] = vp
        l[:,4] = qs
        l[:,5] = qp
        return l,model_parameters
    

def get_mu(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VS'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        return mu
    
def get_M0(mu,slip,patchsize):
        A = (patchsize*1e3)**2
        M0 = (A*slip*mu)*1e7  # dyne.cm 
        return M0
        
receiver = np.array([-50,400])

STRIKE = 90 
NPTS = 512
DT = 0.1
A = (patchsize*1e3)**2
DR = np.linalg.norm(receiver  - Xsrcs[:,:2],axis=1)
AZ = np.arctan2(receiver[1]-Xsrcs[:,1],receiver[0]-Xsrcs[:,0])*(180/np.pi)

for i,x in enumerate(Xsrcs):
    strike,dip,rake,depth = STRIKE,DIP[i],RAKE[i],DEPTH[i]
    
    U = SLIP[i] 
    
    dr,az = DR[i],AZ[i]
    model_array,model_dict = build_model('model_test')
    mu = get_mu(model_dict)
    M0 = (mu*A*U)*1e7
    
    
    M = Mt(strike,dip,rake,M0)
    Tr = DURATION[i]
    

    
    
    print(model_array)
    
    
    
    source = SourceModel(sdep=depth, srcType="dc",source_mechanism =M)
    print(source._source_mechanism)
    
    model = SeisModel(model=model_array)
    
    config = Config(
            model=model,
            source=source,
            npt=NPTS,
            dt=DT,
            receiver_distance=[dr])
    
    
    gf = calculate_gf(config)
    
    
    
    
    source_time_function=generate_source_time_function(dura=Tr, rise=0.5, delta=gf[0][0].stats.delta)
    sync_result = calculate_sync(gf, config, az, source_time_function)
    
    
    
    
    l = sync_result[i].integrate()
    l.plot(outfile= f'{i}_az_{dr}',automerge=False)
