# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 20:52:58 2024

@author: joanv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

name = 'Gorkha'
shape = (9,18)
patchsize = 10
def array_formatter(array,shape_geom):
        return np.flip(array.reshape(shape_geom,order='F'),axis=0).flatten()
Cd = np.loadtxt(f'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/cov_okada_observables/1000_samples/{name}/cov_txt/cov_matrix_{name}_displacement_nsamples_1000.txt')
#Cd = np.loadtxt('C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/cov_okada_observables/1000_samples/Tohoku/kinematic/cov_txt/cov_matrix_Tohoku_displacement_nsamples_1000.txt')
print(Cd.shape)
G = np.loadtxt(f'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/GF/static/displacement/{name}/kinematic/{name}_kinematic_Static_GF_displacement.txt')
print(G.shape)
d = np.loadtxt(f'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/{name}_100_kinematic_model_displacement.txt')
print(d.shape)
d = d.flatten(order='F')


inv_Cd = np.linalg.inv(Cd)

df = pd.read_csv(f'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')
df = df.drop(df.columns[0],axis=1)

Uparallel = array_formatter(df['U_parallel'].values,shape)
Uparallel  = Uparallel + np.random.normal(0, 1, size=Uparallel.shape)
Uperp = array_formatter(df['U_perp'].values,shape)
Uperp = array_formatter(df['U_perp'].values,shape) + np.random.normal(0, 1, size=Uperp.shape)

mprior = np.concatenate((Uparallel,Uperp))
Cm_prior = 1e-6*np.eye(mprior.shape[0])
inv_Cm_prior = np.linalg.inv(Cm_prior)


C = np.linalg.inv(np.matmul(G.T,np.matmul(inv_Cd,G)) + inv_Cm_prior)

D = np.matmul(G.T,np.matmul(inv_Cd,d)) + np.matmul(inv_Cm_prior,mprior)


m = np.matmul(C,D)




mean_slip = np.sqrt(m[:len(m)//2]**2+m[len(m)//2:]**2)

nrows,ncols = shape[0],shape[1]
x = np.arange(patchsize/2,ncols*patchsize,patchsize)
ydown = np.arange(-(nrows-1/2)*patchsize,0,patchsize)




fig, ax = plt.subplots()
param = mean_slip.reshape(nrows,ncols)
im = ax.pcolormesh(x,ydown,param,edgecolors='k',cmap='YlOrRd')

#ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
ax.margins(0)

        
fig.colorbar(im, ax=ax,shrink=0.38,label='slip(m)')
ax.set_title(f'{name} Inverted Static Slip Model from linear-Gaussian inversion',fontsize=10)
ax.set_ylabel('Down-dip distance (km)',fontsize=10)
ax.set_xlabel('Along-strike distance (km)',fontsize=10)
ax.set_aspect('equal', 'box')
ax.tick_params(labelsize=10)
fig.savefig(f'{name} Inverted Static Slip Model from linear-Gaussian inversion',dpi=600)



