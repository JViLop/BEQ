# -*- coding: utf-8 -*-
"""
Created on Thu Feb  1 09:41:29 2024

@author: joanv
"""

import numpy as np 
import os
import h5py
import matplotlib.pyplot as plt

# f = h5py.File('static.gf.h5','r')  
# print(f.keys())
# print(f['static.gf'])
# f1 = h5py.File('green.h5','r')
# print(f1.keys())
# print(f1['green'])

f = h5py.File('step_final.h5')
print(f['ParameterSets'].keys())
dipslip =  f['ParameterSets']['dipslip']
strikeslip =f['ParameterSets']['strikeslip']

mean_ds = np.mean(dipslip,axis=0)
mean_ss = np.mean(strikeslip,axis=0)

mean_slip = np.sqrt(mean_ds**2+mean_ss**2)
patchsize =10
nrows,ncols = 9,18
x = np.arange(patchsize/2,ncols*patchsize,patchsize)
y = np.arange(patchsize/2,nrows*patchsize,patchsize)
ydown = np.arange(-(nrows-1/2)*patchsize,0,patchsize)



X, Y = np.meshgrid(x, y)



fig, ax = plt.subplots()
param = mean_slip.reshape(nrows,ncols)
im = ax.pcolormesh(x,ydown,param,edgecolors='k',cmap='YlOrRd')

#ax.plot(X.flat, Y.flat, '.', color='k',markersize=0.5)
ax.margins(0)

        
fig.colorbar(im, ax=ax,shrink=0.38,label='slip(m)')
ax.set_title('Gorkha Inverted Static Slip Model with $Cd$ diagonal',fontsize=10)
ax.set_ylabel('Down-dip distance (km)',fontsize=10)
ax.set_xlabel('Along-strike distance (km)',fontsize=10)
ax.set_aspect('equal', 'box')
ax.tick_params(labelsize=10)
fig.savefig('Gorkha Inverted Static Slip Model with $Cd$ diagonal.png',dpi=600)

# pwd = os.getcwd()
# name ='data'
# suffix = '.txt'

# file_dir = os.path.join(pwd,'files')
# data_dir = os.path.join(file_dir,'data')
# filein = os.path.join(data_dir,name+suffix)

# h5file = h5py.File(name=name+'.h5', mode='w')
# data = np.loadtxt(filein, dtype='float32')
# h5file.create_dataset(name=name, data=data)
# h5file.close()
# name ='green'
# filein = os.path.join(data_dir,name+suffix)

# h5file = h5py.File(name=name+'.h5', mode='w')
# data = np.loadtxt(filein, dtype='float32')
# h5file.create_dataset(name=name, data=data)
# h5file.close()
# name ='cd'
# filein = os.path.join(data_dir,name+suffix)

# h5file = h5py.File(name=name+'.h5', mode='w')
# data = np.loadtxt(filein, dtype='float32')
# h5file.create_dataset(name=name, data=data)
# h5file.close()