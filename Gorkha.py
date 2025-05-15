# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:04:35 2023

@author: joanv
"""

# from utils.geom_reader import geom_reader
# from utils.model_reader import model_reader

# #geom_reader('Gorkha',header = [2,3,4,5,6],nrows_skip=1)     
# model_reader('Gorkha', 'kinematic', 'h5', (9,18), 10,header=[2,3,4,5,6],nrows_skip=1)

data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Gorkha/model/kinematic/step_012.h5'

import h5py
import numpy as np
import matplotlib.pyplot as plt
f = h5py.File(data,'r')

Np= 162
nramp = 0 
d = np.array(f['Sample Set']).T
# Tr = d[2*Np+nramp:3*Np+nramp,:]
# Vr = d[3*Np+nramp:4*Np+nramp,:]
# U_dip = d[:Np+nramp,:]
# U_stk = d[Np+nramp:2*Np+nramp,:]
# U = np.sqrt(U_dip**2 + U_stk**2)


# slip_velocity = U/Tr

# mean_slip_velocity = np.mean(slip_velocity,axis=1)
# mean_Vr  = np.mean(Vr,axis=1)
# #plt.scatter(mean_Vr ,mean_slip_velocity)
# plt.scatter(slip_velocity,Vr)
# plt.xlabel('Slip Velocity')
# plt.ylabel('Vr')
# plt.title('Gorkha')
from utils.model_reader import EQ_model

#model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,'depth':3,'strike':4,'dip':5},nrows_skip=1,nramp=0,rake=107,sampling=True, nsamples=50000)
# model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,'depth':3,'strike':4,'dip':5},nrows_skip=1,nramp=0,rake=107,sampling=True, nsamples=10000)
#model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,'depth':3,'strike':4,'dip':5},nrows_skip=1,nramp=0,rake=107)
model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,'depth':3,'strike':4,'dip':5},nrows_skip=1,nramp=0,rake=107,sampling=True, nsamples = 500)

# model.plot_corr()
# model.plot_corr_matrix()
# model.geometry_txt()

# model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
# model.plot_distribution('mean','U_parallel',(8,10),'$U_{\perp}(m)$',padding=10)
# model.plot_distribution('mean','U_perp',(8,10),'$U_{||}(m)$',padding=10)
# model.plot_distribution('std','std_U_parallel',(8,10),'$\sigma{(U_{\perp})}$',padding=10)
# model.plot_distribution('std','std_U_perp',(8,10),'$\sigma{(U_{||})}$',padding=10)
# model.plot_distribution('skew','skew_U_parallel',(8,10),'$skewness {(U_{\perp})}$',padding=10)
# model.plot_distribution('skew','skew_U_perp',(8,10),'$skewness {(U_{||})}$',padding=10)

# model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
# model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)
# model.plot_distribution('std','std_Tr',(8,10),'$\sigma{(Tr)}$',padding=10)
# model.plot_distribution('std','std_Vr',(8,10),'$\sigma{(Vr)}$',padding=10)
# model.plot_distribution('skew','skew_Tr',(8,10),'$skewness {(Tr)}$',padding=10)
# model.plot_distribution('skew','skew_Vr',(8,10),'$skewness {(Vr)}$',padding=10)
