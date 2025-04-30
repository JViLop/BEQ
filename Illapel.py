# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:04:09 2023

@author: joanv
"""

# from utils.geom_reader import geom_reader
# from utils.model_reader import model_reader
# before class 'EQ_model()'
# #geom_reader('Illapel')        
# model_reader('Illapel','kinematic','h5',(10,17),18)    


  
data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Illapel/model/kinematic/step_052-001.h5'

import h5py
import numpy as np
import matplotlib.pyplot as plt
f = h5py.File(data,'r')


Tr = np.array(f['ParameterSets']['risetime']).T
# Vr = np.array(f['ParameterSets']['rupturevelocity']).T
# U = np.sqrt(np.array(f['ParameterSets']['dipslip']).T**2 + np.array(f['ParameterSets']['strikeslip']).T**2)


# slip_velocity = U/Tr

# mean_slip_velocity = np.mean(slip_velocity,axis=1)
# mean_Vr  = np.mean(Vr,axis=1)
# #plt.scatter(mean_Vr ,mean_slip_velocity)
# # plt.scatter(slip_velocity,Vr)
# plt.xlabel('Slip Velocity')
# plt.ylabel('Vr')
# plt.title('Illapel')
# plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
# plt.scatter(Vr,slip_velocity)
# plt.scatter(np.array([[1,2],[3,4]]),np.array([[-1,-1],[3,3]]))
from utils.model_reader import EQ_model

#model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0,sampling=True, nsamples=50000)
# model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0,sampling=True, nsamples=10000)
model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0,sampling=True, nsamples=500)
#model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0)
#model.plot_corr()
# model.plot_corr_matrix()

# model.geometry_txt()

# model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
# model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
# model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
# model.plot_distribution('std','std_U_perp',(8,10),'$\sigma{(U_{\perp})}$',padding=10)
# model.plot_distribution('std','std_U_parallel',(8,10),'$\sigma{(U_{||})}$',padding=10)
# model.plot_distribution('skew','skew_U_perp',(8,10),'$skewness {(U_{\perp})}$',padding=10)
# model.plot_distribution('skew','skew_U_parallel',(8,10),'$skewness {(U_{||})}$',padding=10)

# model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
# model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)
# model.plot_distribution('std','std_Tr',(8,10),'$\sigma{(Tr)}$',padding=10)
# model.plot_distribution('std','std_Vr',(8,10),'$\sigma{(Vr)}$',padding=10)
# model.plot_distribution('skew','skew_Tr',(8,10),'$skewness {(Tr)}$',padding=10)
# model.plot_distribution('skew','skew_Vr',(8,10),'$skewness {(Vr)}$',padding=10)
