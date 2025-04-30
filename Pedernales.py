# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:04:44 2023

@author: joanv
"""

# from utils.geom_reader import geom_reader
# from utils.model_reader import model_reader

       
# model_reader('Pedernales','kinematic','h5',(8,10),15)   

  
# data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Pedernales/model/kinematic/step_053.h5'

# import h5py
# import numpy as np
# f = h5py.File(data,'r')
# data = f['Sample Set']




from utils.model_reader import EQ_model

# model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99,sampling=True, nsamples=100)
model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99,sampling = True, nsamples=100)

#model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99)
#model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99,sampling=True, nsamples=50000)
# model.plot_corr()
# model.plot_corr_matrix()

#model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)


# model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10),
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
