# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:04:09 2023

@author: joanv
"""

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
