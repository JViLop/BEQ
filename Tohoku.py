# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 14:03:47 2023

@author: joanv
"""

# from utils.geom_reader import geom_reader
# from utils.model_reader import model_reader

# #geom_reader('Tohoku',nrows_skip=2)        
# #model_reader('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000)) 
# model_reader('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000))   

from utils.model_reader import EQ_model
import matplotlib.pyplot as plt
import time

# start = time.time()

# static_model = EQ_model('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000))
# static_model = EQ_model('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000),sampling=True, nsamples=100)
#static_model = EQ_model('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000),sampling=True, nsamples=10000)
# kinematic_model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=100)
#static_model.geometry_checker()
# static_model = EQ_model('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000))

#static_model = EQ_model('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000))
#static_model.plot_corr()
#plt.close()
#static_model.plot_corr_matrix()
#kinematic_model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000))
kinematic_model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=500)

# kinematic_model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=100)
#kinematic_model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=50000)

#kinematic_model.plot_corr()
#kinematic_model.plot_corr_matrix()
# end = time.time()
# print(end-start)

#Static model 
#static_model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
#static_model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
# static_model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
# static_model.plot_distribution('std','std_U_perp',(8,10),'$\sigma{(U_{\perp})}$',padding=10)
# static_model.plot_distribution('std','std_U_parallel',(8,10),'$\sigma{(U_{||})}$',padding=10)
# static_model.plot_distribution('skew','skew_U_perp',(8,10),'$skewness {(U_{\perp})}$',padding=10)
# static_model.plot_distribution('skew','skew_U_perp',(8,10),'$skewness {(U_{||})}$',padding=10)


# # Kinematic model

#kinematic_model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
# kinematic_model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
# kinematic_model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
# kinematic_model.plot_distribution('std','std_U_perp',(8,10),'$\sigma{(U_{\perp})}$',padding=10)
# kinematic_model.plot_distribution('std','std_U_parallel',(8,10),'$\sigma{(U_{||})}$',padding=10)
# kinematic_model.plot_distribution('skew','skew_U_perp',(8,10),'$skewness {(U_{\perp})}$',padding=10)
# kinematic_model.plot_distribution('skew','skew_U_parallel',(8,10),'$skewness {(U_{||})}$',padding=10)

# kinematic_model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
# kinematic_model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)
# kinematic_model.plot_distribution('std','std_Tr',(8,10),'$\sigma{(Tr)}$',padding=10)
# kinematic_model.plot_distribution('std','std_Vr',(8,10),'$\sigma{(Vr)}$',padding=10)
# kinematic_model.plot_distribution('skew','skew_Tr',(8,10),'$skewness {(Tr)}$',padding=10)
# kinematic_model.plot_distribution('skew','skew_Vr',(8,10),'$skewness {(Vr)}$',padding=10)
