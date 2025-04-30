# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 15:27:47 2024

@author: joanv
"""

from utils.model_reader import EQ_model

model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,
           'depth':3,'strike':4,'dip':5},nrows_skip=1,nramp=0,rake=107,sampling=True, nsamples=1000)
model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0,sampling=True, nsamples=1000)
model = EQ_model('Iquique','kinematic','h5',(11,12),17,nramp=3,sampling=True, nsamples=1000)
model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99,sampling=True, nsamples=1000)
static_model = EQ_model('Tohoku','static','binary',(9,24),29,nrows_skip=2,model_file_shape=(432,1000000),sampling=True, nsamples=1000)
kinematic_model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=1000)

