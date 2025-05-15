# -*- coding: utf-8 -*-
"""
Created on Sat Apr 13 23:00:19 2024

@author: joanv
"""
import numpy as np


def reader(model_file_dir):
    data = np.fromfile(model_file_dir,'double')
    return data

name = ['Iquique','Pedernales','Gorkha','Illapel']
path_file = lambda x : f'INPUT/{x}/model/kinematic/all_samples/bin_data/{x}_kinematic_n_all.dat'

dat = reader(path_file(name[3]))


name = ['Tohuku']
path_file = lambda x : f'INPUT/Tohoku/model/kinematic/all_samples/bin_data/Tohoku_kinematic_n_all.bin'

dat = reader(path_file('static'))
