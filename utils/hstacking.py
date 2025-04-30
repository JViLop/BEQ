# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 20:18:25 2023

@author: joanv
"""

import numpy as np


def stacking(d,keys,n):
    if n!=0:
        x = np.hstack((stacking(d,keys,n-1),d[keys[n+1]]))
    else:
        x = np.hstack((d[keys[n]],d[keys[n+1]]))
    return x
    

def hstacking(d):   
    keys = ['dipslip','strikeslip','rupturevelocity','risetime']
    n = len(keys)
    return stacking(d,keys,n-2)