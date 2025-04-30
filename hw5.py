# -*- coding: utf-8 -*-
"""
Created on Tue Oct 15 21:38:32 2024

@author: joanv
"""

import numpy as np

import matplotlib.pyplot as plt
sigma3 = np.array([0.1,0.30,0.50,1.0]) 
sigma1 = np.array([0.54,1.63,2.72,5.45])

angle = 60*np.pi/180
sigma_n = np.abs(0.5*((sigma3+sigma1) + (sigma1 - sigma3)*np.cos(2*angle)))
sigma_s = np.abs(0.5*((sigma1 - sigma3)*np.sin(2*angle)))

psi = np.arctan2(sigma_s,sigma_n)
plt.plot(sigma_n, sigma_s,marker='o')
plt.xlabel('$\sigma_n$')
plt.ylabel(r'$\tau$')
plt.xlim(0,2.5)
plt.ylim(0,2.5)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.draw()
