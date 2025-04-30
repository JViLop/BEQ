# -*- coding: utf-8 -*-
"""
Created on Tue May 28 16:26:45 2024

@author: joanv
"""

import matplotlib.pyplot as plt
import numpy as np

stn = '0920'
f = open(f'mudpy_set/{stn}.KIN','r')
D = f.readlines()


time = np.array([float(d.split()[2]) for d in D[6:-1] ])
n = len(time)

N = np.array([float(d.split()[3]) for d in D[6:-1] ])

E = np.array([float(d.split()[4]) for d in D[6:-1] ])
Z = np.array([float(d.split()[5]) for d in D[6:-1] ])

f1 = np.argmax(E) - 125
f2 = np.argmax(E) + 100
time = np.arange(0,f2-f1)

N = N[f1:f2]
E = E[f1:f2]
Z = Z[f1:f2]
plt.plot(time,N,label='N')
plt.plot(time,E,label='E')
plt.plot(time,Z,label='Z')
plt.axhline(y=0,ls='dashed',color='black',lw=1)
plt.xlabel('Time (s)')

plt.xlabel('Time (s)')
plt.ylabel('Displacement (m)')
plt.title(f'{stn} GEONET')
plt.legend()
plt.savefig(f'mudpy_set/{stn}_geonet.png',dpi=600)