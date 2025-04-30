# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 14:36:32 2024

@author: joanv
"""

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os 
from matplotlib.colors import TwoSlopeNorm

  

### Illapel ###
data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Illapel/model/kinematic/step_052-001.h5'

f = h5py.File(data,'r')

Tr = np.array(f['ParameterSets']['risetime']).flatten()
Vr = np.array(f['ParameterSets']['rupturevelocity']).flatten()
U_dip = np.array(f['ParameterSets']['dipslip']).flatten()
U_strike =  np.array(f['ParameterSets']['strikeslip']).flatten()



fig, axs = plt.subplots(nrows=2, ncols=2,dpi=500)

axs[0, 0].hist(U_dip, bins=100)
axs[0,0].axvline(x = np.quantile(U_dip,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(U_dip,q=0.05),1)))

axs[0,0].axvline(x = np.quantile(U_dip,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(U_dip,q=0.95),1)))
axs[0,0].text(round(np.quantile(U_dip,q=0.95),1),1e5,'$Q_{95\%}$')
axs[0,0].text(round(np.quantile(U_dip,q=0.05),1),1e5,'$Q_{5\%}$')
axs[0, 0].set_title('$U_{dip}$ (m)')
axs[0,0].legend(fontsize = 7)
axs[0, 1].hist(U_strike, bins=100)
axs[0,1].axvline(x = np.quantile(U_strike,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(U_strike,q=0.05),1)))
axs[0,1].axvline(x = np.quantile(U_strike,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(U_strike,q=0.95),1)))
axs[0,1].text(round(np.quantile(U_strike,q=0.95),1),1e5,'$Q_{95\%}$')
axs[0,1].text(round(np.quantile(U_strike,q=0.05),1),1e5,'$Q_{5\%}$')
axs[0, 1].set_title('$U_{strike}$ (m)')
axs[0,1].legend(fontsize = 7)
axs[1, 0].hist(Tr, bins=100)
axs[1,0].axvline(x = np.quantile(Tr,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(Tr,q=0.05),1)))
axs[1,0].axvline(x = np.quantile(Tr,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(Tr,q=0.95),1)))
axs[1,0].text(round(np.quantile(Tr,q=0.95),1),1e5,'$Q_{95\%}$')
axs[1,0].text(round(np.quantile(Tr,q=0.05),1),1e5,'$Q_{5\%}$')
axs[1, 0].set_title('$T_{r}$ (s)')
axs[1,0].legend(fontsize = 7)
axs[1, 1].hist(Vr, bins=100)
axs[1,1].axvline(x = np.quantile(Vr,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(Vr,q=0.05),1)))
axs[1,1].axvline(x = np.quantile(Vr,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(Vr,q=0.95),1)))
axs[1,1].text(round(np.quantile(Vr,q=0.95),1),1e5,'$Q_{95\%}$')
axs[1,1].text(round(np.quantile(Vr,q=0.05),1),1e5,'$Q_{5\%}$')
axs[1, 1].set_title('$V_{r}$ (km/s)')
axs[1,1].legend(fontsize = 7)
fig.tight_layout()
fig.suptitle('Illapel',fontweight='bold')
plt.show()


### Iquique ###

data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Iquique/model/kinematic/step_056.h5'

f = h5py.File(data,'r')
Np= 132
nramp = 3 
d = np.array(f['Sample Set']).T
Tr = d[2*Np+nramp:3*Np+nramp,:]
Vr = d[3*Np+nramp:4*Np+nramp,:]
U_dip = d[:Np,:]
U_stk = d[Np:2*Np,:]
U = np.sqrt(U_dip**2 + U_stk**2)
slip_velocity = U/Tr


plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('Slip Velocity')
plt.title('Iquique')
plt.show()
plt.close()

plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('Tr')
plt.title('Iquique')
plt.show()
plt.close()
plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('U')
plt.title('Iquique')
plt.show()
plt.close()
plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
plt.xlabel('Tr')
plt.ylabel('U')
plt.title('Iquique')
plt.show()

# ### Gorkha ###
data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Gorkha/model/kinematic/step_012.h5'
f = h5py.File(data,'r')

Np= 162
nramp = 0 
d = np.array(f['Sample Set']).T
Tr = d[2*Np+nramp:3*Np+nramp,:]
Vr = d[3*Np+nramp:4*Np+nramp,:]
U_dip = d[:Np+nramp,:]
U_stk = d[Np+nramp:2*Np+nramp,:]
U = np.sqrt(U_dip**2 + U_stk**2)


slip_velocity = U/Tr

plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('Slip Velocity')
plt.title('Gorkha')
plt.show()
plt.close()

plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('Tr')
plt.title('Gorkha')
plt.show()
plt.close()
plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('U')
plt.title('Gorkha')
plt.show()
plt.close()
plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
plt.xlabel('Tr')
plt.ylabel('U')
plt.title('Gorkha')
plt.show()

'''
# ### Tohoku ###
dir_file = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Tohoku/model/catmip-theta44.bin'
data = np.fromfile(dir_file,'double').reshape((1000000,866))

plt.hist2d(Vr.flatten(),slip_velocity.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('Slip Velocity')
plt.title('Tohoku')
plt.show()
plt.close()

plt.hist2d(Vr.flatten(),Tr.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('Tr')
plt.title('Tohoku')
plt.show()

plt.close()
plt.hist2d(Vr.flatten(),U.flatten(),bins=100)
plt.xlabel('Vr')
plt.ylabel('U')
plt.title('Tohoku')
plt.show()
plt.close()

plt.hist2d(Tr.flatten(),U.flatten(),bins=100)
plt.xlabel('Tr')
plt.ylabel('U')
plt.title('Tohoku')
plt.show()

'''