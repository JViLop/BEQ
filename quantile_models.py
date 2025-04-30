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

  

def plot_model_hist(name, U_dip,U_strike, Tr,Vr,Hypo_dip,Hypo_strike):
    fig, axs = plt.subplots(nrows=3, ncols=2)
    
    axs[0, 0].hist(U_dip, bins=100)
    axs[0,0].axvline(x = np.quantile(U_dip,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(U_dip,q=0.05),1)))
    
    axs[0,0].axvline(x = np.quantile(U_dip,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(U_dip,q=0.95),1)))
    axs[0,0].text(round(np.quantile(U_dip,q=0.95),1),1e5,'$Q_{95\%}$')
    axs[0,0].text(round(np.quantile(U_dip,q=0.05),1),1e5,'$Q_{5\%}$')
    axs[0, 0].set_title('Along-dip/Rake-parallel $U_{||}$',fontsize=8)
    axs[0, 0].set_xlabel('$U_{||}$ (m)',fontsize=6)
    axs[0, 0].set_ylabel('count',fontsize=6)
    axs[0,0].legend(fontsize = 7,labelspacing=0.1)
    axs[0,0].tick_params(labelsize=7)
    axs[0, 1].hist(U_strike, bins=100)
    axs[0,1].axvline(x = np.quantile(U_strike,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(U_strike,q=0.05),1)))
    axs[0,1].axvline(x = np.quantile(U_strike,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(U_strike,q=0.95),1)))
    axs[0,1].text(round(np.quantile(U_strike,q=0.95),1),1e5,'$Q_{95\%}$')
    axs[0,1].text(round(np.quantile(U_strike,q=0.05),1),1e5,'$Q_{5\%}$')
    axs[0, 1].set_title('Along-strike/Rake-perpendicular $U_{\perp}$',fontsize=8)
    axs[0, 1].set_xlabel('$U_{\perp}$ (m)',fontsize=6)
    axs[0, 1].set_ylabel('count',fontsize=6)
    axs[0,1].tick_params(labelsize=6)
    axs[0,1].legend(fontsize = 7,labelspacing=0.1)
    axs[1, 0].hist(Tr, bins=100)
    axs[1,0].axvline(x = np.quantile(Tr,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(Tr,q=0.05),1)))
    axs[1,0].axvline(x = np.quantile(Tr,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(Tr,q=0.95),1)))
    axs[1,0].text(round(np.quantile(Tr,q=0.95),1),1e5,'$Q_{95\%}$')
    axs[1,0].text(round(np.quantile(Tr,q=0.05),1),1e5,'$Q_{5\%}$')
    axs[1, 0].set_title('Rupture Duration',fontsize=8)
    axs[1, 0].set_xlabel('$T_{r}$ (s)',fontsize=6)
    axs[1, 0].set_ylabel('count',fontsize=6)
    axs[1,0].tick_params(labelsize=6)
    axs[1,0].legend(fontsize = 7,labelspacing=0.1)
    axs[1, 1].hist(Vr, bins=100)
    axs[1,1].axvline(x = np.quantile(Vr,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(Vr,q=0.05),1)))
    axs[1,1].axvline(x = np.quantile(Vr,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(Vr,q=0.95),1)))
    axs[1,1].text(round(np.quantile(Vr,q=0.95),1),1e5,'$Q_{95\%}$')
    axs[1,1].text(round(np.quantile(Vr,q=0.05),1),1e5,'$Q_{5\%}$')
    axs[1, 1].set_title('Rupture Speed',fontsize=8)
    axs[1, 1].set_xlabel('$V_{r}$ (km/s)',fontsize=6)
    axs[1, 1].set_ylabel('count',fontsize=6)
    axs[1,1].tick_params(labelsize=6)
    axs[1,1].legend(fontsize = 7,labelspacing=0.1)
    
    axs[2, 0].hist(Hypo_dip, bins=100)
    axs[2,0].axvline(x = np.quantile(Hypo_dip,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(Hypo_dip,q=0.05),1)))
    axs[2,0].axvline(x = np.quantile(Hypo_dip,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(Hypo_dip,q=0.95),1)))
    axs[2,0].text(round(np.quantile(Hypo_dip,q=0.95),1),1e3,'$Q_{95\%}$')
    axs[2,0].text(round(np.quantile(Hypo_dip,q=0.05),1),1e3,'$Q_{5\%}$')
    axs[2,0].set_title('Along-dip/Rake-perpendicular Hypocenter ',fontsize=8)
    axs[2,0].set_xlabel('$x_{\perp}$ (km)',fontsize=6)
    axs[2,0].set_ylabel('count',fontsize=6)
    axs[2,0].tick_params(labelsize=6)
    axs[2,0].legend(fontsize = 7,labelspacing=0.1)
    
    axs[2, 1].hist(Hypo_strike, bins=100)
    axs[2, 1].axvline(x = np.quantile(Hypo_strike,q=0.05),color='black',linestyle='--',label='$Q_{5\%}$'+'$ = %s $'%(round(np.quantile(Hypo_strike,q=0.05),1)))
    axs[2, 1].axvline(x = np.quantile(Hypo_strike,q=0.95),color='red',linestyle='--',label='$Q_{95\%}$'+'$ = %s $'%(round(np.quantile(Hypo_strike,q=0.95),1)))
    axs[2, 1].text(round(np.quantile(Hypo_strike,q=0.95),1),1e3,'$Q_{95\%}$')
    axs[2, 1].text(round(np.quantile(Hypo_strike,q=0.05),1),1e3,'$Q_{5\%}$')
    axs[2, 1].set_title('Along-strike/Rake-parallel Hypocenter' ,fontsize=8)
    axs[2, 1].set_xlabel('$x_{||}$ (km)',fontsize=6)
    axs[2, 1].set_ylabel('count',fontsize=6)
    axs[2,1].tick_params(labelsize=6)
    axs[2,1].legend(fontsize = 7,labelspacing=0.1)
    plt.subplots_adjust(wspace=0.4,hspace = 0.9)
    fig.suptitle(f'{name}',fontweight='bold',y = 1.0)
    fig_dir = os.path.join(os.getcwd(),'posterior_prob')
    os.makedirs(fig_dir,exist_ok=True)
    file_dir = os.path.join(fig_dir,f'{name}.png')
    fig.savefig(file_dir,dpi=500)
    plt.show()
    

### Illapel ###
data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Illapel/model/kinematic/step_052-001.h5'

f = h5py.File(data,'r')

Tr = np.array(f['ParameterSets']['risetime']).flatten()
Vr = np.array(f['ParameterSets']['rupturevelocity']).flatten()
U_dip = np.array(f['ParameterSets']['dipslip']).flatten()
U_strike =  np.array(f['ParameterSets']['strikeslip']).flatten()


Hypo_as=np.array(f['ParameterSets']['hypo_as']).flatten()
Hypo_dd=np.array(f['ParameterSets']['hypo_dd']).flatten()

plot_model_hist('Illapel', U_dip,U_strike, Tr,Vr,Hypo_dd,Hypo_as)


### Iquique ###

data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Iquique/model/kinematic/step_056.h5'

f = h5py.File(data,'r')
Np= 132
nramp = 3 
d = np.array(f['Sample Set']).T
Tr = d[2*Np+nramp:3*Np+nramp,:].flatten()
Vr = d[3*Np+nramp:4*Np+nramp,:].flatten()
U_dip = d[:Np,:].flatten()
U_strike = d[Np:2*Np,:].flatten()

Hypo_as=d[4*Np+nramp:4*Np+nramp+1,:].flatten() # along-strike
Hypo_dd=d[4*Np+nramp+1:4*Np+nramp+2,:].flatten() # down-dip

plot_model_hist('Iquique',  U_strike,U_dip, Tr,Vr,Hypo_dd,Hypo_as)


# ### Gorkha ###
data = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Gorkha/model/kinematic/step_012.h5'
f = h5py.File(data,'r')

Np= 162
nramp = 0 
d = np.array(f['Sample Set']).T
Tr = d[2*Np+nramp:3*Np+nramp,:].flatten()
Vr = d[3*Np+nramp:4*Np+nramp,:].flatten()
U_dip = d[:Np+nramp,:].flatten()
U_strike = d[Np+nramp:2*Np+nramp,:].flatten()

Hypo_as=d[4*Np+nramp:4*Np+nramp+1,:].flatten() # along-strike
Hypo_dd=d[4*Np+nramp+1:4*Np+nramp+2,:].flatten() # down-dip
plot_model_hist('Gorkha', U_dip,U_strike, Tr,Vr,Hypo_dd,Hypo_as)
# ### Pedernales ###
data = 'EQ/Pedernales/model/kinematic/step_053.h5'
f = h5py.File(data,'r')

Np= 80
nramp = 9
data = np.array(f['Sample Set']).T


zeros_id = np.where(~data.any(axis=0))[0]
data = np.delete(data,zeros_id,axis=1)
zero_id_tr = np.where(data[2*Np+nramp:3*Np+nramp,:]==0)[1]
data = np.delete(data,zero_id_tr,axis=1)
zero_id_vr = np.where(data[3*Np+nramp:4*Np+nramp,:]==0)[1]
data = np.delete(data,zero_id_vr,axis=1)
zero_id_hyp = np.where(data[4*Np+nramp+1:,:]==0)[1]
data = np.delete(data,zero_id_hyp,axis=1)
zero_id_hyp = np.where(data[4*Np+nramp:4*Np+nramp+1,:]==0)[1]
data = np.delete(data,zero_id_hyp,axis=1)
d = data
Tr = d[2*Np+nramp:3*Np+nramp,:].flatten()
Vr = d[3*Np+nramp:4*Np+nramp,:].flatten()
U_dip = d[:Np+nramp,:].flatten()
U_strike = d[Np+nramp:2*Np+nramp,:].flatten()
Hypo_as=d[4*Np+nramp:4*Np+nramp+1,:].flatten() # along-strike
Hypo_dd=d[4*Np+nramp+1:4*Np+nramp+2,:].flatten() # down-dip
plot_model_hist('Pedernales', U_strike,U_dip, Tr,Vr,Hypo_dd,Hypo_as)

# ### Tohoku ###
dir_file = 'EQ/Tohoku/model/kinematic/step_52/catmip-theta52.bin'
data = np.fromfile(dir_file,'double').reshape((866,int(1e6)))

Np= 24*9
nramp = 0 
d = data
Tr = d[2*Np+nramp:3*Np+nramp,:].flatten()
Vr = d[3*Np+nramp:4*Np+nramp,:].flatten()
U_dip = d[:Np+nramp,:].flatten()
U_strike = d[Np+nramp:2*Np+nramp,:].flatten()
Hypo_as=d[4*Np+nramp:4*Np+nramp+1,:].flatten() # along-strike
Hypo_dd=d[4*Np+nramp+1:4*Np+nramp+2,:].flatten() # down-dip
plot_model_hist('Tohoku', U_strike,U_dip, Tr,Vr,Hypo_dd,Hypo_as)



