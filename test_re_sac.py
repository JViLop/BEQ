# -*- coding: utf-8 -*-
"""
Created on Mon Mar  4 11:08:00 2024

@author: joanv
"""

import re
import obspy as obs
import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# st_01 = obs.read('testwm/single/p0180_0000_E.SAC')
# st_01.plot()
# st_02 = obs.read('testwm/all/p0180_0000_E.SAC')
# st_02.plot()

## test single patch ###

# dir_ad = 'test_single_patch/adZ.SAC'
# dir_as = 'test_single_patch/asZ.SAC'

# st_all = obs.Stream()
# st_ad = obs.read(dir_ad)
# st_as = obs.read(dir_as)
# # st_ad[0].stats.starttime = st_ad[0].stats.starttime + 6000
# # st_ad[0].trim(starttime = sttime,pad=True,fill_value=0)
# st_ad.plot()
# st_as.plot()
# x1 = st_ad[0].data*1e-2*70
# x2 = st_as[0].data*1e-2*10
# x = x1 +x2
# st_all = st_ad + st_as
# st_all.plot()
# st_all.stack(stack_type = 'linear')
# st_all[0].normalize(norm=1/2)
# st_all.plot()


dir_mean_data = 'Tohoku_mean_kinematic_model.csv'
df = pd.read_csv(dir_mean_data)
shape_geom = (9,24)
Udd =  np.flip(df['U_parallel'].values.reshape(shape_geom,order='F'),axis=0).flatten()
Uas =  np.flip(df['U_perp'].values.reshape(shape_geom,order='F'),axis=0).flatten()

# l = ['p0158_0002_E.SAC', 'p0068_0002_E.SAC', 'p0158_0001_E.SAC']
# x = '0002'
# comp = 'E'
# fmt = re.compile(r'p\d{4}_%s_%s.\w+' % (x, comp))
# result = re.findall(fmt, l[1])


# st1 = obs.read('p0006_0000_N.SAC')
# st2 = obs.read('p0007_0000_N.SAC')


# f = st1[0].stats.sampling_rate

# # st = obs.read('Tohoku_100_full_r_0000_Z_ad.mseed')
# # st[:3].differentiate()
# # st[:3].plot()
# # # st[87].plot()
# # # st[87].stats.starttime = st[87].stats.starttime  + 100
# # # st[87].plot()
# # # st.stack(stack_type='linear')
# # # st[0].normalize(norm=1/len(st))
# # st[65:69].plot()

# st0 = obs.read('p0000_0000_Z.SAC')
# st0.plot()
# print(st)
# st1[0].trim(endtime = st1[0].stats.endtime + 120,pad=True,fill_value=0)
# print(st1)
# st1.plot()
# st0 = obs.read('p0008_0000_E.SAC')
# st0[0].trim(starttime= st0[0].stats.starttime - 10,endtime= st0[0].stats.endtime + 100,pad=True,fill_value=0)
# st0.plot()


# files = os.listdir('wavemod_test')
# st_all = obs.Stream()
# for file in files:
#     file_dir = os.path.join('wavemod_test',file)
#     st = obs.read(file_dir)
#     st_all+=st
# st_all.stack(stack_type='linear')
# st_all[0].normalize(norm=1/len(st_all))
# st_all.integrate().integrate()
# st_all.plot()



# st = obs.read('p0007_0000_N.SAC')



st_global = obs.Stream()
acc = False
vel = False
comps = {'0':'E','1':'N','2':'Z'}
k = 2



for num in [77]:
    st_total = obs.Stream()
    for comp in 'ENZ':
        files = os.listdir('traces_test/traces{}/{}'.format(num,comp))
        st_as = obs.Stream()
        st_dd = obs.Stream()
        file_dir_as = os.path.join('traces_test/traces{}/{}'.format(num,comp),files[1])
        file_dir_dd = os.path.join('traces_test/traces{}/{}'.format(num,comp),files[0])
        st = obs.read(file_dir_as)
        st.sort()
        data_as = []
        for i in range(len(Uas)):
                tr = st[i]
                # tr = st[i].integrate()
                tr.data = tr.data*(1e-2)*Uas[i]
                st_as.append(tr)   
                data_as.append(tr.data)
        
        data_as = np.array(data_as)
        n_as = len(st_as)    
        st_as.stack(stack_type='linear')
        st_as[0].normalize(norm=1/n_as)
        if vel:
            st_as.differentiate()
        elif acc:
            st_as.differentiate().differentiate()  
        st = obs.read(file_dir_dd)
        st.sort()
        data_dd = []
        for i in range(len(Udd)):  
                tr = st[i]
                # tr = st[i].integrate()
                tr.data = tr.data*(1e-2)*Udd[i]
                st_dd.append(tr)     
                data_dd.append(tr.data)
        data_dd = np.array(data_dd) 
        st_dd.stack(stack_type='linear')
        st_dd[0].normalize(norm=1/n_as)
        if vel:
            st_dd.differentiate()
        elif acc:
            st_dd.differentiate().differentiate()
            
        
            
        st_all = st_as + st_dd
        st_all.stack(stack_type = 'linear')
        st_all[0].normalize(norm=1/2)
 

            
        st_total +=st_all
    
    st_total[0].plot()
    st_total[1].plot()
    st_total[2].plot()


# i = 1

# st_total[i].detrend("linear")
# st_total[i].taper(max_percentage=0.05, type='hann')
# st_total[i].filter('bandpass',freqmin=0.001,freqmax=4,corners=2,zerophase=True)  

# # st_total[i].filter('lowpass',freq=4)  
# st_total[i].plot()  
#     st_global.append(st_total[k])
# st_global.plot(outfile = '{}_disp.jpg'.format(comps[str(k)]),automerge=False)
# st_all.stack(stack_type='linear')
# st_all[0].normalize(norm=1/len(st_all))
# st_all.integrate().integrate()
# st_all.plot()
