# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:10:25 2024

@author: joanv
"""

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


def array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()

def proj_ysrc_coords(dip,nrows):
    dip = dip[:nrows]
    proj_dysrc = -patchsize*np.cos(dip*np.pi/180) # in meters
    proj_ysrc = np.zeros_like(proj_dysrc)
    for i in range(len(proj_ysrc)):
        proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
    ysrc = np.flip(proj_ysrc)    
    return ysrc

df = pd.read_csv('Tohoku_mean_kinematic_model.csv')

nrows,ncols = 9,24
shape = (nrows,ncols)
patchsize = 29


Uperp = array_formatter(df['U_perp'].values,shape)
Uparallel =  array_formatter(df['U_parallel'].values,shape)
SLIP = array_formatter(df['Slip'].values,shape)
RAKE = array_formatter(np.arctan2(Uparallel,Uperp)*(180/np.pi),shape)
DIP = array_formatter(df['dip'].values,shape)
DEPTH = array_formatter(df['depth'].values,shape)
DURATION = array_formatter(df['Tr'].values,shape)


xsrc = np.arange((1/2)*patchsize,ncols*patchsize, patchsize)
ysrc = proj_ysrc_coords(DIP,nrows)
Xsrc,Ysrc= np.meshgrid(xsrc,ysrc)
        
xsrc_flat = Xsrc.flatten()
ysrc_flat = Ysrc.flatten()
zsrc_flat = DEPTH  # it mush be in km (not meters)
Xsrcs  = np.array([ysrc_flat,xsrc_flat,zsrc_flat]).T



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

patchsize = 29

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


xstn = np.array([50,0,-50,-100,-150,-200,-250,-300])
ystn = np.array([0,50,100,150,200,250,300,350,400,450,500,550,600])
Xstn ,Ystn = np.meshgrid(xstn,ystn)
xflat,yflat = Xstn.flatten(),Ystn.flatten()
receivers = np.column_stack((xflat,yflat))



for num in [13]:
    st_total = obs.Stream()
    for comp in 'ENZ':
        files = os.listdir('traces_single_src/test{}/{}'.format(num,comp))
        st_all = obs.Stream()
        file_dir_all = os.path.join('traces_single_src/test{}/{}'.format(num,comp),files[0])
        st = obs.read(file_dir_all)
        st.sort()
        data_all = []
        for i in range(len(Uas)):
                tr = st[i]
                # tr = st[i].integrate()
                tr.data = tr.data*(1e-2)
                st_all.append(tr)   
                data_all.append(tr.data)
        
        data_all = np.array(data_all)
        n_all = len(st_all)    
        st_all.stack(stack_type='linear')
        st_all[0].normalize(norm=1/n_all)
        if vel:
            st_all.differentiate()
        elif acc:
            st_all.differentiate().differentiate()  
        
        st_total += st_all
    
    st_total.plot()




xstn = np.array([50,0,-50,-100,-150,-200,-250,-300])
ystn = np.array([0,50,100,150,200,250,300,350,400,450,500,550,600])
Xstn,Ystn = np.meshgrid(xstn,ystn)
xflat,yflat = Xstn.flatten(order='F'),Ystn.flatten(order='F')
xflat_u,yflat_u = Xstn.flatten(),Ystn.flatten()
receivers_u = np.column_stack((xflat_u,yflat_u))
receivers_u = [list(r) for r in receivers_u]
receivers = np.column_stack((xflat,yflat))
receivers = [list(r) for r in receivers]
order = [receivers_u.index(r) for r in receivers ]

ncols,nrows = len(ystn),len(xstn)
G = np.zeros((nrows*ncols,333))

comp = 'N'
for k,r in enumerate(receivers):
    
    
    x,y = r[0],r[1]
    receiver = np.array([x,y])

    # for comp in 'ZRT':
    #     strm = obs.read(f'traces_{comp}.mseed')
    #     traces = 0
    #     for i in range(216):
    #         traces += strm[i].data/100
            
    #     plt.plot(traces,label=comp)
    #     plt.legend()
    

    path = f'./traces_single_src/all/{comp}/'
    i = order[k]

    strm = obs.read(path + 'Tohoku_100_full_r_0000_%s_%.3i_all.mseed'%(comp,i))
    strm.sort()
    st = obs.Stream()
    data_all = []
    for i in range(len(strm)):
            tr = strm[i]
            # tr = st[i].integrate()
            tr.data = tr.data*(1e-2)
            st.append(tr)   
            data_all.append(tr.data)
    
    dt = tr.stats.delta
    n = tr.stats.npts
    time = np.arange(0,n*dt,dt)
    data_all = np.array(data_all)
    n_all = len(st)    
    st.stack(stack_type='linear')
    st[0].normalize(norm=1/n_all) 
    d = st[0].data
    G[k,:] = d
    



G = G.reshape(nrows,ncols,333)
import numpy as np
from matplotlib import pyplot as plt, animation
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = True


G0 = G[:,:,::10]
nframes = G0.shape[2]
fig, ax = plt.subplots()
xstn = np.array([50,0,-50,-100,-150,-200,-250,-300])
ystn = np.array([0,50,100,150,200,250,300,350,400,450,500,550,600])

cax = ax.pcolormesh(ystn, xstn, G0[:, :, 0], vmin = min(np.min(G0),-np.max(G0)),vmax =  max(-np.min(G0),np.max(G0)),cmap='bwr')
ax.set_title(f'{comp} displacement')
fig.colorbar(cax)

def animate(i):
    cax.set_array(G0[:, :, i].flatten())

anim = animation.FuncAnimation(fig, animate,frames =nframes,interval = 50)
anim.save(filename=f"{comp}_singlesrc.mp4", writer="ffmpeg")
plt.show()
    
 
