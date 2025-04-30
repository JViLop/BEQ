# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:11:57 2025

@author: joanv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from pylab import figure,setp



import h5py

import os 
from matplotlib.colors import TwoSlopeNorm

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']

nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patches = [29,17,18,10,15]
geoms = [(nrows[i],ncols[i]) for i in range(len(names))]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,90,90,107,360-99]
ramps = [0,3,0,0,9]
factors =[3,3,3,2,2]

def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,nramps,factors):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['nramp'] = nramps[i]
        model[name]['factor'] = factors[i]
    return model

def set_stn(c,factor,patch,control=0):
      
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn
  
def proj_ysrc_coords(patch,dip):
      proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
      proj_ysrc = np.zeros_like(proj_dysrc)
      for i in range(len(proj_ysrc)):
          proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
      ysrc = np.flip(proj_ysrc)
        
      return ysrc
       
        
models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,ramps,factors)

observable = 'displacement'
nsamples = 100
model_type = 'kinematic'

working_dir = os.getcwd()  




# for q,name in enumerate(['Tohoku']):
#     print(q)
#     working_dir = os.getcwd()        
#     mean_csv_file_dir = os.path.join(working_dir,f'{name}_mean_{model_type}_model.csv')
#     df = pd.read_csv(mean_csv_file_dir)
#     df= df.drop(df.columns[0],axis=1)    
    
#     nrows, ncols = models[name]['geom'][0], models[name]['geom'][1]
#     patch = models[name]['patch']
#     nparam = models[name]['nparam']
    
       

#     xstn,ystn = set_stn(2,4,patch,control=1)
    
#     nrows_d,ncols_d = len(ystn),len(xstn)
#     Xstn,Ystn = np.meshgrid(xstn,ystn)
#     xstn_flat,ystn_flat =  Xstn.flatten(),Ystn.flatten()

#     h5file_dir = os.path.join(working_dir,f'EDKS_{name}_{observable}_nsamples_{nsamples}.h5')
    
#     f = h5py.File(h5file_dir,'r')

#     dset = np.array(f[observable])
#     # covariance matrix

#     mean = np.mean(dset,axis=0)
#     corr = np.corrcoef(dset.transpose())
#     nparameters = dset.shape[1]
    
#     mean_dx = mean[:nparameters//3]
#     mean_dy = mean[nparameters//3:2*nparameters//3]
#     mean_dz = mean[2*nparameters//3:]
    
#     mean_d = {'x':mean_dx,'y':mean_dy,'z':mean_dz}
    
   
#     # covariance of each observable (x,y,z) 
    

#     corr1 = corr[:nparameters//3,:nparameters//3]
#     corr2 = corr[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
#     corr3 = corr[2*nparameters//3:,2*nparameters//3:]
    
#     corr_d = {'x':corr1, 'y':corr2,'z':corr3}

#     std_d = {'x':np.sqrt(np.diag(corr1)), 'y':np.sqrt(np.diag(corr2)),'z':np.sqrt(np.diag(corr3))}

    
#     d = np.sqrt(mean_dx**2 + mean_dy**2 + mean_dz**2)
#     max_d_patch = np.argmin(abs(d-max(d)))
#     d = np.reshape(d,(nrows_d,ncols_d))
#     col0_max_d,row0_max_d = max_d_patch%ncols_d,max_d_patch//ncols_d
    
#     d_tn = np.reshape(mean_d['y'],(nrows_d,ncols_d))[:,col0_max_d]
#     d_z = np.reshape(mean_d['z'],(nrows_d,ncols_d))[:,col0_max_d]
#     err_tn = np.reshape(std_d['y'],(nrows_d,ncols_d))[:,col0_max_d]
#     err_z =  np.reshape(std_d['y'],(nrows_d,ncols_d))[:,col0_max_d]
    
#     ystn_d,dy,err_dy,dz,err_dz = ystn, d_tn,err_tn,d_z,err_z


names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']

nrows = [9,11,10,9,8]
ncols = [24,12,17,18,10]
patches = [29,17,18,10,15]
geoms = [(nrows[i],ncols[i]) for i in range(len(names))]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,90,90,107,360-99]
ramps = [0,3,0,0,9]
factors =[3,3,3,2,2]
strikes = [194,-13.58,4,293,27.05]
scales = [0.5,0.8e-1,1e-1,0.75e-1,0.75e-1]
shrinks = [0.65,0.6,0.6,0.6,0.65]
def model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,nramps,factors,strikes,scales,shrinks):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i]
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i]
        model[name]['nramp'] = nramps[i]
        model[name]['factor'] = factors[i]
        model[name]['factor'] = factors[i]
        model[name]['strike'] = strikes[i]
        model[name]['scale'] = scales[i]
        model[name]['shrink'] = shrinks[i]
    return model



models = model_dict(names,geoms,patches,arrow_sizes,nparams,rakes,ramps,factors,strikes,scales,shrinks)
name = 'Tohoku'
model_type = 'kinematic'
samples = 100
offset = 75

# file_m = f'{name}_{samples}_{model_type}_curve_m.csv'
# df_m = pd.read_csv(file_m)
# ysrc,m,err_m  = df_m['ysrc'].values,df_m['m'].values, df_m['err_m'].values

# file_disp = f'{name}_{samples}_{model_type}_curve_displacement.csv'
# df_d = pd.read_csv(file_disp)
# ystn_d,dy,err_dy,dz,err_dz = df_d['ystn_d'].values,df_d['tn'].values, df_d['err_tn'].values,df_d['v'].values, df_d['err_v'].values

# file_stress = f'{name}_{samples}_{model_type}_curve_stress.csv'
# df_s = pd.read_csv(file_stress)
# ystn_s,sn,err_sn,sy,err_sy = df_s['ytn_s'].values,df_s['n'].values, df_s['err_n'].values,df_s['ad'].values, df_s['err_ad'].values
label_stress = ['normal','along-dip']
label_displacement = ['trench-normal','vertical']
color_stress = ['deepskyblue','blue']
color_displacement = ['red','deeppink']
lns = ['solid','solid']
markers = ['o','s']

fig = plt.figure(figsize=(20,10))
gs0 = fig.add_gridspec(2, 3,wspace=0.3,hspace=0.2)

gs00 = gs0[0].subgridspec(2, 1,height_ratios=[1,1],hspace=0)
gs01 = gs0[1].subgridspec(2,1,height_ratios=[1,1],hspace=0)
gs02 = gs0[2].subgridspec(2,1,height_ratios=[1,1],hspace=0)
gs10 = gs0[3].subgridspec(2,1,height_ratios=[1,1],hspace=0)
gs11 = gs0[4].subgridspec(2,1,height_ratios=[1,1],hspace=0)

grs = [gs00,gs01,gs02,gs10,gs11]
nsamples = 100
for k,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        if name=='Pedernales':
          samples = 5000
        else:
          samples = 100
          
        file_m = f'summary_curves/{name}_{samples}_{model_type}_curve_m.csv'
        df_m = pd.read_csv(file_m)
        ysrc,m,err_m  = df_m['ysrc'].values,df_m['m'].values, df_m['err_m'].values
        
        
      
        mean_csv_file_dir = os.path.join(working_dir,f'{name}_mean_{model_type}_model.csv')
        df = pd.read_csv(mean_csv_file_dir)
        df= df.drop(df.columns[0],axis=1)    
        
        nrows, ncols = models[name]['geom'][0], models[name]['geom'][1]
        patch = models[name]['patch']
        nparam = models[name]['nparam']
        
           

        xstn,ystn = set_stn(2,4,patch,control=1)
        
        nrows_d,ncols_d = len(ystn),len(xstn)
        Xstn,Ystn = np.meshgrid(xstn,ystn)
        xstn_flat,ystn_flat =  Xstn.flatten(),Ystn.flatten()

        h5file_dir = os.path.join(working_dir,f'EDKS_{name}_{observable}_nsamples_{nsamples}.h5')
        
        f = h5py.File(h5file_dir,'r')

        dset = np.array(f[observable])
        # covariance matrix

        mean = np.mean(dset,axis=0)
        cov = np.cov(dset.transpose())
        nparameters = dset.shape[1]
        
        mean_dx = mean[:nparameters//3]
        mean_dy = mean[nparameters//3:2*nparameters//3]
        mean_dz = mean[2*nparameters//3:]
        
        mean_d = {'x':mean_dx,'y':mean_dy,'z':mean_dz}
        
       
        # covariance of each observable (x,y,z) 
        

        cov1 = cov[:nparameters//3,:nparameters//3]
        cov2 = cov[nparameters//3:2*nparameters//3,nparameters//3:2*nparameters//3]
        cov3 = cov[2*nparameters//3:,2*nparameters//3:]
        
        cov_d = {'x':cov1, 'y':cov2,'z':cov3}

        std_d = {'x':np.sqrt(np.diag(cov1)), 'y':np.sqrt(np.diag(cov2)),'z':np.sqrt(np.diag(cov3))}

        

        slip = np.flip(df['Slip'].values.reshape(nrows,ncols,order='F'),axis=0).flatten()
        xS = np.arange(patch/2 , ncols*patch,patch)
        yS = np.arange(-(nrows-1/2)*patch,0,patch)
           # shift accordingly at surface
        yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
        XS, YS = np.meshgrid(xS,yS)
        xsrc_flat,ysrc_flat = XS.flatten(),YS.flatten() 
        max_slip_patch = np.argmin(abs(slip-max(slip)))
        
        x_target = xsrc_flat[max_slip_patch]
        y_target = ysrc_flat[max_slip_patch]

        r_target = np.array([x_target,y_target])

        r_all = np.column_stack((xstn_flat,ystn_flat))
        dr = np.sqrt((r_all[:,0]  - r_target[0])**2 + (r_all[:,1]  - r_target[1])**2)
        id_max_patch_primed = np.argmin(dr)
        
        d = np.sqrt(mean_dx**2 + mean_dy**2 + mean_dz**2)
        #max_d_patch = np.argmin(abs(d-max(d)))
        max_d_patch = id_max_patch_primed
        d = np.reshape(d,(nrows_d,ncols_d))
        col0_max_d,row0_max_d = max_d_patch%ncols_d,max_d_patch//ncols_d
        
        print(col0_max_d,row0_max_d)
        d_tn = np.reshape(mean_d['y'],(nrows_d,ncols_d))[:,col0_max_d]
        d_z = np.reshape(mean_d['z'],(nrows_d,ncols_d))[:,col0_max_d]
        err_tn = np.reshape(std_d['y'],(nrows_d,ncols_d))[:,col0_max_d]
        err_z =  np.reshape(std_d['z'],(nrows_d,ncols_d))[:,col0_max_d]
        
        if name == 'Gorkha':
            ystn-=offset
        ystn_d,dy,err_dy,dz,err_dz = ystn, d_tn,err_tn,d_z,err_z
        
        #file_disp = f'summary_curves/{name}_{samples}_{model_type}_curve_displacement.csv'
        #df_d = pd.read_csv(file_disp)
        #ystn_d,dy,err_dy,dz,err_dz = df_d['ystn_d'].values,df_d['tn'].values, df_d['err_tn'].values,df_d['v'].values, df_d['err_v'].values

        file_stress = f'summary_curves/{name}_{samples}_{model_type}_curve_stress.csv'
        df_s = pd.read_csv(file_stress)
        ystn_s,sn,err_sn,sy,err_sy = df_s['ytn_s'].values,df_s['n'].values, df_s['err_n'].values,df_s['ad'].values, df_s['err_ad'].values

        axprops = dict()
        ax1 = fig.add_subplot(grs[k][0,0],**axprops)

        for n,(d,errd) in enumerate(zip([dy,dz],[err_dy,err_dz])):

          ax1.errorbar(ystn_d,d,yerr = errd,capthick = 0.2,capsize=2,linestyle = lns[n],marker=markers[n],markeredgecolor='black',markeredgewidth = 0.2,linewidth=1.25,ms=2, color=color_displacement[n],label=label_displacement[n])
          if name !='Gorkha':
              ax1.axvline(x = 0,ls='--',linewidth=0.6)

          ax1.axhline(y = 0,ls='--',linewidth=0.4,color='red')
          ax1.set_ylabel('Displacement (m)',fontsize=9)
          ax1.tick_params(axis='y',labelsize=9,labelcolor='red')
          ax1.legend(fontsize=7.5,loc='upper left',handlelength= 1.5,markerscale = 1.0)
          ax1.text(-0.1,1.15,'abcde'[k],transform=ax1.transAxes,fontsize=18,fontweight='bold')
          axprops['sharex'] = ax1
          #ax1.tick_params(axis='x',labelsize=0)

          # force x axes to remain in register, even with toolbar navigatio
        ax2 = ax1.twinx()
        for n,(s,errs) in enumerate(zip([sn,sy],[err_sn,err_sy])):
          # ax2.plot(self.ystn_s*1e-3,self.s_norm(i)['s'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_stress[i],label=label_stress[i])
          # ax2.fill_between(self.ystn_s*1e-3,self.s_norm(i)['s_down'],self.s_norm(i)['s_up'],facecolor=color_stress[i],alpha=0.2)
          ax2.errorbar(ystn_s*1e-3,s,yerr = errs,capthick = 0.2,capsize=2,linestyle = lns[n],marker=markers[n],markeredgecolor='black',markeredgewidth = 0.2,linewidth=1.25,ms=2, color=color_stress[n],label=label_stress[n])
          ax2.set_ylabel('Stress Change (MPa)',fontsize=9)
          if name !='Gorkha':
              ax2.axvline(x = 0,ls='--',linewidth=0.6)
          ax2.axhline(y = 0,ls='--',linewidth=0.4,color='blue')
          ax2.legend(fontsize=7.5,loc='upper right',handlelength= 1.5,markerscale = 1.0 )
          ax2.tick_params(axis='y',labelsize=9,labelcolor='blue')
          #ax2.tick_params(axis='x',labelsize=0)

        ax3 = fig.add_subplot(grs[k][1,0],**axprops)
        # ax3.errorbar(self.ysrc*1e-3,self.m()['m'],yerr = self.m()['err_m'],drawstyle='steps-mid',capthick = 0.2,capsize=1.25,marker='.',linewidth=0.6,ms=2, color='black')
        y0 = ysrc*1e-3
        dysrcleft = np.abs(y0[0] - y0[1])/2
        dysrcright = np.abs(y0[-1] -y0[-2])/2
        ysrcleft = np.array([y0[0]-dysrcleft,y0[0]])
        ysrcright = np.array([y0[-1],y0[-1] + dysrcright])
        mleft = np.array([m[0],m[0]])
        mright = np.array([m[-1],m[-1]])


    # Create the step function
        ax3.text(0.025,0.825,f'{name}',fontweight='bold',transform=ax3.transAxes)
        ax3.plot(ysrcleft, mleft,linewidth=0.8,color='black')
        ax3.plot(ysrcright, mright,linewidth=0.8,color='black')
        ax3.fill_between(ysrcleft, mleft-err_m[0],mleft+err_m[0],step='pre',facecolor='grey',alpha=0.75)
        ax3.fill_between(ysrcright, mright-err_m[-1],mright+err_m[-1],step='post',facecolor='grey',alpha=0.75)
        ax3.plot(ysrc*1e-3,m,drawstyle='steps-mid',marker = '.',linewidth=1.25,ms=2.5, color='black')
        ax3.fill_between(ysrc*1e-3,m + err_m ,m - err_m,step='mid',facecolor='grey',alpha=0.75)
        if name !='Gorkha':
          ax3.axvline(x = 0,ls='--',linewidth=0.6)
        ax3.axhline(y = 0,ls='--',linewidth=0.4,color='black')
        ax3.set_ylabel('Slip (m) ',fontsize=9)
        ax3.annotate('Trench',(0,0.5),(2,0.5),fontsize=9,rotation = 90)
        ax3.set_xlabel('Distance from trench (km)',fontsize=9)
        if name !='Gorkha':
            ax3.axvline(x = 0,ls='--',linewidth=0.6)
        ax3.axhline(y = 0,ls='--',linewidth=0.4,color='black')
        ax3.set_ylabel('Slip (m) ',fontsize=9)
        if name =='Gorkha':
            ax3.text(7.5-offset, 3, "Trench",
            ha="center", va="center", rotation=0, size=8,
            bbox=dict(boxstyle="rarrow,pad=0.3",
                      fc="white", ec="black", lw=0.6))
            # ax3.annotate((0,0.5),(2,0.5),ha="center", va="center", rotation=0, size=4,
            # bbox=dict(boxstyle="rarrow,pad=0.3",fc="blue", ec="black", lw=2))
        else:
            ax3.annotate('Trench',(0,0.5),(2,0.5),fontsize=9,rotation = 90)
        ax3.set_xlabel('Distance from trench (km)',fontsize=9)
        ax3.tick_params(labelsize=9)
        for ax in ax1, ax2:
          setp(ax.get_xticklabels(), visible=False)

fig.savefig('summary_curves_all_events.pdf')

