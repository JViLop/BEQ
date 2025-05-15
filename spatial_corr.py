# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:18:29 2025

@author: joanv
"""

import matplotlib.pyplot as plt
import numpy
from utils.model_reader import EQ_model


CORR = {}
SCALED_DR ={}
model = EQ_model('Iquique','kinematic','h5',(11,12),17,nramp=3, sampling=True,nsamples=5000)
model.plot_spatial_corr(0.025)
model.spatial_corr_max()
model.corr_matrix()
CORR['Iquique'] = model.corr_wrt_max
SCALED_DR['Iquique'] = model.scaled_dr 



model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0,sampling=True, nsamples=5000)
model.plot_spatial_corr(-0.19)
model.spatial_corr_max()
model.corr_matrix()
CORR['Illapel'] = model.corr_wrt_max
SCALED_DR['Illapel'] = model.scaled_dr
 


model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99,sampling = True, nsamples=5000)
model.plot_spatial_corr(-0.1)
model.spatial_corr_max()
model.corr_matrix()
CORR['Pedernales'] = model.corr_wrt_max
SCALED_DR['Pedernales'] = model.scaled_dr  

model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=5000)
model.plot_spatial_corr(-0.425)
model.spatial_corr_max()
model.corr_matrix()
CORR['Tohoku'] = model.corr_wrt_max
SCALED_DR['Tohoku'] = model.scaled_dr 

model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,'depth':3,'strike':4,'dip':5},nrows_skip=1,nramp=0,rake=107,sampling=True, nsamples = 5000)
model.plot_spatial_corr(-0.32)
model.spatial_corr_max()
model.corr_matrix()
CORR['Gorkha'] = model.corr_wrt_max
SCALED_DR['Gorkha'] = model.scaled_dr 


markers = ['o','D','*','^','s']
labels = {'U':'$U$ (m)',
              'Tr':'$T_r$ (s)',
              'Vr':'$V_r$ (km/s)',
              'U_Tr':'$U/T_r$ (m/s)'}
colors = ['red','blue','green','skyblue','orange']
alphas = [0.9,0.85,0.8,0.75,0.7]
sizes = [9,8,7,6,5]


# all events into 4 subplots

fig, axes = plt.subplots(2,2,figsize=(6,4))
for k,param in enumerate(['U','Tr','Vr','U_Tr']):
    title = labels[param]
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Gorkha','Pedernales']):
    #for i,name in enumerate(['Gorkha']):     
         corr = CORR[name][param] 
         scaled_dr = SCALED_DR[name]
         if k ==0:
             
             
             axes[k//2][k%2].scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i],label=f'{name}')
             axes[k//2][k%2].legend(fontsize=4,ncols = 2,handletextpad =0.3,labelspacing=0.3)
         else:
             axes[k//2][k%2].scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i])
        
         axes[k//2][k%2].axhline(y = 0,linestyle='dashed',color='grey',linewidth = 0.1)

         axes[k//2][k%2].margins(0)        
         axes[k//2][k%2].set_ylim(-1.1,1.1)
         axes[k//2][k%2].set_xlim(-0.05,1.1)
         axes[k//2][k%2].set_yticks([-1.0,-0.5,0,0.5,1])
         axes[k//2][k%2].set_ylabel('Correlation',fontsize=8)
         axes[k//2][k%2].set_xlabel('Distance from maximum-slip patch (km)',fontsize=8)
         # axes[i].set_aspect('equal', 'box')
         axes[k//2][k%2].set_title(title,fontsize=10)
         axes[k//2][k%2].tick_params(labelsize=8)
         
plt.subplots_adjust(wspace=0.3,hspace=0.75)
fig.savefig('spatial_correlations_with_distance.pdf')



# all events into 4 subpanels each with 5 subsubpanel for all events
markers = ['o','D','*','^','s']
labels = {'U':'$U$ (m)',
              'Tr':'$T_r$ (s)',
              'Vr':'$V_r$ (km/s)',
              'U_Tr':'$U/T_r$ (m/s)'}
colors = ['red','blue','green','skyblue','orange']
alphas = [0.9,0.85,0.8,0.75,0.7]
sizes = [12,12,12,12,12]
fig = plt.figure(figsize=(10,12))
gs0 = fig.add_gridspec(2,2,wspace=0.3)

gs00 = gs0[0].subgridspec(5,1,hspace=0.0)
gs01 = gs0[1].subgridspec(5,1,hspace=0.0)
gs02 = gs0[2].subgridspec(5,1,hspace=0.0)
gs03 = gs0[3].subgridspec(5,1,hspace=0.0)


gs = [gs00,gs01,gs02,gs03]
for k,(param,g) in enumerate(zip(['U','Tr','Vr','U_Tr'],gs)):
    title = labels[param]
    
    for i,name in enumerate(['Tohoku','Illapel','Iquique','Pedernales','Gorkha']):
        ax = fig.add_subplot(g[i])
        corr = CORR[name][param] 
        scaled_dr = SCALED_DR[name]
        if i==0 and k==0:
            ax.scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i],label=f'{name}')
            ax.set_ylabel('Correlation',fontsize=8)
            ax.set_xlabel('Normalized distance from maximum-slip patch (km)',fontsize=8)
            # axes[i].set_aspect('equal', 'box')
            ax.legend(fontsize=8,ncols = 2,handletextpad =0.3,labelspacing=0.3)
            ax.set_title(title,fontsize=12)
            ax.tick_params(labelsize=8)
        elif i==0:
            ax.scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i])
            ax.set_title(title,fontsize=10)
        elif i ==4 and k==0:
            ax.scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i],label=f'{name}')
            ax.set_xlabel('Normalized distance from maximum-slip patch (km)',fontsize=8)
            ax.legend(fontsize=8,ncols = 2,handletextpad =0.3,labelspacing=0.3)
        
        elif i==4:
            ax.scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i],label=f'{name}')
            ax.set_xlabel('Normalized distance from maximum-slip patch (km)',fontsize=8)
        elif k==0:
            ax.scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i],label=f'{name}')
            ax.legend(fontsize=8,ncols = 2,handletextpad =0.3,labelspacing=0.3)
        else:
            ax.scatter(scaled_dr, corr,s=sizes[i],marker=markers[i],edgecolors='k',linewidth=0.25,c=colors[i],alpha=alphas[i])

        ax.axhline(y = 0,linestyle='dashed',color='grey',linewidth = 0.1)

        ax.margins(0)        
        ax.set_ylim(-1.2,1.2)
        ax.set_xlim(-0.05,1.1)

        
plt.subplots_adjust(wspace=0.3,hspace=0.2)
fig.savefig('events_spatial_correlations_with_distance.pdf')
