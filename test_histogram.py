# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:07:00 2024

@author: joanv
"""

import h5py 
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import matplotlib.patches as m_patches
name = 'Tohoku'
parameters = 866
samples  = 100
shape = (9,24)
th = 0.2

def ave_potency(name,parameters,shape,samples,dir_folder= 'histograms',rake=90, RotAngle=None,th=None):
    bin_data = np.fromfile(f'{dir_folder}/{name}_kinematic_n_{samples}.dat','float').reshape((parameters,int(samples)))

    # data = h5py.File('Tohoku_Stress change_nsamples_100.h5')

    # npatches = data['Stress change'].shape[1]//3
    # stress = data['Stress change']

    data = h5py.File(f'{dir_folder}/{name}_Stress_change_nsamples_{samples}.h5')
    #data = h5py.File(f'{dir_folder}/{name}_stress_nsamples_{samples}.h5')


    npatches = data['Stress change'].shape[1]//3
    stress = data['Stress change']
    
    # just valid for Pedernales 
    
    geometry_d = pd.read_csv('INPUT/Pedernales/model/kinematic/all_samples/mean/Pedernales_mean_kinematic_model.csv')
    # av_stress_dip = []
    # av_stress_combined = []
    av_stress_amplitude = []
    for i in range(samples):
        
        Udip = bin_data[npatches:2*npatches,i]
        Ustrike = bin_data[:npatches,i]
        if name =='Gorkha':
              Udip = bin_data[:npatches,i] 
              Ustrike = bin_data[npatches:2*npatches,i]
        if name == 'Pedernales':
            ss = np.zeros(npatches)
            ds = np.zeros(npatches)
            strike = geometry_d['strike']*((np.pi)/180)
            dip = geometry_d['dip']*((np.pi)/180)
            for p in range(npatches):
                RotAngle2 = RotAngle*((np.pi) / 180.)
                
                rotation = np.arctan2(np.tan(strike[p]) - np.tan(RotAngle2), 
                                    np.cos(dip[p])*(1.+np.tan(RotAngle2)*np.tan(strike[p])))
    
                # If RotAngle within ]90, 270], change rotation
                if RotAngle > 90. and RotAngle<=270.:
                    rotation += np.pi
    
                rp = Ustrike[p]    # rake-perpendicular
                ar = Udip[p] # rake-parallel
    
                ss[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                ds[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
            Ustrike = ss
            Udip= ds    
        else:
            rake0 = (rake - 90)*(np.pi/180)
            Ustrike =  Ustrike*np.cos(rake0)  - Udip*np.sin(rake0) 
            Udip = Ustrike*np.sin(rake0)  + Udip*np.cos(rake0)      
            
        Udip = np.flip(Udip.reshape(shape,order='F'),axis=0).flatten()
        Ustrike = np.flip(Ustrike.reshape(shape,order='F'),axis=0).flatten()
        U = np.sqrt(Udip**2 + Ustrike**2)
        
    
        Umax = max(U)
        ind_Uth = np.where(U>th*Umax)
        Ustrike = Ustrike[ind_Uth]
        Udip = Udip[ind_Uth]
        
        dip_stress = -stress[i,2*npatches:] 
        strike_stress = -stress[i,npatches:2*npatches] 
        dip_stress = dip_stress[ind_Uth]/(2*30e3)  # added to account for potency
        strike_stress = strike_stress[ind_Uth]/(2*30e3)
        
        dot_prod_dip = np.dot(Udip,dip_stress)
        dot_prod_strike = np.dot(Ustrike,strike_stress)
        # denominator = np.sum(Udip)
        denominator_amplitude = np.sqrt((np.sum(Udip))**2 + (np.sum(Ustrike))**2)
        # stress_drop_dip = dot_prod_dip/denominator
        # stress_drop_combined = (dot_prod_dip + dot_prod_strike)/denominator
        
        try:
            
            if th==0:
                stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
            else:
                # uncomment below for spatial-average 
                #stress_drop_amplitude = np.sqrt((np.dot(dip_stress,np.ones_like(dip_stress)))**2 + (np.dot(strike_stress,np.ones_like(strike_stress)))**2)/len(Ustrike)
                stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude

        # av_stress_dip.append(stress_drop_dip)
        # av_stress_combined.append(stress_drop_combined)
            av_stress_amplitude.append(stress_drop_amplitude)
            
        except:
            continue
 
    return av_stress_amplitude



def ave_stress(name,parameters,shape,samples,dir_folder= 'histograms',rake=90, RotAngle=None,th=None):
    bin_data = np.fromfile(f'{dir_folder}/{name}_kinematic_n_{samples}.dat','float').reshape((parameters,int(samples)))

    # data = h5py.File('Tohoku_Stress change_nsamples_100.h5')

    # npatches = data['Stress change'].shape[1]//3
    # stress = data['Stress change']

    data = h5py.File(f'{dir_folder}/{name}_Stress_change_nsamples_{samples}.h5')
    #data = h5py.File(f'{dir_folder}/{name}_stress_nsamples_{samples}.h5')


    npatches = data['Stress change'].shape[1]//3
    stress = data['Stress change']
    
    # just valid for Pedernales 
    
    geometry_d = pd.read_csv('INPUT/Pedernales/model/kinematic/all_samples/mean/Pedernales_mean_kinematic_model.csv')
    # av_stress_dip = []
    # av_stress_combined = []
    av_stress_amplitude = []
    for i in range(samples):
        
        Udip = bin_data[npatches:2*npatches,i]
        Ustrike = bin_data[:npatches,i]
        if name =='Gorkha':
              Udip = bin_data[:npatches,i] 
              Ustrike = bin_data[npatches:2*npatches,i]
        if name == 'Pedernales':
            ss = np.zeros(npatches)
            ds = np.zeros(npatches)
            strike = geometry_d['strike']*((np.pi)/180)
            dip = geometry_d['dip']*((np.pi)/180)
            for p in range(npatches):
                RotAngle2 = RotAngle*((np.pi) / 180.)
                
                rotation = np.arctan2(np.tan(strike[p]) - np.tan(RotAngle2), 
                                    np.cos(dip[p])*(1.+np.tan(RotAngle2)*np.tan(strike[p])))
    
                # If RotAngle within ]90, 270], change rotation
                if RotAngle > 90. and RotAngle<=270.:
                    rotation += np.pi
    
                rp = Ustrike[p]    # rake-perpendicular
                ar = Udip[p] # rake-parallel
    
                ss[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                ds[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
            Ustrike = ss
            Udip= ds    
        else:
            rake0 = (rake - 90)*(np.pi/180)
            Ustrike =  Ustrike*np.cos(rake0)  - Udip*np.sin(rake0) 
            Udip = Ustrike*np.sin(rake0)  + Udip*np.cos(rake0)      
            
        Udip = np.flip(Udip.reshape(shape,order='F'),axis=0).flatten()
        Ustrike = np.flip(Ustrike.reshape(shape,order='F'),axis=0).flatten()
        U = np.sqrt(Udip**2 + Ustrike**2)
        
    
        Umax = max(U)
        ind_Uth = np.where(U>th*Umax)
        Ustrike = Ustrike[ind_Uth]
        Udip = Udip[ind_Uth]
        
        dip_stress = -stress[i,2*npatches:] 
        strike_stress = -stress[i,npatches:2*npatches] 
        dip_stress = dip_stress[ind_Uth]
        strike_stress = strike_stress[ind_Uth]
        
        dot_prod_dip = np.dot(Udip,dip_stress)
        dot_prod_strike = np.dot(Ustrike,strike_stress)
        # denominator = np.sum(Udip)
        denominator_amplitude = np.sqrt((np.sum(Udip))**2 + (np.sum(Ustrike))**2)
        # stress_drop_dip = dot_prod_dip/denominator
        # stress_drop_combined = (dot_prod_dip + dot_prod_strike)/denominator
        
        try:
            if th==0:
                stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
            else:
                stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
                #stress_drop_amplitude = np.sqrt((np.dot(dip_stress,np.ones_like(dip_stress)))**2 + (np.dot(strike_stress,np.ones_like(strike_stress)))**2)/len(Ustrike)
        # av_stress_dip.append(stress_drop_dip)
        # av_stress_combined.append(stress_drop_combined)
            av_stress_amplitude.append(stress_drop_amplitude)
            
        except:
            continue
 
    return av_stress_amplitude

def plot_stress_hist(*args,rake=90,RotAngle =None,thr = [0.1,0.2],bin_lims=[12,20]):
    # av_stress =  ave_stress(name,parameters,samples,shape)
    # av10_stress = ave_stress(name,parameters,samples,shape,th = thrs[0])
    # av20_stress = ave_stress(name,parameters,samples,shape,th = thrs[1])
    nbins = np.arange(bin_lims[0],bin_lims[1]+0.05,0.05)
    av_stress =  ave_stress(*args,th=0,rake=rake,RotAngle=RotAngle)
    av10_stress = ave_stress(*args,th = thr[0],rake=rake,RotAngle=RotAngle)
    av20_stress = ave_stress(*args,th = thr[1],rake=rake,RotAngle=RotAngle)
    
    
    fig, ax = plt.subplots(figsize=(6,4),dpi=900)
    # the histogram of the data
    # n, bins, patches = axes[0].hist(av_stress_dip, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \Delta u_{dip}\cdot\Delta \sigma_{dip}\,dS}{\int \Delta u_{dip} \,dS}$')
    # n, bins, patches = axes[1].hist(av_stress_combined, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,dS}{\int \Delta u_{dip} \,dS}$')
    n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.6, color='blue',ec = 'black',lw=1,density=False,label=r'$\Sigma$')
    n, bins, patches = ax.hist(av10_stress, nbins,color='red', alpha = 0.25, ec = 'black',lw=0.4,density=False,label=r'$\Sigma_{0.1}$')
    n, bins, patches = ax.hist(av20_stress, nbins,color='orange',alpha = 0.25,ec = 'black',lw=0.4, density=False,label=r'$\Sigma_{0.2}$')
    # n, bins, patches = ax.hist(av10_stress, nbins,histtype = 'step',fill = False, alpha = 0.4, ec = 'red',lw=1,density=False,label=r'$\Sigma_{0.1}$')
    # n, bins, patches = ax.hist(av20_stress, nbins,histtype = 'step',fill=False,alpha = 0.4,ec = 'black',lw=1, density=False,label=r'$\Sigma_{0.2}$')
    
    ax.set_title(f'{args[0]} Stress drop',fontweight='bold',fontsize=13)
    # add a 'best fit' line
    ax.set_xlabel('stress drop (MPa)')
    
    ax.set_ylabel('Count')
    ax.set_xlabel('Stress drop (MPa)')
    ax.legend(fontsize=11)

    # ax.legend(title=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,d\Sigma}{|\int \mathbf{\Delta u}\,d\Sigma|}$')
    plt.tight_layout()
    plt.savefig(f'histograms/{args[0]}_av_stress_histograms_{args[3]}.png')
    plt.show()

# the name of h5 files should be changed as `stress` only
# plot_stress_hist('Tohoku', 866,(9,24),1000)
# plot_stress_hist('Gorkha', 650,(9,18),1000,rake=107,bin_lims = [2,8])
# plot_stress_hist('Iquique', 533,(11,12),1000,bin_lims = [5,13])
# plot_stress_hist('Illapel', 682,(10,17),1000,bin_lims = [5,17])
# plot_stress_hist('Pedernales', 331,(8,10),1000,RotAngle=360-99,bin_lims = [1,5])



 


# EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges':[12e3/(2*30),18.5e3/(2*30)],'rake':90,'RotAngle':None,'color':'blue','loc':(13.2e3/(2*30),825)},
#        'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges':[4e3/(2*30),8e3/(2*30)],'rake':107,'RotAngle':None,'color':'red','loc':(14.2e3/(2*30),825)},
#        'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges':[5e3/(2*30),13e3/(2*30)],'rake':90,'RotAngle':None,'color':'orange','loc':(15.3e3/(2*30),825)},
#        'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges':[5e3/(2*30),17e3/(2*30)],'rake':90,'RotAngle':None,'color':'seagreen','loc':(16.45e3/(2*30),825)},
#        'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges':[4e3/(2*30),10e3/(2*30)],'rake':None,'RotAngle':360-99,'color':'brown','loc':(17.3e3/(2*30),825)}}


EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges':[1,18.5e3/(2*30)],'rake':90,'RotAngle':None,'color':'blue','loc':(270,2300)},
       'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges':[0,8e3/(2*30)],'rake':107,'RotAngle':None,'color':'red','loc':(300,2300)},
       'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges':[0,13e3/(2*30)],'rake':90,'RotAngle':None,'color':'orange','loc':(330,2300)},
       'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges':[0,17e3/(2*30)],'rake':90,'RotAngle':None,'color':'seagreen','loc':(360,2300)},
       'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges':[0,10e3/(2*30)],'rake':None,'RotAngle':360-99,'color':'brown','loc':(385,2300)}}

EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges':[1,18.5e3/(2*30)],'rake':90,'RotAngle':None,'color':'blue','loc':(280,825)},
       'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges':[0,8e3/(2*30)],'rake':107,'RotAngle':None,'color':'red','loc':(310,825)},
       'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges':[0,13e3/(2*30)],'rake':90,'RotAngle':None,'color':'orange','loc':(338,825)},
       'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges':[0,17e3/(2*30)],'rake':90,'RotAngle':None,'color':'seagreen','loc':(365,825)},
       'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges':[0,10e3/(2*30)],'rake':None,'RotAngle':360-99,'color':'brown','loc':(387,825)}}



def plot_all_hist_potency(EQs,samples= 1000,thr = [0.1,0.2],folder = 'stress_samples_5000'):
    # av_stress =  ave_stress(name,parameters,samples,shape)
    # av10_stress = ave_stress(name,parameters,samples,shape,th = thrs[0])
    # av20_stress = ave_stress(name,parameters,samples,shape,th = thrs[1])
    
    values = {'all':{},'10':{},'20':{}}
    stds =  {'all':{},'10':{},'20':{}}
    fig, ax = plt.subplots(figsize=(10,5),dpi=900)
    for name in list(EQs.keys()):
        nbins = np.arange(EQs[name]['edges'][0],EQs[name]['edges'][1]+0.1e3/(2*30),0.1e3/(2*30))
        av_stress =  ave_potency(name,EQs[name]['n'],EQs[name]['shape'],samples,th=0,rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder)
        av10_stress = ave_potency(name,EQs[name]['n'],EQs[name]['shape'],samples,th = thr[0],rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder)
        av20_stress = ave_potency(name,EQs[name]['n'],EQs[name]['shape'],samples,th = thr[1],rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder)
        
        av_stress = np.array(av_stress)*1e6
        av10_stress = np.array(av10_stress)*1e6
        av20_stress = np.array(av20_stress)*1e6
        values['all'][name] = np.nanmean(np.array(av_stress))
        values['10'][name] = np.nanmean(np.array(av10_stress))
        values['20'][name] = np.nanmean(np.array(av20_stress))
        
        stds['all'][name] = np.nanstd(np.array(av_stress))
        stds['10'][name] = np.nanstd(np.array(av10_stress))
        stds['20'][name] = np.nanstd(np.array(av20_stress))
        
  
        if name == 'Pedernales':

            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label='$\overline{\Delta\sigma}_{E}$')
            #n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label = '$\overline{\Delta\sigma}_{20\%}$')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label = '$\overline{\Delta\sigma}_{10\%}$')
            
            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\epsilon}_{E}$')
            n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = EQs[name]['color'],histtype='stepfilled',lw=0.2,density=False,label = '$\overline{\epsilon}_{10\%}$')
            n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label = '$\overline{\epsilon}_{20\%}$')
        else:
            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label=' ')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label=' ')
            #n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label=' ')


            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label=' ')
            n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, histtype='stepfilled',ec = EQs[name]['color'],lw=0.2,density=False,label=' ')
            n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label=' ')

            
        # n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = 'black',lw=0.25,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'black',lw=0.2, density=False,label=r'$\Sigma_{0.2}$')
        # n, bins, patches = ax.hist(av10_stress, nbins,histtype = 'step',fill = False, alpha = 0.4, ec = 'red',lw=1,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,histtype = 'step',fill=False,alpha = 0.4,ec = 'black',lw=1, density=False,label=r'$\Sigma_{0.2}$')
        
        ax.text(EQs[name]['loc'][0],EQs[name]['loc'][1],name,fontsize = 7,fontweight='bold')
        # add a 'best fit' line
        ax.set_xlabel('Potency Density (microstrain)')
        
        ax.set_ylabel('Count')
        ax.set_xlim([2e3/(2*30),26e3/(2*30)])
        ax.legend(frameon=False,handlelength=2.750,handleheight=1.25,handletextpad = 0.2,columnspacing =0.5,labelspacing=0.8,ncols=5,loc="upper right", fontsize=10)
        
        # the histogram of the data
        # n, bins, patches = axes[0].hist(av_stress_dip, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \Delta u_{dip}\cdot\Delta \sigma_{dip}\,dS}{\int \Delta u_{dip} \,dS}$')
        # n, bins, patches = axes[1].hist(av_stress_combined, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,dS}{\int \Delta u_{dip} \,dS}$')
        
    '''
        if name == 'Pedernales':

            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label='$\bar{\Delta\sigma}_{E}$')
            n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label = '$\bar{\Delta\sigma}_{20\%}$')
            n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label = '$\bar{\Delta\sigma}_{10\%}$')
            # n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'none',lw=0.2, density=False,label = '$\overline{\Delta\sigma}_{20\%}$')
           # n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = 'none',lw=0.2,density=False,label = '$\overline{\Delta\sigma}_{10\%}$')
           # n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\Delta\sigma}_{E}$')
        else:
            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label=' ')
            n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label=' ')
            n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label=' ')


            #n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'none',lw=0.2, density=False,label=' ')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,alpha = 0.4, histtype='step',ec = EQs[name]['color'],lw=0.2,density=False,label=' ')
            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label=' ')
        # n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = 'black',lw=0.25,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'black',lw=0.2, density=False,label=r'$\Sigma_{0.2}$')
        # n, bins, patches = ax.hist(av10_stress, nbins,histtype = 'step',fill = False, alpha = 0.4, ec = 'red',lw=1,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,histtype = 'step',fill=False,alpha = 0.4,ec = 'black',lw=1, density=False,label=r'$\Sigma_{0.2}$')
        
        ax.text(EQs[name]['loc'][0],EQs[name]['loc'][1],name,fontsize = 9.5)
        # add a 'best fit' line
        ax.set_xlabel('Static Stress Drop (MPa)')
        
        ax.set_ylabel('Count')
        ax.legend(handlelength=1,handletextpad = 0.1,columnspacing =0.1,labelspacing=0.175,ncols=5,loc="upper right", fontsize=8)
        
    '''
       
        # ax.legend(title=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,d\Sigma}{|\int \mathbf{\Delta u}\,d\Sigma|}$')
    print('potency')
    print(values)
    print(stds)
    #fig.suptitle('Average Stress Drop Across Events',y= 0.96,fontweight='bold',fontsize=14.25)
    plt.savefig(f'{folder}/all_average_strain_histograms_label.png')
    plt.show()


plot_all_hist_potency(EQs,samples = 5000)



EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges':[12,18.5],'rake':90,'RotAngle':None,'color':'blue','loc':(13.0,825)},
       'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges':[4,8],'rake':107,'RotAngle':None,'color':'red','loc':(14.1,825)},
       'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges':[5,13],'rake':90,'RotAngle':None,'color':'orange','loc':(15.2,825)},
       'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges':[5,17],'rake':90,'RotAngle':None,'color':'seagreen','loc':(16.35,825)},
       'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges':[4,10],'rake':None,'RotAngle':360-99,'color':'brown','loc':(17.2,825)}}



def plot_all_hist_stress(EQs,samples= 1000,thr = [0.1,0.2],folder = 'stress_samples_5000'):
    # av_stress =  ave_stress(name,parameters,samples,shape)
    # av10_stress = ave_stress(name,parameters,samples,shape,th = thrs[0])
    # av20_stress = ave_stress(name,parameters,samples,shape,th = thrs[1])
    
    values = {'all':{},'10':{},'20':{}}
    stds =  {'all':{},'10':{},'20':{}}
    fig, ax = plt.subplots(figsize=(10,5),dpi=900)
    for name in list(EQs.keys()):
        nbins = np.arange(EQs[name]['edges'][0],EQs[name]['edges'][1]+0.1,0.1)
        av_stress =  ave_stress(name,EQs[name]['n'],EQs[name]['shape'],samples,th=0,rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder)
        av10_stress = ave_stress(name,EQs[name]['n'],EQs[name]['shape'],samples,th = thr[0],rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder)
        av20_stress = ave_stress(name,EQs[name]['n'],EQs[name]['shape'],samples,th = thr[1],rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder)
        
        values['all'][name] = np.nanmean(np.array(av_stress))
        values['10'][name] = np.nanmean(np.array(av10_stress))
        values['20'][name] = np.nanmean(np.array(av20_stress))
        
        stds['all'][name] = np.nanstd(np.array(av_stress))
        stds['10'][name] = np.nanstd(np.array(av10_stress))
        stds['20'][name] = np.nanstd(np.array(av20_stress))
        
  
        if name == 'Pedernales':

            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label='$\overline{\Delta\sigma}_{E}$')
            #n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label = '$\overline{\Delta\sigma}_{20\%}$')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label = '$\overline{\Delta\sigma}_{10\%}$')
            
            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\Delta\sigma}_{E}$')
            n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = EQs[name]['color'],histtype='stepfilled',lw=0.2,density=False,label = '$\overline{\Delta\sigma}_{10\%}$')
            n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label = '$\overline{\Delta\sigma}_{20\%}$')
        else:
            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label=' ')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label=' ')
            #n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label=' ')


            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label=' ')
            n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, histtype='stepfilled',ec = EQs[name]['color'],lw=0.2,density=False,label=' ')
            n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label=' ')

            
        # n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = 'black',lw=0.25,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'black',lw=0.2, density=False,label=r'$\Sigma_{0.2}$')
        # n, bins, patches = ax.hist(av10_stress, nbins,histtype = 'step',fill = False, alpha = 0.4, ec = 'red',lw=1,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,histtype = 'step',fill=False,alpha = 0.4,ec = 'black',lw=1, density=False,label=r'$\Sigma_{0.2}$')
        
        ax.text(EQs[name]['loc'][0],EQs[name]['loc'][1],name,fontsize = 7,fontweight='bold')
        # add a 'best fit' line
        ax.set_xlabel('Static Stress Drop (MPa)')
        
        ax.set_ylabel('Count')
        ax.legend(frameon=False,handlelength=2.750,handleheight=1.25,handletextpad = 0.2,columnspacing =0.5,labelspacing=0.8,ncols=5,loc="upper right", fontsize=10)
        
        # the histogram of the data
        # n, bins, patches = axes[0].hist(av_stress_dip, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \Delta u_{dip}\cdot\Delta \sigma_{dip}\,dS}{\int \Delta u_{dip} \,dS}$')
        # n, bins, patches = axes[1].hist(av_stress_combined, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,dS}{\int \Delta u_{dip} \,dS}$')
        
    '''
        if name == 'Pedernales':

            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label='$\bar{\Delta\sigma}_{E}$')
            n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label = '$\bar{\Delta\sigma}_{20\%}$')
            n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label = '$\bar{\Delta\sigma}_{10\%}$')
            # n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'none',lw=0.2, density=False,label = '$\overline{\Delta\sigma}_{20\%}$')
           # n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = 'none',lw=0.2,density=False,label = '$\overline{\Delta\sigma}_{10\%}$')
           # n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\Delta\sigma}_{E}$')
        else:
            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label=' ')
            n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label=' ')
            n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label=' ')


            #n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'none',lw=0.2, density=False,label=' ')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,alpha = 0.4, histtype='step',ec = EQs[name]['color'],lw=0.2,density=False,label=' ')
            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label=' ')
        # n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = 'black',lw=0.25,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = 'black',lw=0.2, density=False,label=r'$\Sigma_{0.2}$')
        # n, bins, patches = ax.hist(av10_stress, nbins,histtype = 'step',fill = False, alpha = 0.4, ec = 'red',lw=1,density=False,label=r'$\Sigma_{0.1}$')
        # n, bins, patches = ax.hist(av20_stress, nbins,histtype = 'step',fill=False,alpha = 0.4,ec = 'black',lw=1, density=False,label=r'$\Sigma_{0.2}$')
        
        ax.text(EQs[name]['loc'][0],EQs[name]['loc'][1],name,fontsize = 9.5)
        # add a 'best fit' line
        ax.set_xlabel('Static Stress Drop (MPa)')
        
        ax.set_ylabel('Count')
        ax.legend(handlelength=1,handletextpad = 0.1,columnspacing =0.1,labelspacing=0.175,ncols=5,loc="upper right", fontsize=8)
        
    '''
    print('stress')   
        # ax.legend(title=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,d\Sigma}{|\int \mathbf{\Delta u}\,d\Sigma|}$')
    print(values)
    print(stds)
    #fig.suptitle('Average Stress Drop Across Events',y= 0.96,fontweight='bold',fontsize=14.25)
    plt.savefig(f'{folder}/all_average_stress_histograms_label.png')
    plt.show()


plot_all_hist_stress(EQs,samples = 5000)
# fig, axes = plt.subplots(3,1,figsize = (7,10),dpi=400)
# nbins = 40
# # the histogram of the data
# n, bins, patches = axes[0].hist(av_stress_dip, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \Delta u_{dip}\cdot\Delta \sigma_{dip}\,dS}{\int \Delta u_{dip} \,dS}$')
# n, bins, patches = axes[1].hist(av_stress_combined, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,dS}{\int \Delta u_{dip} \,dS}$')
# n, bins, patches = axes[2].hist(av_stress_amplitude, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,dS}{|\int \Delta u \,dS|}$')

# axes[0].set_title(f'Kinematic {name} for {samples} samples')
# # add a 'best fit' line
# axes[0].set_xlabel('stress drop (MPa)')

# axes[0].set_ylabel('Count')
# axes[1].set_xlabel('stress drop (MPa)')
# axes[1].set_ylabel('Count')
# axes[0].legend()
# axes[1].legend()
# axes[2].legend()
# plt.show()
# plt.savefig(f'histogram/{name}_av_stress_histograms_{samples}.png')
