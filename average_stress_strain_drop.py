# -*- coding: utf-8 -*-
"""
Created on Sun Jan 19 11:59:33 2025

@author: joanv
"""



import h5py 
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import matplotlib.patches as m_patches
import os 
import re


def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)
    
def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()

    
working_dir = os.getcwd()




def build_model(names,geoms,patches,arrow_sizes,nparams,rakes):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['geom'] =  geoms[i] 
        model[name]['patch'] = patches[i]
        model[name]['arrow_size'] = arrow_sizes[i]
        model[name]['nparam'] = nparams[i]
        model[name]['rake'] = rakes[i] 
        model[name]['z_offset'] = z_offsets[i]
        model[name]['nramp'] = nramps[i]
    return model


def get_mu(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VP'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        return mu


def layered_model(name):
    f = open(f'1d_models/model_{name}', 'r')
    content = f.readlines()
    earth_model = content[12:]
    keys = re.findall('\w+',content[11],flags=re.DOTALL)
    fmt = re.compile('\d+.\d+')
    layers = [list(map(float,re.findall(fmt,m))) for m in earth_model]
    layers = np.array(layers).T
    layers0 = np.zeros(layers.shape)
    layers0[0,:] = layers[0,:]
    layers0[1,:] = layers[2,:]
    layers0[2,:] = layers[1,:]
    layers0[3,:] = layers[3,:]
    layers0[4,:] = layers[5,:]
    layers0[5,:] = layers[4,:]
    
    model_dict = dict(zip(keys,layers0))
    model_df = pd.DataFrame(model_dict)
    
    return model_df
    

names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,90,90,0,107]
z_offsets = [2.8,0,0,0,0]
nramps = [0,3,0,9,0]

models = build_model(names,geoms,patches,arrow_sizes,nparams,rakes)

km = 1000


def observables(name,parameters,shape,samples,dir_folder= 'histograms',rake=90, RotAngle=None,th=None):


    df = pd.read_csv(f'INPUT/{name}/model/kinematic/{samples}_samples/mean/{name}_mean_kinematic_model.csv')


    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]

    shape = (nrows,ncols)
    
    # layered model
    model_df = layered_model(name)
    
    
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0).flatten() 
    MU =  np.ones_like(DEPTH)
    for j in range(len(MU)):
        mu = get_mu(dict(model_df),DEPTH[j])
        MU[j] = mu
    
        
    # binary data with sampled kinematic models     
    bin_data = np.fromfile(f'INPUT/{name}/model/kinematic/{samples}_samples/bin_data/{name}_kinematic_n_{samples}.dat','float').reshape((parameters,int(samples)))
    
    # HDF5 data with stress change ensemble 
    data = h5py.File(f'stress_change_{samples}_samples/{name}_Stress change_nsamples_{samples}.h5')


    npatches = data['Stress change'].shape[1]//3
    stress = data['Stress change']
    
    # only used by Pedernales   
    geometry_d = pd.read_csv(f'INPUT/Pedernales/model/kinematic/{samples}_samples/mean/Pedernales_mean_kinematic_model.csv')


    # allocate average stress drop arrays
    av_stress_amplitude = np.zeros(samples)
    av_stress_amplitude_layered = np.zeros(samples)
    av_potency_amplitude = np.zeros(samples)
    
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
        
        
        # layered 
        dip_stress_layered = -stress[i,2*npatches:]*MU/(30e9)
        strike_stress_layered = -stress[i,npatches:2*npatches]*MU/(30e9) 
        dip_stress_layered = dip_stress_layered[ind_Uth]
        strike_stress_layered = strike_stress_layered[ind_Uth]
        
        dot_prod_dip_layered = np.dot(Udip,dip_stress_layered)
        dot_prod_strike_layered = np.dot(Ustrike,strike_stress_layered)
        
        
        
        # halfspace
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
        
        
        # potency 
        dip_stress = -stress[i,2*npatches:] 
        strike_stress = -stress[i,npatches:2*npatches] 
        dip_potency = dip_stress[ind_Uth]/(2*30e3)
        strike_potency = strike_stress[ind_Uth]/(2*30e3)
        dot_prod_dip_potency = np.dot(Udip,dip_potency)
        dot_prod_strike_potency = np.dot(Ustrike,strike_potency)
        

        
        try:
  
            stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
            potency_amplitude = (dot_prod_dip_potency + dot_prod_strike_potency)/denominator_amplitude
            stress_drop_amplitude_layered = (dot_prod_dip_layered + dot_prod_strike_layered)/denominator_amplitude
      

                
            # for spatial averages: uncomment below 
            #stress_drop_amplitude = np.sqrt((np.dot(dip_stress,np.ones_like(dip_stress)))**2 + (np.dot(strike_stress,np.ones_like(strike_stress)))**2)/len(Ustrike)
    
   
            av_stress_amplitude[i] = stress_drop_amplitude
            av_stress_amplitude_layered[i] = stress_drop_amplitude_layered
            av_potency_amplitude[i] = potency_amplitude
            
        except:
            continue
 
    return av_stress_amplitude, av_stress_amplitude_layered,av_potency_amplitude




EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges_stress':[4,27],'edges_potency':[1,18.5e3/(2*30)],'rake':90,'RotAngle':None,'color':'blue','loc_stress':(15.0,590),'loc_potency':(255,740)},
       'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges_stress':[0,9],'edges_potency':[0,8e3/(2*30)],'rake':107,'RotAngle':None,'color':'red','loc_stress':(16.9,590),'loc_potency':(285,740)},
       'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges_stress':[0,18],'edges_potency':[0,13e3/(2*30)],'rake':90,'RotAngle':None,'color':'orange','loc_stress':(18.6,590),'loc_potency':(316,740)},
       'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges_stress':[0,18],'edges_potency':[0,17e3/(2*30)],'rake':90,'RotAngle':None,'color':'seagreen','loc_stress':(20.6,590),'loc_potency':(347,740)},
       'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges_stress':[0,13],'edges_potency':[0,10e3/(2*30)],'rake':None,'RotAngle':360-99,'color':'brown','loc_stress':(22.3,590),'loc_potency':(375,740)}}



def plot_all_hist_observables(EQs,samples= 1000,thr = [0.1,0.2],folder_out = 'stress_change_5000_samples'):

    stress = {'all':{},'10':{},'20':{}}
    std_stress =  {'all':{},'10':{},'20':{}}
    
    stress_layered = {'all':{},'10':{},'20':{}}
    std_stress_layered =  {'all':{},'10':{},'20':{}}
    
    potency = {'all':{},'10':{},'20':{}}
    std_potency =  {'all':{},'10':{},'20':{}}
    
    fig, axes = plt.subplots(2,1,figsize=(8,8))
    
    
    for name in list(EQs.keys()):
        nbins_potency = np.arange(EQs[name]['edges_potency'][0],EQs[name]['edges_potency'][1]+0.1e3/(2*30),0.1e3/(2*30))
        nbins_stress = np.arange(EQs[name]['edges_stress'][0],EQs[name]['edges_stress'][1]+0.1,0.1)
        
        av_stress, av_stress_layered, av_potency =  observables(name,EQs[name]['n'],EQs[name]['shape'],samples,th=0,rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder_out)
        av10_stress, av10_stress_layered, av10_potency= observables(name,EQs[name]['n'],EQs[name]['shape'],samples,th = thr[0],rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder_out)
        av20_stress, av20_stress_layered, av20_potency = observables(name,EQs[name]['n'],EQs[name]['shape'],samples,th = thr[1],rake=EQs[name]['rake'],RotAngle=EQs[name]['RotAngle'],dir_folder=folder_out)
        
        av_stress,av_stress_layered, av_potency = np.array(av_stress),np.array(av_stress_layered),np.array(av_potency)*1e6
        av10_stress,av10_stress_layered, av10_potency = np.array(av10_stress),np.array(av10_stress_layered),np.array(av10_potency)*1e6
        av20_stress, av20_stress_layered, av20_potency = np.array(av20_stress),np.array(av20_stress_layered),np.array(av20_potency)*1e6
        
        
        
        stress['all'][name] = np.nanmean(np.array(av_stress))
        stress['10'][name] = np.nanmean(np.array(av10_stress))
        stress['20'][name] = np.nanmean(np.array(av20_stress))
        
        std_stress['all'][name] = np.nanstd(np.array(av_stress))
        std_stress['10'][name] = np.nanstd(np.array(av10_stress))
        std_stress['20'][name] = np.nanstd(np.array(av20_stress))
        
        
        stress_layered['all'][name] = np.nanmean(np.array(av_stress_layered))
        stress_layered['10'][name] = np.nanmean(np.array(av10_stress_layered))
        stress_layered['20'][name] = np.nanmean(np.array(av20_stress_layered))
        
        std_stress_layered['all'][name] = np.nanstd(np.array(av_stress_layered))
        std_stress_layered['10'][name] = np.nanstd(np.array(av10_stress_layered))
        std_stress_layered['20'][name] = np.nanstd(np.array(av20_stress_layered))
        
        

        
        potency['all'][name] = np.nanmean(np.array(av_potency))
        potency['10'][name] = np.nanmean(np.array(av10_potency))
        potency['20'][name] = np.nanmean(np.array(av20_potency))
        
        std_potency['all'][name] = np.nanstd(np.array(av_potency))
        std_potency['10'][name] = np.nanstd(np.array(av10_potency))
        std_potency['20'][name] = np.nanstd(np.array(av20_potency))
        
        
        df_stress = pd.DataFrame(stress)
        df_stress_layered = pd.DataFrame(stress_layered)
        df_potency = pd.DataFrame(potency)

        if name == 'Pedernales':

            n, bins, patches = axes[0].hist(av_stress_layered, nbins_stress,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\Delta \sigma}_{E}$')
            n, bins, patches = axes[0].hist(av10_stress_layered, nbins_stress,color=EQs[name]['color'],alpha = 0.4, ec = EQs[name]['color'],histtype='stepfilled',lw=0.2,density=False,label = '$\overline{\Delta \sigma}_{E,{10\%}}$')
            n, bins, patches = axes[0].hist(av20_stress_layered, nbins_stress,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label = '$\overline{\Delta \sigma}_{E,{20\%}}$')
            
            n, bins, patches = axes[1].hist(av_potency, nbins_potency,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\epsilon}_{E}$')
            n, bins, patches = axes[1].hist(av10_potency, nbins_potency,color=EQs[name]['color'],alpha = 0.4, ec = EQs[name]['color'],histtype='stepfilled',lw=0.2,density=False,label = '$\overline{\epsilon}_{E,{10\%}}$')
            n, bins, patches = axes[1].hist(av20_potency, nbins_potency,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label = '$\overline{\epsilon}_{E,{20\%}}$')
        
        else:

            n, bins, patches = axes[0].hist(av_stress_layered, nbins_stress,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label=' ')
            n, bins, patches = axes[0].hist(av10_stress_layered, nbins_stress,color=EQs[name]['color'],alpha = 0.4, histtype='stepfilled',ec = EQs[name]['color'],lw=0.2,density=False,label=' ')
            n, bins, patches = axes[0].hist(av20_stress_layered, nbins_stress,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label=' ')

            n, bins, patches = axes[1].hist(av_potency, nbins_potency,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label=' ')
            n, bins, patches = axes[1].hist(av10_potency, nbins_potency,color=EQs[name]['color'],alpha = 0.4, histtype='stepfilled',ec = EQs[name]['color'],lw=0.2,density=False,label=' ')
            n, bins, patches = axes[1].hist(av20_potency, nbins_potency,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label=' ')
       
        

        axes[0].text(EQs[name]['loc_stress'][0],EQs[name]['loc_stress'][1],name,fontsize = 7,fontweight='bold')   
        axes[1].text(EQs[name]['loc_potency'][0],EQs[name]['loc_potency'][1],name,fontsize = 7,fontweight='bold')


        axes[0].text(0.010,1.05,'a',transform=axes[0].transAxes,fontsize=16,fontweight='bold')
        axes[1].text(0.010,1.05,'b',transform=axes[1].transAxes,fontsize=16,fontweight='bold')

        # add a 'best fit' line
        
        axes[0].set_xlabel('Static Stress Drop (MPa)')
        axes[1].set_xlabel('Potency Density (microstrain)')
        
        axes[0].set_ylabel('Count')
        axes[1].set_ylabel('Count')
        axes[0].set_xlim([2,26])
        axes[1].set_xlim([2e3/(2*30),27e3/(2*30)])
        axes[0].legend(frameon=False,handlelength=2.750,handleheight=1.25,handletextpad = 0.2,columnspacing =0.5,labelspacing=0.8,ncols=5,loc="upper right", fontsize=9)
        axes[1].legend(frameon=False,handlelength=2.750,handleheight=1.25,handletextpad = 0.2,columnspacing =0.5,labelspacing=0.8,ncols=5,loc="upper right", fontsize=9)


    print('Stress (Homogeneous halfspace)\n:', df_stress)
    print('Potency (Homogeneous halfspace)\n:', df_potency)
    print('Stress (Layered):', df_stress_layered)
    plt.savefig(f'{folder_out}/all_stress_potency_energy_averages_{samples}.pdf')
    
    # plt.show()


plot_all_hist_observables(EQs,samples = 5000)


