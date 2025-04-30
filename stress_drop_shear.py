# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 22:44:21 2024

@author: joanv
"""


import h5py 
import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd
import matplotlib.patches as m_patches
import os 
import re

name = 'Tohoku'
parameters = 866
samples  = 100
shape = (9,24)
th = 0.2

working_dir = os.getcwd()


def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)
    
def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()
    
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


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,90,90,0,107]
z_offsets = [2.8,0,0,0,0]
nramps = [0,3,0,9,0]
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

models = build_model(names,geoms,patches,arrow_sizes,nparams,rakes)

km = 1000


model_type = 'kinematic'
nsamples = 100

    
def ave_stress(name,parameters,shape,samples,dir_folder= 'histograms',rake=90, RotAngle=None,th=None):
    z_offset = models[name]['z_offset']

    earth_model_dir = os.path.join(working_dir,f'1d_models/model_{name}')
    f = open(earth_model_dir, "r")
    content = f.readlines()
    df = pd.read_csv(f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')

    #nlayers =int(sys.argv[4])

    nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
    nparam = models[name]['nparam']

    shape = (nrows,ncols)
    patch = models[name]['patch']
    nramp = models[name]['nramp']
    earth_model = content[12:]
    keys = re.findall('\w+',content[11],flags=re.DOTALL)
    fmt = re.compile('\d+.\d+')
    layers = [list(map(float,re.findall(fmt,m))) for m in earth_model]
    #print(layers)
    layers = np.array(layers).T
    layers0 = np.zeros(layers.shape)
    layers0[0,:] = layers[0,:]
    layers0[1,:] = layers[2,:]
    layers0[2,:] = layers[1,:]
    layers0[3,:] = layers[3,:]
    layers0[4,:] = layers[5,:]
    layers0[5,:] = layers[4,:]
    model_dict = dict(zip(keys,layers0))
    #print(layers0)
    
    model_df = pd.DataFrame(model_dict)
    
    
    DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0).flatten() 
    MU =  np.ones_like(DEPTH)

    for j in range(len(MU)):
        mu = get_mu(dict(model_df),DEPTH[j])
        MU[j] = mu

 # convert to GPA
    bin_data = np.fromfile(f'{dir_folder}/{name}_kinematic_n_{samples}.dat','float').reshape((parameters,int(samples)))

    # data = h5py.File('Tohoku_Stress change_nsamples_100.h5')

    # npatches = data['Stress change'].shape[1]//3
    # stress = data['Stress change']

    data = h5py.File(f'{dir_folder}/{name}_Stress_change_nsamples_{samples}.h5')
    #data = h5py.File(f'{dir_folder}/{name}_stress_nsamples_{samples}.h5')


    npatches = data['Stress change'].shape[1]//3
    stress = data['Stress change']
    
    # just valid for Pedernales 
    
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
            strike = df['strike']*((np.pi)/180)
            dip = df['dip']*((np.pi)/180)
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
        
        dip_stress = -stress[i,2*npatches:]*MU/(30e9)
        strike_stress = -stress[i,npatches:2*npatches]*MU/(30e9) 
        dip_stress = dip_stress[ind_Uth]
        strike_stress = strike_stress[ind_Uth]
        
        dot_prod_dip = np.dot(Udip,dip_stress)
        dot_prod_strike = np.dot(Ustrike,strike_stress)
        # denominator = np.sum(Udip)
     
        
        denominator_amplitude = np.sqrt((np.sum(Udip))**2 + (np.sum(Ustrike))**2)
        
        
        # stress_drop_dip = dot_prod_dip/denominator
        # stress_drop_combined = (dot_prod_dip + dot_prod_strike)/denominator
        try:
            # corrected for spatial average
            if th==0:
                stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
                
            
            else:
                stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
                #stress_drop_amplitude = np.abs(np.dot(dip_stress,np.ones_like(dip_stress)))/len(Ustrike)
                
                # uncomment the one below to compute spatial average
                #stress_drop_amplitude = np.sqrt((np.dot(dip_stress,np.ones_like(dip_stress)))**2 + (np.dot(strike_stress,np.ones_like(strike_stress)))**2)/len(Ustrike)
            # before
            #stress_drop_amplitude = (dot_prod_dip + dot_prod_strike)/denominator_amplitude
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
    
    
    fig, ax = plt.subplots(figsize=(6,4),dpi=1200)
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


'''
EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges':[4,25],'rake':90,'RotAngle':None,'color':'blue','loc':(15.75,1325)},
       'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges':[0,9],'rake':107,'RotAngle':None,'color':'red','loc':(17.4,1325)},
       'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges':[0,18],'rake':90,'RotAngle':None,'color':'orange','loc':(19.25,1325)},
       'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges':[0,18],'rake':90,'RotAngle':None,'color':'seagreen','loc':(21.25,1325)},
       'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges':[0,13],'rake':None,'RotAngle':360-99,'color':'brown','loc':(22.75,1325)}}
'''
EQs = {'Tohoku':{'Mw':9.0,'n':866,'shape':(9,24),'edges':[4,25],'rake':90,'RotAngle':None,'color':'blue','loc':(16.55,660)},
       'Gorkha':{'Mw':7.8,'n':650,'shape':(9,18),'edges':[0,9],'rake':107,'RotAngle':None,'color':'red','loc':(18.3,660)},
       'Iquique':{'Mw':8.1,'n':533,'shape':(11,12),'edges':[0,18],'rake':90,'RotAngle':None,'color':'orange','loc':(19.9,660)},
       'Illapel':{'Mw':8.3,'n':682,'shape':(10,17),'edges':[0,18],'rake':90,'RotAngle':None,'color':'seagreen','loc':(21.6,660)},
       'Pedernales':{'Mw':7.8,'n':331,'shape':(8,10),'edges':[0,13],'rake':None,'RotAngle':360-99,'color':'brown','loc':(23.0,660)}}
def plot_all_hist(EQs,samples= 1000,thr = [0.1,0.2],folder = 'stress_samples_5000'):
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
        
  
        
        
        # the histogram of the data
        # n, bins, patches = axes[0].hist(av_stress_dip, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \Delta u_{dip}\cdot\Delta \sigma_{dip}\,dS}{\int \Delta u_{dip} \,dS}$')
        # n, bins, patches = axes[1].hist(av_stress_combined, nbins, density=False,label=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,dS}{\int \Delta u_{dip} \,dS}$')
        
    
        if name == 'Pedernales':

            #n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=1.5,density=False,label='$\overline{\Delta\sigma}_{E}$')
            #n, bins, patches = ax.hist(av20_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=0.5, density=False,label = '$\overline{\Delta\sigma}_{20\%}$')
            #n, bins, patches = ax.hist(av10_stress, nbins,color=None,histtype='step',ec = EQs[name]['color'],lw=1,density=False,label = '$\overline{\Delta\sigma}_{10\%}$')
            
            n, bins, patches = ax.hist(av_stress, nbins,alpha = 0.75, color=EQs[name]['color'],histtype='stepfilled',ec = 'black',lw=0.6,density=False,label='$\overline{\Delta\sigma}_{E}$')
            n, bins, patches = ax.hist(av10_stress, nbins,color=EQs[name]['color'],alpha = 0.4, ec = EQs[name]['color'],histtype='stepfilled',lw=0.2,density=False,label = '$\overline{\Delta\sigma}_{E_{10\%}}$')
            n, bins, patches = ax.hist(av20_stress, nbins,color=EQs[name]['color'],alpha = 0.2,ec = EQs[name]['color'],histtype='stepfilled',lw=0.2, density=False,label = '$\overline{\Delta\sigma}_{E_{20\%}}$')
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
        ax.set_xlim([2,26])
        ax.legend(frameon=False,handlelength=2.750,handleheight=1.25,handletextpad = 0.2,columnspacing =0.5,labelspacing=0.8,ncols=5,loc="upper right", fontsize=10)
        
    
        # ax.legend(title=r'$\overline{\Delta \sigma} = \frac{\int \mathbf{\Delta u}\cdot \mathbf{\Delta \sigma}\,d\Sigma}{|\int \mathbf{\Delta u}\,d\Sigma|}$')
    print('layered stress')
    print(values)
    print(stds)
    fig.suptitle('Average Stress Drop Across Events',y= 0.96,fontweight='bold',fontsize=14.25)
    plt.savefig(f'{folder}/all_average_stress_histograms_corrected_depth.png')
    plt.show()


plot_all_hist(EQs,samples = 5000)
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
