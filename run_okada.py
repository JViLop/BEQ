

import sys
import os 
import pandas as pd 


from utils.displacement import Displacement
from utils.curves import Curves
from utils.ensemble import Observable,Ensemble_Displacement,Ensemble_Stress
from utils.stress_change import Stress

pwd = os.getcwd()




def model_dict(names,patches,sizes,nparameters,strikes,ylocs, dts, smaxs,scales, larrows):
    model = dict()
    for i,name in enumerate(names):
        model[name] = dict()
        model[name]['patch'] = patches[i]
        model[name]['size'] = sizes[i]
        model[name]['nparameters'] = nparameters[i]
        model[name]['strike'] = strikes[i]
        model[name]['yloc'] = ylocs[i]
        model[name]['dt'] = dts[i]
        model[name]['smax'] = smaxs[i]
        model[name]['scale'] = scales[i]
        model[name]['larrow'] = larrows[i]

    return model

names = ['Tohoku','Iquique','Illapel','Gorkha','Pedernales']
patches = [29,17,18,10,15]
nparameters = [866,533,682,650,331]
sizes = [(9,24),(11,12),(10,17),(9,18),(8,10)]
strikes = [90]*5
ylocs = [0.70,0.93,0.93,0.93,0.93]
dts = [8,8,8,3,3]
scales = [0.05,0.05,0.05,0.05,0.075]
larrows = [1]*5
smaxs = [35,15,15,15,15]


models = model_dict(names,patches,sizes,nparameters,strikes,ylocs,dts,smaxs,scales,larrows)

names = ['Iquique','Pedernales']
nsamples = 100

for name in names:
    event = Observable(pwd,name,'kinematic',None,False,models[name]['patch'],models[name]['size'],observable='RuptureTime',samples=nsamples,
                         new_station=None,xstation=None,ystation=None,strikeOkada = models[name]['strike'],factor=None)
    event.plot_rupture(dt=models[name]['dt'],yloc=models[name]['yloc'],smax=models[name]['smax'])
    event.save_T0()
    if name == 'Pedernales':
        
        ens_d = Ensemble_Displacement(pwd,name,'kinematic',None, False,models[name]['patch'],models[name]['size'],nparameters =models[name]['nparameters'],RotAngle=360-99,samples  = nsamples)
        ens_stress = Ensemble_Stress(pwd,name,'kinematic',None,False,models[name]['patch'],models[name]['size'],nparameters =models[name]['nparameters'],new_station=True,xstation=ens_d.xsrc,ystation=ens_d.ysrc,RotAngle=360-99,samples = nsamples)
    elif name == 'Gorkha':
        ens_d = Ensemble_Displacement(pwd,name,'kinematic',None, False,models[name]['patch'],models[name]['size'],nparameters =models[name]['nparameters'],rake=107,samples  = nsamples)
        ens_stress = Ensemble_Stress(pwd,name,'kinematic',None,False,models[name]['patch'],models[name]['size'],nparameters =models[name]['nparameters'],rake=107,new_station=True,xstation=ens_d.xsrc,ystation=ens_d.ysrc,samples = nsamples)
    else:
        ens_d = Ensemble_Displacement(pwd,name,'kinematic',None, False,models[name]['patch'],models[name]['size'],nparameters =models[name]['nparameters'],samples  = nsamples)
        ens_stress = Ensemble_Stress(pwd,name,'kinematic',None,False,models[name]['patch'],models[name]['size'],nparameters =models[name]['nparameters'],new_station=True,xstation=ens_d.xsrc,ystation=ens_d.ysrc,samples = nsamples)
    
    d = Displacement(pwd,name,'kinematic',None, False,models[name]['patch'],models[name]['size'],larrow=models[name]['larrow'],scale=models[name]['scale'],GF_on=False,samples=nsamples)
    
    stress = Stress(pwd,name,'kinematic',None, False,models[name]['patch'],models[name]['size'],new_station=True,xstation=d.xsrc,ystation=d.ysrc,samples=nsamples)
    
    Curves(pwd,name,'kinematic',None, False,models[name]['patch'],models[name]['size'],1,2,d.xsrc,d.ysrc,d.xstn,d.ystn,d.xsrc,d.ysrc,samples = nsamples)
    
    
