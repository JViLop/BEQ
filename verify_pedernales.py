# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 20:46:03 2024

@author: joanv
"""

import h5py 
import numpy as np
import matplotlib.pyplot as plt
import struct
path = 'EQ/Pedernales/model/kinematic/step_053.h5'

f = h5py.File(path)

cov = np.array(f['Covariance'])
data = np.array(f['Sample Set']).T
print(data.shape)
Np = 80
nramp = 9



zero_id = np.where(~data.any(axis=0))[0]
data = np.delete(data,zero_id,axis=1)


U1 = data[:,:Np]
U2 = data[:,Np:2*Np]
U = np.sqrt(U1**2 + U2**2)
Tr = data[:,2*Np+nramp:3*Np+nramp]
Vr = data[:,3*Np+nramp:4*Np+nramp]



# U1 = np.delete(U1,zero_id,axis=0)
# U2 = np.delete(U2,zero_id,axis=0)
U = np.sqrt(U1**2 + U2**2)
# Tr = np.delete(Tr,zero_id,axis=0)
# Vr = np.delete(Vr,zero_id,axis=0)
counts, bins = np.histogram(Tr,bins=100,density=False)
plt.hist(bins[:-1], bins,density=False,weights=counts)
plt.show()
counts, bins = np.histogram(Vr,bins=100,density=False)
plt.hist(bins[:-1], bins,density=False,weights=counts)
plt.show()

fig,axes = plt.subplots(2,1,dpi=500)
counts, bins = np.histogram(np.mean(Vr,axis=0),bins=20,density=False)
axes[0].hist(bins[:-1], bins,density=False,weights=counts)
axes[0].set_title('Pedernales Vr')
counts, bins = np.histogram(np.mean(Tr,axis=0),bins=20,density=False)
axes[1].hist(bins[:-1], bins,density=False,weights=counts)
axes[1].set_title('Pedernales Tr')
plt.tight_layout()
plt.show()


nsamples = 100
data2 = np.fromfile(f'INPUT/Pedernales/model/kinematic/{nsamples}_samples/bin_data/Pedernales_kinematic_n_{nsamples}.dat')
data2 = data2.reshape((331,nsamples))

U1_2 = data2[:Np,:]
U2_2 = data2[Np:2*Np,:]
U_2 = np.sqrt(U1_2**2 + U2_2**2)
Tr_2 = data2[2*Np+nramp:3*Np+nramp,:]
Vr_2 = data2[3*Np+nramp:4*Np+nramp,:]

fig,axes = plt.subplots(2,1,dpi=500)
counts, bins = np.histogram(Vr_2,bins=100,density=False)
axes[0].hist(bins[:-1], bins,density=False,weights=counts)
axes[0].set_title('Pedernales Vr')
counts, bins = np.histogram(Tr_2,bins=100,density=False)
axes[1].hist(bins[:-1], bins,density=False,weights=counts)
axes[1].set_title('Pedernales Tr')
plt.tight_layout()
plt.show()




rng = np.random.default_rng() 
nsamples = 100
data_sampled = np.array(data).T
sampled_data = rng.choice(data_sampled,nsamples,axis=1)

dat = sampled_data.flatten()

filename  = f'Pedernales_{nsamples}.dat'   
out_file = open(filename,"wb")
s = struct.pack('d'*len(dat), *dat)
out_file.write(s)
out_file.close() 


data3 = np.fromfile(f'Pedernales_{nsamples}.dat')
data3 = data3.reshape((331,nsamples))

U1_3 = data3[:Np,:]
U2_3 = data3[Np:2*Np,:]
U_3 = np.sqrt(U1_3**2 + U2_3**2)
Tr_3 = data3[2*Np+nramp:3*Np+nramp,:]
Vr_3 = data3[3*Np+nramp:4*Np+nramp,:]
fig,axes = plt.subplots(2,1,dpi=500)
counts, bins = np.histogram(Vr_3,bins=100,density=False)
axes[0].hist(bins[:-1], bins,density=False,weights=counts)
axes[0].set_title('Pedernales Vr')
counts, bins = np.histogram(Tr_3,bins=100,density=False)
axes[1].hist(bins[:-1], bins,density=False,weights=counts)
axes[1].set_title('Pedernales Tr')
plt.tight_layout()
plt.show()



### Compare with Iquique ###

   
f2 = h5py.File('EQ/Iquique/model/kinematic/step_056.h5')

def stacking(d,keys,n):
    if n!=0:
        x = np.hstack((stacking(d,keys,n-1),d[keys[n+1]]))
    else:
        x = np.hstack((d[keys[n]],d[keys[n+1]]))
    return x
    

def hstacking(d):   
    keys = ['strikeslip','dipslip','risetime','rupturevelocity','hypo_as','hypo_dd']
    n = len(keys)
    return stacking(d,keys,n-2)




data4 = np.array(f2['Sample Set'])

print(data4.shape)
Np = 11*12
nramp = 3
U1 = data4[:,:Np]
U2 = data4[:,Np:2*Np]
U = np.sqrt(U1**2 + U2**2)
Tr = data4[:,2*Np+nramp:3*Np+nramp]
Vr = data4[:,3*Np+nramp:4*Np+nramp]

# counts, bins = np.histogram(Tr,bins=100,density=False)
# plt.hist(bins[:-1], bins,density=False,weights=counts)
# plt.show()
# counts, bins = np.histogram(Vr,bins=100,density=False)
# plt.hist(bins[:-1], bins,density=False,weights=counts)
# plt.show()


fig,axes = plt.subplots(2,1,dpi=500)
counts, bins = np.histogram(Vr,bins=100,density=False)
axes[0].hist(bins[:-1], bins,density=False,weights=counts)
axes[0].set_title('Iquique Vr')
counts, bins = np.histogram(Tr,bins=100,density=False)
axes[1].hist(bins[:-1], bins,density=False,weights=counts)
axes[1].set_title('Iquique Tr')
plt.tight_layout()
plt.show()




path_sampled = 'INPUT/Iquique/model/kinematic/all_samples/bin_data/Iquique_kinematic_n_all.dat'
data5 = np.fromfile(path_sampled)


data5 = data5.reshape((533,72000))

U1_5 = data5[:Np,:]
U2_5 = data5[Np:2*Np,:]
U_5 = np.sqrt(U1_5**2 + U2_5**2)
Tr_5 = data5[2*Np+nramp:3*Np+nramp,:]
Vr_5 = data5[3*Np+nramp:4*Np+nramp,:]
fig,axes = plt.subplots(2,1,dpi=500)
counts, bins = np.histogram(Vr_5,bins=100,density=False)
axes[0].hist(bins[:-1], bins,density=False,weights=counts)
axes[0].set_title('Iquique Vr')
counts, bins = np.histogram(Tr_5,bins=100,density=False)
axes[1].hist(bins[:-1], bins,density=False,weights=counts)
axes[1].set_title('Iquique Tr')
plt.tight_layout()
plt.show()





name = 'Pedernales'
names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,331,650]
rakes = [90,0,0,0,0]
z_offsets = [2.8,0,0,0,0]
nramps = [0,3,0,9,0]
def build_model(names,geoms,patches,arrow_sizes,nparams,rakes,nramps):
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

models = build_model(names,geoms,patches,arrow_sizes,nparams,rakes,nramps)

path_sampled = f'INPUT/{name}/model/kinematic/all_samples/bin_data/{name}_kinematic_n_all.dat'
data = np.fromfile(path_sampled)

nparam = models[name]['nparam']
nramp = models[name]['nramp']

data = data.reshape((nparam,len(data)//nparam))
if name=='Pedernales':
    zero_id = np.where(~data.any(axis=0))[0]
    
    data = np.delete(data,zero_id,axis=1)
    zero_id_tr = np.where(data[2*Np+nramp:3*Np+nramp,:]==0)[1]
    data = np.delete(data,zero_id_tr,axis=1)
    zero_id_vr = np.where(data[3*Np+nramp:4*Np+nramp,:]==0)[1]
    data = np.delete(data,zero_id_vr,axis=1)
    zero_id_hyp = np.where(data[4*Np+nramp+1:,:]==0)[1]
    data = np.delete(data,zero_id_hyp,axis=1)
    zero_id_hyp = np.where(data[4*Np+nramp:4*Np+nramp+1,:]==0)[1]
    data = np.delete(data,zero_id_hyp,axis=1)

Np = models[name]['geom'][0]*models[name]['geom'][1]
U1 = data[:Np,:]
U2 = data[Np:2*Np,:]
U = np.sqrt(U1**2 + U2**2)
Tr = data[2*Np+nramp:3*Np+nramp,:]
Vr = data[3*Np+nramp:4*Np+nramp,:]
hyp_strike = data[4*Np+nramp:4*Np+nramp+1,:]
hyp_dip = data[4*Np+nramp+1:,:]

plt.rc('axes', titlesize=10)
fig,axes = plt.subplots(3,2,dpi=500)

counts, bins = np.histogram(U1,bins=100,density=False)
axes[0][0].hist(bins[:-1], bins,density=False,weights=counts)
axes[0][0].set_title('U1 (m)')
axes[0][0].text(np.min(U1),-max(counts)/5,'%s'%(round(np.min(U1),2)),fontsize=5)
axes[0][0].text(np.max(U1),-max(counts)/5,'%s'%(round(np.max(U1),2)),fontsize=5)
axes[0][0].axvline(np.min(U1),linestyle='dashed',color='k',linewidth=1)
axes[0][0].axvline(np.max(U1),linestyle='dashed',color='k',linewidth=1)


counts, bins = np.histogram(U2,bins=100,density=False)
axes[0][1].hist(bins[:-1], bins,density=False,weights=counts)
axes[0][1].set_title('U2 (m)')
axes[0][1].text(np.min(U2),-max(counts)/5,'%s'%(round(np.min(U2),2)),fontsize=5)
axes[0][1].text(np.max(U2),-max(counts)/5,'%s'%(round(np.max(U2),2)),fontsize=5)
axes[0][1].axvline(np.min(U2),linestyle='dashed',color='k',linewidth=1)
axes[0][1].axvline(np.max(U2),linestyle='dashed',color='k',linewidth=1)

counts, bins = np.histogram(Tr,bins=100,density=False)
axes[1][0].hist(bins[:-1], bins,density=False,weights=counts)
axes[1][0].set_title('Tr (s)')
axes[1][0].text(np.min(Tr),-max(counts)/5,'%s'%(round(np.min(Tr),2)),fontsize=5)
axes[1][0].text(np.max(Tr),-max(counts)/5,'%s'%(round(np.max(Tr),2)),fontsize=5)
axes[1][0].axvline(np.min(Tr),linestyle='dashed',color='k',linewidth=1)
axes[1][0].axvline(np.max(Tr),linestyle='dashed',color='k',linewidth=1)

counts, bins = np.histogram(Vr,bins=100,density=False)
axes[1][1].hist(bins[:-1], bins,density=False,weights=counts)
axes[1][1].set_title('Vr (km/s)')
axes[1][1].text(np.min(Vr),-max(counts)/5,'%s'%(round(np.min(Vr),2)),fontsize=5)
axes[1][1].text(np.max(Vr),-max(counts)/5,'%s'%(round(np.max(Vr),2)),fontsize=5)
axes[1][1].axvline(np.min(Vr),linestyle='dashed',color='k',linewidth=1)
axes[1][1].axvline(np.max(Vr),linestyle='dashed',color='k',linewidth=1)

counts, bins = np.histogram(hyp_strike,bins=100,density=False)
axes[2][0].hist(bins[:-1], bins,density=False,weights=counts)
axes[2][0].set_title('Hypocenter (along-strike) (km)')
axes[2][0].text(np.min(hyp_strike),-max(counts)/5,'%s'%(round(np.min(hyp_strike),2)),fontsize=5)
axes[2][0].text(np.max(hyp_strike),-max(counts)/5,'%s'%(round(np.max(hyp_strike),2)),fontsize=5)
axes[2][0].axvline(np.min(hyp_strike),linestyle='dashed',color='k',linewidth=1)
axes[2][0].axvline(np.max(hyp_strike),linestyle='dashed',color='k',linewidth=1)

counts, bins = np.histogram(hyp_dip,bins=100,density=False)
axes[2][1].hist(bins[:-1], bins,density=False,weights=counts)
axes[2][1].set_title('Hypocenter (along-dip) (km)')
axes[2][1].text(np.min(hyp_dip),-max(counts)/5,'%s'%(round(np.min(hyp_dip),2)),fontsize=5)
axes[2][1].text(np.max(hyp_dip),-max(counts)/5,'%s'%(round(np.max(hyp_dip),2)),fontsize=5)
axes[2][1].axvline(np.min(hyp_dip),linestyle='dashed',color='k',linewidth=1)
axes[2][1].axvline(np.max(hyp_dip),linestyle='dashed',color='k',linewidth=1)

fig.suptitle(f'{name}',fontweight='bold')
plt.tight_layout()
plt.show()





'''
    else:
        self.sampled_data = data
        dat = self.sampled_data.flatten()
        try: 
            os.makedirs(self.output_bin_data_dir)    
        except FileExistsError:
           pass   
        
        filename  = f"{self.name}_{self.model_type}_n_{self.nsamples}.dat"    
        output_file = os.path.join(self.output_bin_data_dir,filename)
        out_file = open(output_file,"wb")
        s = struct.pack('d'*len(dat), *dat)
        out_file.write(s)
        out_file.close()
'''

# Tr_2 = data2[:,2*Np+nramp:3*Np+nramp].T
# counts, bins = np.histogram(Tr_2,bins=100,density=False)
# plt.hist(bins[:-1], bins,density=False,weights=counts,color='red')
# plt.show()


# Vr_2 = data2[:,3*Np+nramp:4*Np+nramp].T
# counts, bins = np.histogram(Vr_2,bins=100,density=False)
# plt.hist(bins[:-1], bins,density=False,weights=counts,color='blue')
# plt.show()
