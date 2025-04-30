import numpy as np
import h5py

N = h5py.File('H5_10/Gorkha_N_nmodels_10.h5','r')['N']

E = h5py.File('H5_10/Gorkha_E_nmodels_10.h5','r')['E']

Z = h5py.File('H5_10/Gorkha_Z_nmodels_10.h5','r')['Z']




mean_E = np.mean(E,axis = 0)
print(max(mean_E[:,-1]))

mean_N = np.mean(N,axis = 0)
print(max(mean_N[:,-1]))

mean_Z = np.mean(Z,axis = 0)
print(max(mean_Z[:,-1]))
