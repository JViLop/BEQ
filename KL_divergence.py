# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:22:33 2024

@author: joanv
"""

import numpy as np
import multiprocessing as mp
from scipy.spatial import KDTree
#from sklearn.neighbors import NearestNeighbors


def knn_distance(point, sample, k):
    """Euclidean distance from `point` to it's `k`-Nearest
    Neighbour in `sample`

    This function works for points in arbitrary dimensional spaces.
    """
    # Compute all euclidean distances
    norms = np.linalg.norm(sample - point, axis=1)
    # Return the k-th nearest
    return np.sort(norms)[k]


def verify_sample_shapes(s1, s2, k):
    # Expects [N, D]
    assert len(s1.shape) == len(s2.shape) == 2
    # Check dimensionality of sample is identical
    assert s1.shape[1] == s2.shape[1]


def naive_estimator(s1, s2, k=1):
    """KL-Divergence estimator using brute-force (numpy) k-NN
    s1: (N_1,D) Sample drawn from distribution P
    s2: (N_2,D) Sample drawn from distribution Q
    k: Number of neighbours considered (default 1)
    return: estimated D(P|Q)
    """
    verify_sample_shapes(s1, s2, k)

    n, m = len(s1), len(s2)
    D = np.log(m / (n - 1))
    d = float(s1.shape[1])

    for p1 in s1:
        nu = knn_distance(p1, s2, k - 1)  # -1 because 'p1' is not in 's2'
        rho = knn_distance(p1, s1, k)
        if rho==0:
             continue
    
        D += (d / n) * np.log(nu / rho)
    
    return D


samples = 500
step = '52'
file_dir = f'INPUT/Tohoku/model/kinematic/step_{step}/{samples}_samples/bin_data/Tohoku_kinematic_n_{samples}.dat'
nparameters = 866
data = np.fromfile(file_dir,'float').reshape((nparameters,int(samples)))
Q = data.T

step = '37'
file_dir = f'INPUT/Tohoku/model/kinematic/step_{step}/{samples}_samples/bin_data/Tohoku_kinematic_n_{samples}.dat'
nparameters = 866
data = np.fromfile(file_dir,'float').reshape((nparameters,int(samples)))
P = data.T
n = len(Q)

D = naive_estimator(Q,P,k=1)
print(D)

step = '44'
file_dir = f'INPUT/Tohoku/model/kinematic/step_{step}/{samples}_samples/bin_data/Tohoku_kinematic_n_{samples}.dat'
nparameters = 866
data = np.fromfile(file_dir,'float').reshape((nparameters,int(samples)))
P = data.T


D = naive_estimator(Q,P,k=1)
print(D)


# def kl_estimator_parallel(p1,s1,s2,k,q):
#     """KL-Divergence estimator using brute-force (numpy) k-NN
#     s1: (N_1,D) Sample drawn from distribution P
#     s2: (N_2,D) Sample drawn from distribution Q
#     k: Number of neighbours considered (default 1)
#     return: estimated D(P|Q)
#     """
    
#     n= len(s1)
#     d = float(s2.shape[1])

    
#     nu = knn_distance(p1, s2, k - 1)  # -1 because 'p1' is not in 's2'
#     rho = knn_distance(p1, s1, k)
#     # if rho==0:
#     #     D = 0.0
#     # else:
         
#     D = (d / n) * np.log(nu / rho)
#     D = 0
#     q.put(D)

# if __name__ == "__main__":
    
#     qout = mp.Queue()

#     inputs = [(q, Q, P, 1, qout) for q in Q]
#     processes = [mp.Process(target=kl_estimator_parallel, args=inp)
#                   for inp in inputs]
    
#     for p in processes:
#         p.start()
    
#     for p in processes:
#         p.join()
    
    # result = [qout.get() for p in processes]
    # print(result)
    # result = np.array(result)
    
    # D = np.log(len(P) / (len(Q) - 1))
    # D = D + np.sum(result)
    # print(D)

