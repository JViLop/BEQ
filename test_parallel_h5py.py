# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 12:42:50 2025

@author: joanv
"""

from mpi4py import MPI
import h5py
import numpy as np

from mpi4py import MPI

def myFun(x):
    return x+2 # simple example, the real one would be complicated

comm = MPI.COMM_WORLD
rank = comm.Get_rank() # get your process ID
data = [] # init the data    

if rank == 0: # The master is the only process that reads the file
    data = np.ones(1000) # something read from file

# Divide the data among processes
data = comm.scatter(data, root=0)

result = []
for item in data:
    result.append(myFun(item))

# Send the results back to the master processes
newData = comm.gather(result,root=0)
