import sys
 
# adding csi functions to the system path

from okadafull import displacement, stress, strain

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


        

dr = 10
km = 1000
# station 
xstn = np.arange(-50,50 + dr,dr)*km
ystn = np.arange(-50,50 + dr,dr)*km

Xstn,Ystn = np.meshgrid(xstn,ystn)
xstn_flat,ystn_flat = Xstn.flatten(),Ystn.flatten()
zstn_flat = np.zeros_like(xstn_flat)


# source

xsrc = np.array([0])
ysrc = np.array([0])
zsrc = np.array([20])*km

width = np.array([2])*km
length = np.array([2])*km
strike_rad = np.array([90])*(np.pi/180)
dip_rad = np.array([90])*(np.pi/180)


ss = np.array([5])
ds = np.array([0])
ts = np.array([0])

nu = 0.25
Displacement = displacement(
                                    xstn_flat,
                                    ystn_flat,
                                    zstn_flat, 
                                    xsrc, 
                                    ysrc, 
                                    zsrc, 
                                    width, 
                                    length, 
                                    strike_rad, 
                                    dip_rad, 
                                    ss, 
                                    ds, 
                                    ts, nu=nu)

Displacement_df = pd.DataFrame(Displacement,columns=['x','y','z'])


