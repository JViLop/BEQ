# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:13:57 2024

@author: joanv
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ncols =24
nrows= 9 

name = 'Tohoku_geometry.csv'
df = pd.read_csv(name)

lon0 = df['lon'].values
lat0 = df['lat'].values
lon =lon0.reshape(nrows,ncols,order='F')
lat = lat0.reshape(nrows,ncols,order='F')

dncols = 4
dnrows = 2

LON = np.pad(lon,[(dnrows,dnrows),(dncols,dncols)])

LAT = np.pad(lat,[(dnrows,dnrows),(dncols,dncols)])
print(LON)
strike = 293.0*(np.pi/180)

strike_lon = np.cos(np.pi/2-strike) 
strike_lat = np.sin(np.pi/2-strike)

dip_lon = np.sin(np.pi/2-strike)
dip_lat = -np.cos(np.pi/2-strike)


dstrike = abs(lon0[1]-lon0[0])
ddip = abs(lat0[1]-lat0[0])
dr = np.sqrt(dstrike**2+ddip**2)

# need to add new rows/columns based on immediately close patches making up edges 
# not a fixed point

# second-to-last along strike row
# for i in range(1,2):
#     for j in range(0, 21):
#         try:
#             lon = LON[i+1][j] + 0*dr*(j-dncols)*strike_lon + dr*(-1)*dip_lon
#             lat = LAT[i+1][j] + 0*dr*(j-dncols)*strike_lat + dr*(-1)*dip_lat
#             LON[i][j] = lon
#             LAT[i][j] = lat
#         except:
#             pass
# # last along strike row
# for i in range(0,1):
#     for j in range(0, 21):
#         try:
#             lon = LON[i+1][j] + 0*dr*(j-dncols)*strike_lon + dr*(-1)*dip_lon
#             lat = LAT[i+1][j] + 0*dr*(j-dncols)*strike_lat + dr*(-1)*dip_lat
#             LON[i][j] = lon
#             LAT[i][j] = lat
#         except:
#             pass
i = 1

lon = LON[i+1] + dr*(-1)*dip_lon
lat = LAT[i+1] + + dr*(-1)*dip_lat
LON[i]= lon
LAT[i] = lat     

i = 0

lon = LON[i+1] + dr*(-1)*dip_lon
lat = LAT[i+1] + dr*(-1)*dip_lat
LON[i] = lon
LAT[i] = lat        
# last along strike row
# for i in range(11,12):
#     for j in range(0, 21):
#         try:
#             lon = LON[i-1][j] + 0*dr*(j-dncols)*strike_lon + dr*(1)*dip_lon
#             lat = LAT[i-1][j] + 0*dr*(j-dncols)*strike_lat + dr*(1)*dip_lat
#             LON[i][j] = lon
#             LAT[i][j] = lat
#         except:
#             pass    

i = 11
lon = LON[i-1] +  dr*(1)*dip_lon
lat = LAT[i-1]  + dr*(1)*dip_lat
LON[i]= lon
LAT[i]= lat    

i = 12
lon = LON[i-1] + dr*(1)*dip_lon
lat = LAT[i-1] +  dr*(1)*dip_lat
LON[i]= lon
LAT[i]= lat        

# last along strike row
# for i in range(12,13):
#     for j in range(0, 21):
#         try:
#             lon = LON[i-1][j] + 0*dr*(j-dncols)*strike_lon + dr*(1)*dip_lon
#             lat = LAT[i-1][j] + 0*dr*(j-dncols)*strike_lat + dr*(1)*dip_lat
#             LON[i][j] = lon
#             LAT[i][j] = lat
#         except:
#             pass            

j = 3
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j+1] + dr*(-1)*strike_lon
lat = LAT[:,j+1] + dr*(-1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat     

j = 2
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j+1] + dr*(-1)*strike_lon
lat = LAT[:,j+1] + dr*(-1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat  
j = 1
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j+1] + dr*(-1)*strike_lon
lat = LAT[:,j+1] + dr*(-1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat     

j = 0
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j+1] + dr*(-1)*strike_lon
lat = LAT[:,j+1] + dr*(-1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat     

j = 22
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j-1] + dr*(1)*strike_lon
lat = LAT[:,j-1] + dr*(1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat     

j = 23
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j-1] + dr*(1)*strike_lon
lat = LAT[:,j-1] + dr*(1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat     



j = 24
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j-1] + dr*(1)*strike_lon
lat = LAT[:,j-1] + dr*(1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat     

j = 25
LON[:,j] = np.zeros(len(LON[:,j]))
lon = LON[:,j-1] + dr*(1)*strike_lon
lat = LAT[:,j-1] + dr*(1)*strike_lat
LON[:,j]= lon
LAT[:,j] = lat 


lon = LON.flatten()
lat = LAT.flatten()
ids =df.iloc[:,0].values
plt.scatter(lon,lat,color='red')
plt.xlim(84,87)
plt.ylim(27,28.8)

plt.show()
LON = np.flip(LON)
LAT = np.flip(LAT)
LON = np.flip(LON,axis=1)
LAT = np.flip(LAT,axis=1)
lat = LAT.flatten()
lon = LON.flatten()
dict_d = {'lon':lon,'lat':lat}
df_d = pd.DataFrame(dict_d)
df_d.to_csv('Gorkha_displacement.csv')
# for i,id0 in enumerate(ids):
#     plt.annotate(str(int(id0)),(lon[i],lat[i]),(lon[i],lat[i]),fontsize=8)
# plt.axis('equal')




def extending_grid(name, strike,ncols,nrows,dncols,dnrows):
    df = pd.read_csv(name)

    lon0 = df['lon'].values
    lat0 = df['lat'].values
    lon =lon0.reshape(nrows,ncols,order='F')
    lat = lat0.reshape(nrows,ncols,order='F')

    LON = np.pad(lon,[(dnrows,dnrows),(dncols,dncols)])
    LAT = np.pad(lat,[(dnrows,dnrows),(dncols,dncols)])
    
    strike = strike*(np.pi/180)

    strike_lon = np.cos(np.pi/2-strike) 
    strike_lat = np.sin(np.pi/2-strike)

    dip_lon = np.sin(np.pi/2-strike)
    dip_lat = -np.cos(np.pi/2-strike)

    dstrike = abs(lon0[1]-lon0[0])
    ddip = abs(lat0[1]-lat0[0])
    dr = np.sqrt(dstrike**2+ddip**2)

    
    iref = dnrows
    for i in range(0,iref):
        
        LON[i][dncols:-dncols] = LON[iref][dncols:-dncols] + (iref - i)*dr*(-1)*dip_lon
        LAT[i][dncols:-dncols] = LAT[iref][dncols:-dncols] +  (iref - i)*dr*(-1)*dip_lat
        
        LON[-i-1][dncols:-dncols] = LON[-iref-1][dncols:-dncols] + (iref - i)*dr*(1)*dip_lon
        LAT[-i-1][dncols:-dncols] = LAT[-iref-1][dncols:-dncols] +  (iref - i)*dr*(1)*dip_lat
    
    jref = dncols
    for j in range(0,jref):
        LON[:,j] = LON[:,jref] + (jref - j)*dr*(-1)*strike_lon
        LAT[:,j] = LAT[:,jref] +  (jref - j)*dr*(-1)*strike_lat
        
        LON[:,-j-1] = LON[:,-jref-1] + (jref - j)*dr*(1)*strike_lon
        LAT[:,-j-1] = LAT[:,-jref-1] +  (jref - j)*dr*(1)*strike_lat
        
    
    LON = np.flip(LON)
    LAT = np.flip(LAT)
    LON = np.flip(LON,axis=1)
    LAT = np.flip(LAT,axis=1)
    lat = LAT.flatten()
    lon = LON.flatten()
    dict_d = {'lon':lon,'lat':lat}
    df_d = pd.DataFrame(dict_d)
    df_d.to_csv('Tohoku_displacement.csv')
    
    return LON, LAT
    
lon ,lat = extending_grid('Tohoku_geometry.csv',194,24,9,0,0)

plt.scatter(lon,lat)
# plt.scatter(lon[0][0],lat[0][0],marker='x',color='k')
# plt.scatter(lon[0][-1],lat[0][-1],marker='x',color='k')
# plt.scatter(lon[-1][0],lat[-1][0],marker='x',color='k')
# plt.scatter(lon[-1][-1],lat[-1][-1],marker='x',color='k')
plt.axis('equal')

p1 = (lon[0][0],lat[0][0])
p2 = (lon[0][-1],lat[0][-1])
p3 = (lon[-1][0],lat[-1][0])
p4 = (lon[-1][-1],lat[-1][-1])
dir_strike = np.array([p4[0]-p3[0],p4[1]-p3[1]])
dir_dip = np.array([p4[0]-p2[0],p4[1]-p2[1]])
 
split1= np.arange(0,1 + 0.20,0.20)
split2 = np.arange(0,1 + 0.25,0.25)
branch1_x = p3[0] + split1*dir_strike[0]
branch1_y = p3[1] + split1*dir_strike[1]

branch2_x = p2[0] + split2*dir_dip[0]
branch2_y = p2[1] + split2*dir_dip[1]


branch1 = np.column_stack((branch1_x,branch1_y))

branch2 = np.column_stack((branch2_x,branch2_y))

# plt.scatter(branch1[:,0] - 1*dir_dip[0],branch1[:,1] - 1*dir_dip[1])
# plt.scatter(branch2[:,0],branch2[:,1])
stations = []
for sp in split2:
    st = branch1 - sp*dir_dip
    stations.append(st)
    # plt.scatter(st[:,0],st[:,1])
    
stations = np.array(stations)
stations = stations.reshape(stations.shape[0]*stations.shape[1],2)
for i in range(len(stations)):
    plt.scatter(stations[i,0],stations[i,1],ec='k')
    plt.annotate(f'{i+1}',(stations[i,0],stations[i,1]))


np.savetxt('tohoku_stations_mudpy.txt',stations,fmt='%.6f')



# ncols = 24
# nrows= 9 

# name = 'Tohoku_geometry.csv'
# df = pd.read_csv(name)

# lon0 = df['lon'].values
# lat0 = df['lat'].values
# lon =lon0.reshape(nrows,ncols,order='F')
# lat = lat0.reshape(nrows,ncols,order='F')

# dncols = 6
# dnrows = 2

# LON = np.pad(lon,[(dnrows,dnrows),(dncols,dncols)])

# LAT = np.pad(lat,[(dnrows,dnrows),(dncols,dncols)])

# strike = 194.0*(np.pi/180)

# strike_lon = np.cos(np.pi/2-strike) 
# strike_lat = np.sin(np.pi/2-strike)

# dip_lon = np.sin(np.pi/2-strike)
# dip_lat = -np.cos(np.pi/2-strike)


# dstrike = abs(lon0[1]-lon0[0])
# ddip = abs(lat0[1]-lat0[0])
# dr = np.sqrt(dstrike**2+ddip**2)

# # need to add new rows/columns based on immediately close patches making up edges 
# # not a fixed point

# # second-to-last along strike row
# # for i in range(1,2):
# #     for j in range(0, 21):
# #         try:
# #             lon = LON[i+1][j] + 0*dr*(j-dncols)*strike_lon + dr*(-1)*dip_lon
# #             lat = LAT[i+1][j] + 0*dr*(j-dncols)*strike_lat + dr*(-1)*dip_lat
# #             LON[i][j] = lon
# #             LAT[i][j] = lat
# #         except:
# #             pass
# # # last along strike row
# # for i in range(0,1):
# #     for j in range(0, 21):
# #         try:
# #             lon = LON[i+1][j] + 0*dr*(j-dncols)*strike_lon + dr*(-1)*dip_lon
# #             lat = LAT[i+1][j] + 0*dr*(j-dncols)*strike_lat + dr*(-1)*dip_lat
# #             LON[i][j] = lon
# #             LAT[i][j] = lat
# #         except:
# #             pass
# i = 1

# lon = LON[i+1] + dr*(-1)*dip_lon
# lat = LAT[i+1] + + dr*(-1)*dip_lat
# LON[i]= lon
# LAT[i] = lat     

# i = 0

# lon = LON[i+1] + dr*(-1)*dip_lon
# lat = LAT[i+1] + dr*(-1)*dip_lat
# LON[i] = lon
# LAT[i] = lat        
# # last along strike row
# # for i in range(11,12):
# #     for j in range(0, 21):
# #         try:
# #             lon = LON[i-1][j] + 0*dr*(j-dncols)*strike_lon + dr*(1)*dip_lon
# #             lat = LAT[i-1][j] + 0*dr*(j-dncols)*strike_lat + dr*(1)*dip_lat
# #             LON[i][j] = lon
# #             LAT[i][j] = lat
# #         except:
# #             pass    

# i = 11
# lon = LON[i-1] +  dr*(1)*dip_lon
# lat = LAT[i-1]  + dr*(1)*dip_lat
# LON[i]= lon
# LAT[i]= lat    

# i = 12
# lon = LON[i-1] + dr*(1)*dip_lon
# lat = LAT[i-1] +  dr*(1)*dip_lat
# LON[i]= lon
# LAT[i]= lat        

# # last along strike row
# # for i in range(12,13):
# #     for j in range(0, 21):
# #         try:
# #             lon = LON[i-1][j] + 0*dr*(j-dncols)*strike_lon + dr*(1)*dip_lon
# #             lat = LAT[i-1][j] + 0*dr*(j-dncols)*strike_lat + dr*(1)*dip_lat
# #             LON[i][j] = lon
# #             LAT[i][j] = lat
# #         except:
# #             pass            


# j = 5
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j+1] + dr*(-1)*strike_lon
# lat = LAT[:,j+1] + dr*(-1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 4
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j+1] + dr*(-1)*strike_lon
# lat = LAT[:,j+1] + dr*(-1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat  

# j = 3
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j+1] + dr*(-1)*strike_lon
# lat = LAT[:,j+1] + dr*(-1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 2
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j+1] + dr*(-1)*strike_lon
# lat = LAT[:,j+1] + dr*(-1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat  
# j = 1
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j+1] + dr*(-1)*strike_lon
# lat = LAT[:,j+1] + dr*(-1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 0
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j+1] + dr*(-1)*strike_lon
# lat = LAT[:,j+1] + dr*(-1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     



# j = 30
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j-1] + dr*(1)*strike_lon
# lat = LAT[:,j-1] + dr*(1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 31
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j-1] + dr*(1)*strike_lon
# lat = LAT[:,j-1] + dr*(1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat 

# j = 32
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j-1] + dr*(1)*strike_lon
# lat = LAT[:,j-1] + dr*(1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 33
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j-1] + dr*(1)*strike_lon
# lat = LAT[:,j-1] + dr*(1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 34
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j-1] + dr*(1)*strike_lon
# lat = LAT[:,j-1] + dr*(1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     

# j = 35
# LON[:,j] = np.zeros(len(LON[:,j]))
# lon = LON[:,j-1] + dr*(1)*strike_lon
# lat = LAT[:,j-1] + dr*(1)*strike_lat
# LON[:,j]= lon
# LAT[:,j] = lat     


# lon = LON.flatten()
# lat = LAT.flatten()
# ids =df.iloc[:,0].values
# plt.scatter(lon,lat,color='red')
# plt.xlim(135,148)
# plt.ylim(33,45)

# plt.show()
# # LON = np.flip(LON)
# # LAT = np.flip(LAT)
# # LON = np.flip(LON,axis=1)
# # LAT = np.flip(LAT,axis=1)
# # lat = LAT.flatten()
# # lon = LON.flatten()
# # dict_d = {'lon':lon,'lat':lat}
# # df_d = pd.DataFrame(dict_d)
# # df_d.to_csv('Tohoku_displacement.csv')
