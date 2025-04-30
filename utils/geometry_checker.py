# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 13:41:12 2023

@author: joanv
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

geometry_path = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Pedernales/geometry/Ecuador_8x10_END_lonlat.txt'

df = np.loadtxt(geometry_path,skiprows=1)
head = ['lon','lat','E','N','depth','strike','dip','Area','ID']

df = pd.DataFrame(df,columns = head)
df = df.drop(['E','N','Area'],axis=1)

lat = df['lat'].values

lon = df['lon'].values

ids = df['ID'].values
plt.scatter(lon,lat,color='red')
for i,id0 in enumerate(ids):
    plt.annotate(str(int(id0)),(lon[i],lat[i]),fontsize=8)
plt.title('Pedernales kinematic')
plt.axis('equal')
plt.show() 

# data already reformatted #

pathTohoku = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Tohoku/geometry/geometry_formatted/Tohoku_geometry.csv'
pathIquique = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Iquique/geometry/geometry_formatted/Iquique_geometry.csv'
pathIllapel = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Illapel/geometry/geometry_formatted/Illapel_geometry.csv'
pathGorkha = 'C:/Users/joanv/OneDrive/Escritorio/University_of_Oklahoma/GRA/EQ_source_models/EQ_source_models/EQ/Gorkha/geometry/geometry_formatted/Gorkha_geometry.csv'

def geom_checker(path,name,model_type='kinematic'):
    df = pd.read_csv(path)
    lat = df['lat'].values
    
    lon = df['lon'].values
    
    ids =df.iloc[:,0].values
    plt.scatter(lon,lat,color='red')
    for i,id0 in enumerate(ids):
        plt.annotate(str(int(id0)),(lon[i],lat[i]),(lon[i],lat[i]),fontsize=8)
    plt.title(f'{name} {model_type}')
    plt.axis('equal')
    plt.show()
        
geom_checker(pathTohoku,'Tohoku')
geom_checker(pathIquique,'Iquique')
geom_checker(pathIllapel,'Illapel')
geom_checker(pathGorkha,'Gorkha')