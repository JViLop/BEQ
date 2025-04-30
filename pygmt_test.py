# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 08:41:36 2024

@author: joanv
"""

import pygmt 
import pandas as pd
import numpy as np
import os
boxes = {'Tohoku':[140, 146, 35, 42],
         'Gorkha':[84,86.5,27,29],
         'Iquique':[-73,-69,-21,-18.5],
         'Illapel':[-74,-70,-33,-29],
         'Pedernales':[-81.5,-79,-1,1]}

sizes = {'Tohoku':0.65,
         'Gorkha':0.53,
         'Iquique':0.5,
         'Illapel':0.6,
         'Pedernales':0.9}

strike = {'Tohoku':193,
         'Gorkha':3,
         'Iquique':0.5,
         'Illapel':0.6,
         'Pedernales':0.9}

pwd = os.getcwd()
figure_folder  = os.path.join(pwd,'georeferenced_figs')
os.makedirs(figure_folder,exist_ok=True)
# for event in ['Tohoku','Gorkha','Iquique','Illapel','Pedernales']:
    
#     df = pd.read_csv(f'INPUT/{event}/model/kinematic/all_samples/mean/{event}_mean_kinematic_model.csv')
#     strike = df['strike'].values[0] 
    
#     strike = 90 - strike
#     patchsize =sizes[event]
    
#     df = df[['lon','lat','Slip']]
#     data = np.array(df)
#     region_map = boxes[event]
#     fig = pygmt.Figure()
#     pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
#     fig.basemap(
#         region=region_map,
#         projection="M12c",  # Mercator projection with a width of 12 centimeters
#         frame=["a1", f"+t{event}"],
#     )
#     grid_map = pygmt.datasets.load_earth_relief(
#         resolution="30s",
#         region=region_map,
#     )
#     fig.grdimage(grid=grid_map, cmap="oleron")
    

#     # Plot the downloaded grid with color-coding based on the elevation
    
    
#     # Ellipse
#     # [[lon, lat, direction, width, height]]
    
#     # data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
#     # fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)
    
    
#     # x = np.array([144.857,144.524,144.190])
#     # y = np.array([41.247,41.315,41.383])
    
#     pygmt.makecpt(
#         transparency = 80,
#         cmap="hot",
#         series=[df["Slip"].min(), df["Slip"].max(),0.5],
#         continuous=True,
#         reverse = True
#     )
    
#     # fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
#     fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)
    
    
#     # fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')

#     fig.colorbar(frame="x+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')
#     fig.coast(borders="1/0.5p,black",shorelines=True)
#     fig_dir = os.path.join(figure_folder,f'{event}.png')
#     fig.show()
#     fig.savefig(fig_dir)
    
    
Gorkha_stn = pd.read_csv('Gorkha_displacement.csv')
d = np.loadtxt('Gorkha_all_kinematic_model_displacement.txt')
Gorkha_stn['z'] = d[:,2]
strike = 293
strike =90 - strike
patchsize = 0.37
df = Gorkha_stn[['lon','lat','z']]
data = np.array(df)
region_map = [83.5,87,26.75,29.25]
fig = pygmt.Figure()
pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
fig.basemap(
    region=region_map,
    projection="M12c",  # Mercator projection with a width of 12 centimeters
    frame=["a1", "+tGorkha Uplift"],
)
grid_map = pygmt.datasets.load_earth_relief(
    resolution="30s",
    region=region_map,
)
fig.grdimage(grid=grid_map, cmap="oleron")


# Plot the downloaded grid with color-coding based on the elevation


# Ellipse
# [[lon, lat, direction, width, height]]

# data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
# fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)


# x = np.array([144.857,144.524,144.190])
# y = np.array([41.247,41.315,41.383])

pygmt.makecpt(
    transparency = 80,
    cmap="polar",
    series=[-2,2,0.5],
    continuous=True,
    reverse = False
)

# fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)


# fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')

fig.colorbar(frame="x+lUplift(m)",position='JMB+o0c/1c+w5c/0.2c+h')
fig.coast(borders="1/0.5p,black",shorelines=True)
fig_dir = os.path.join(figure_folder,'Gorkha_z_displacement.png')
fig.show()
fig.savefig(fig_dir)

Gorkha_stn = pd.read_csv('Tohoku_displacement.csv')
d = np.loadtxt('Tohoku_all_kinematic_model_displacement.txt')
Gorkha_stn['z'] = d[:,2]
strike = 194
strike =90 - strike
patchsize = 0.58
df = Gorkha_stn[['lon','lat','z']]
data = np.array(df)
region_map = [138,147,33,44]
fig = pygmt.Figure()
pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
fig.basemap(
    region=region_map,
    projection="M12c",  # Mercator projection with a width of 12 centimeters
    frame=["a1", "+tTohoku Uplift"],
)
grid_map = pygmt.datasets.load_earth_relief(
    resolution="30s",
    region=region_map,
)
fig.grdimage(grid=grid_map, cmap="oleron")


# Plot the downloaded grid with color-coding based on the elevation


# Ellipse
# [[lon, lat, direction, width, height]]

# data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
# fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)


# x = np.array([144.857,144.524,144.190])
# y = np.array([41.247,41.315,41.383])

pygmt.makecpt(
    transparency = 80,
    cmap="polar",
    series=[-8,8,0.5],
    continuous=True,
    reverse = False
)

# fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)


# fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')

fig.colorbar(frame="x+lUplift(m)",position='JMB+o0c/1c+w5c/0.2c+h')
fig.coast(borders="1/0.5p,black",shorelines=True)
fig_dir = os.path.join(figure_folder,'Tohoku_z_displacement.png')
fig.show()
fig.savefig(fig_dir)




for event in ['Tohoku','Gorkha','Iquique','Illapel','Pedernales']:
    
    eq = pd.read_csv(f'{event}_displacement.csv')
    d = np.loadtxt(f'{event}_all_kinematic_model_displacement.txt')
    eq['z'] = d[:,2]
    strike = strike[event]
    strike =90 - strike
    patchsize = 0.37
    df = eq[['lon','lat','z']]
    data = np.array(df)
    region_map = boxes[event]
    fig = pygmt.Figure()
    pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
    fig.basemap(
        region=region_map,
        projection="M12c",  # Mercator projection with a width of 12 centimeters
        frame=["a1", f"+t{event} Uplift"],
    )
    grid_map = pygmt.datasets.load_earth_relief(
        resolution="30s",
        region=region_map,
    )
    fig.grdimage(grid=grid_map, cmap="oleron")
    
    
    # Plot the downloaded grid with color-coding based on the elevation
    
    
    # Ellipse
    # [[lon, lat, direction, width, height]]
    
    # data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
    # fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)
    
    
    # x = np.array([144.857,144.524,144.190])
    # y = np.array([41.247,41.315,41.383])
    
    pygmt.makecpt(
        transparency = 80,
        cmap="polar",
        series=[-max(df['z'].values),max(df['z'].values),0.5],
        continuous=True,
        reverse = False
    )
    
    # fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
    fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)
    
    
    # fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')
    
    fig.colorbar(frame="x+lUplift(m)",position='JMB+o0c/1c+w5c/0.2c+h')
    fig.coast(borders="1/0.5p,black",shorelines=True)
    fig_dir = os.path.join(figure_folder,f'{event}_z_displacement.png')
    fig.show()
    fig.savefig(fig_dir)
    
# Gorkha_stn = pd.read_csv('Tohoku_displacement.csv')
# d = np.loadtxt('Tohoku_all_kinematic_model_displacement.txt')
# Gorkha_stn['z'] = d[:,2]
# strike = 194
# strike =90 - strike
# patchsize = 0.58
# df = Gorkha_stn[['lon','lat','z']]
# data = np.array(df)
# region_map = [138,147,33,44]
# fig = pygmt.Figure()
# pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
# fig.basemap(
#     region=region_map,
#     projection="M12c",  # Mercator projection with a width of 12 centimeters
#     frame=["a1", "+tTohoku Uplift"],
# )
# grid_map = pygmt.datasets.load_earth_relief(
#     resolution="30s",
#     region=region_map,
# )
# fig.grdimage(grid=grid_map, cmap="oleron")


# # Plot the downloaded grid with color-coding based on the elevation


# # Ellipse
# # [[lon, lat, direction, width, height]]

# # data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
# # fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)


# # x = np.array([144.857,144.524,144.190])
# # y = np.array([41.247,41.315,41.383])

# pygmt.makecpt(
#     transparency = 80,
#     cmap="polar",
#     series=[-8,8,0.5],
#     continuous=True,
#     reverse = False
# )

# # fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
# fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)


# # fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')

# fig.colorbar(frame="x+lUplift(m)",position='JMB+o0c/1c+w5c/0.2c+h')
# fig.coast(borders="1/0.5p,black",shorelines=True)
# fig_dir = os.path.join(figure_folder,'Tohoku_z_displacement.png')
# fig.show()
# fig.savefig(fig_dir)
    





