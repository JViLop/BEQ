# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 23:47:49 2024

@author: joanv
"""

""


import pygmt 
import pandas as pd
import numpy as np
import os

pygmt.config(COLOR_BACKGROUND = 'white',
             MAP_GRID_PEN_PRIMARY = '0.3p,dimgrey',
             MAP_ANNOT_OFFSET_PRIMARY = '5p',
             MAP_ANNOT_OFFSET_SECONDARY = '5p', 
             MAP_ANNOT_OBLIQUE = '30',
             FONT_ANNOT_PRIMARY = '8p,Helvetica', 
             FONT_LABEL = '8p',
             MAP_FRAME_WIDTH = '2p',
             MAP_FRAME_PEN = '1.2p',
             MAP_FRAME_TYPE = 'plain',
             MAP_TICK_LENGTH_PRIMARY = '12p',
             MAP_LABEL_OFFSET = '5.5p',
             FORMAT_GEO_MAP = 'F')


boxes = {'Tohoku':[140, 145.5, 35, 42],
         'Gorkha':[84,86.5,27,29],
         'Iquique':[-72.5,-69,-21,-18.5],
         'Illapel':[-74,-70,-32.5,-29.5],
         'Pedernales':[-81.5,-79,-1,1]}



hypocenter = {'Tohoku':[140],
         'Gorkha':[84,86.5,27,29],
         'Iquique':[-73,-69,-21,-18.5],
         'Illapel':[-74,-70,-32.5,-29.5],
         'Pedernales':[-81.5,-79,-1,1]}


sizes = {'Tohoku':1.285,
         'Gorkha':0.867,
         'Iquique':1.1,
         'Illapel':1,
         'Pedernales':1.05}

strike = {'Tohoku':193,
         'Gorkha':3,
         'Iquique':0.5,
         'Illapel':0.6,
         'Pedernales':0.9}

dcw_label = {'Tohoku':'JP',
            'Gorkha':'NP',
            'Iquique':'CL',
            'Illapel':'CL',
            'Pedernales':'EC'
            }

# for event in list(boxes.keys()):
# for event in ['Gorkha']:

for event in ['Pedernales']:
    if event in ['Pedernales','Iquique','Illapel']:
        edge = 'T'
        side = 'L'
        continent = [-90, -50, -50,10]
        w =4
        h = 6.7
    elif event=='Gorkha':
        edge = 'T'
        side = 'R'
        continent = [60,100,10,35]
        w = 5
        h = 3.4
    elif event =='Tohoku':
        edge = 'B'
        side = 'R'
        continent = [110,150,15,50]
        w = 5.8
        h = 6.25
        
    fig = pygmt.Figure()
    
    region = boxes[event]
    fig.coast(
        region = region,
        shorelines = 'thin',
        borders=["1/1,black"],
        land="gray",
        water="lightskyblue",
        map_scale="jBL+o1c/1c+w50k",
        projection = "M20c",
        frame = ["af",f"+t {event}"])

    df0 = pd.read_csv(f'INPUT/{event}/model/kinematic/all_samples/mean/{event}_mean_kinematic_model.csv')
    strike = df0['strike'].values[0] 
    strike = 90 - strike
    patchsize =sizes[event]
    df = df0[['lon','lat','Slip']]
    data0 = np.array(df)
    # fig.plot(x = df['lon'],y = df['lat'],projection="M20c",style="c0.05c", pen="0.01p,black",transparency = 80)
    pygmt.makecpt(
            transparency = 80,
            cmap="hot",
            series=[df["Slip"].min(),7.5,1],
            continuous=True,
            reverse = True
        )

    fig.plot(data=data0,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)
    vectors = df0[['lon','lat','Slip','U_parallel','U_perp']]
    dir_vectors = np.arctan2(vectors['U_parallel'].values,vectors['U_perp'].values)*(180/np.pi) + (90-27.05)          
    vectors = np.array(vectors[['lon','lat']])
    vectors = np.column_stack((vectors,dir_vectors))
    vectors = np.column_stack((vectors,df0[['Slip']].values/2))
    fig.plot(
    region=region,
    frame="a",
    data=vectors,
    style="v0.4c+e",
    pen="1p",
    fill="black",
)

    with fig.inset(position=f"j{edge}{side}+w{w}c/{h}c+o0.2c/0.2c",box="+p4.0p,black"):
        # Use a plotting function to create a figure inside the inset
        # fig.basemap(
        #     region=[-100, -40, -50,30],
        #     frame='a'
        # )
        country = dcw_label[event]
        fig.coast(region=continent,shorelines = ["1/0.5,black"],land="gray",borders=["1/0.5,black"], water="lightskyblue",dcw=[f"{country}+gred"])
        rectangle = [[region[0], region[2], region[1], region[3]]]
        fig.plot(data=rectangle, style="r+s", pen="2p,white")
    # event='Pedernales'
    # df = pd.read_csv(f'INPUT/{event}/model/kinematic/all_samples/mean/{event}_mean_kinematic_model.csv')
    # strike = df['strike'].values[0] 
    # strike = 90 - strike
    # patchsize =sizes[event]
    # df = df[['lon','lat','Slip']]
    # data = np.array(df)
    # fig.plot(x = df['lon'],y = df['lat'],projection="M20c",style="c0.1c", pen="0.01p,black",transparency = 80)
    

        # Use dcw to selectively highlight an area
    

    # event='Pedernales'
    # df = pd.read_csv(f'INPUT/{event}/model/kinematic/all_samples/mean/{event}_mean_kinematic_model.csv')
    # strike = df['strike'].values[0] 
    # strike = 90 - strike
    # patchsize =sizes[event]
    # df = df[['lon','lat','Slip']]
    # data = np.array(df)
    # fig.plot(x = df['lon'],y = df['lat'],projection="M20c",style="c0.05c", pen="0.01p,black",transparency = 80)
    # fig.savefig()
    fig.show()
    # fig.savefig(f'{event}.pdf')


fig = pygmt.Figure()
strike0 = (27.05)*np.pi/180
df = pd.DataFrame(
    data={
        "x": vectors[:,0],
        "y": vectors[:,1],
        "east_velocity": df0['U_perp'].values*np.sin(strike0) - df0['U_parallel'].values*np.cos(strike0),
        "north_velocity": df0['U_perp'].values*np.cos(strike0) + df0['U_parallel'].values*np.sin(strike0),
        "east_sigma": df0['std_U_perp'].values*np.sin(strike0) - df0['std_U_parallel'].values*np.cos(strike0),
        "north_sigma":df0['std_U_perp'].values*np.cos(strike0) + df0['std_U_parallel'].values*np.sin(strike0),
        "correlation_EN": np.arctan2(df0['U_parallel'].values,df0['U_perp'].values)*(180/np.pi) + 90 -27.05
    }
)
fig.basemap(
    region = region,
    map_scale="jBL+o1c/1c+w50k",
    projection = "M20c",
    frame = ["af",f"+t {event}"])
fig.coast(shorelines="0.5p,black",land="gray",water="lightskyblue",map_scale="jBL+o1c/1c+w50k")
pygmt.makecpt(
           transparency = 80,
           cmap="hot",
           series=[df0["Slip"].min(),7.5,1],
           continuous=True,
           reverse = True
       )

fig.plot(data=data0,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)
fig.velo(
    data=df,
    region=region,
    pen="0.6p,black",
    uncertaintyfill="lightblue1",
    line=True,
    spec="r0.15/0.39/18",
    frame=["a"],
    projection="M20c",
    vector="0.5c+p1p+e+gblack",
)
# fig.colorbar(frame="x+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')
with fig.inset(position=f"j{edge}{side}+w{w}c/{h}c+o0.2c/0.2c",box="+p4.0p,black"):
    # Use a plotting function to create a figure inside the inset
    # fig.basemap(
    #     region=[-100, -40, -50,30],
    #     frame='a'
    # )
    country = dcw_label[event]
    fig.coast(region=continent,shorelines = ["1/0.5,black"],land="gray",borders=["1/0.5,black"], water="lightskyblue",dcw=[f"{country}+gred"])
    rectangle = [[region[0], region[2], region[1], region[3]]]
    fig.plot(data=rectangle, style="r+s", pen="2p,white")
# fig.savefig('gorkha_ref.png',dpi=500)

fig.show()








# pwd = os.getcwd()
# figure_folder  = os.path.join(pwd,'georeferenced_figs')
# os.makedirs(figure_folder,exist_ok=True)
# # for event in ['Tohoku','Gorkha','Iquique','Illapel','Pedernales']:
    
# #     df = pd.read_csv(f'INPUT/{event}/model/kinematic/all_samples/mean/{event}_mean_kinematic_model.csv')
# #     strike = df['strike'].values[0] 
    
# #     strike = 90 - strike
# #     patchsize =sizes[event]
    
# #     df = df[['lon','lat','Slip']]
# #     data = np.array(df)
# #     region_map = boxes[event]
# #     fig = pygmt.Figure()
# #     pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
# #     fig.basemap(
# #         region=region_map,
# #         projection="M12c",  # Mercator projection with a width of 12 centimeters
# #         frame=["a1", f"+t{event}"],
# #     )
# #     grid_map = pygmt.datasets.load_earth_relief(
# #         resolution="30s",
# #         region=region_map,
# #     )
# #     fig.grdimage(grid=grid_map, cmap="oleron")
    

# #     # Plot the downloaded grid with color-coding based on the elevation
    
    
# #     # Ellipse
# #     # [[lon, lat, direction, width, height]]
    
# #     # data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
# #     # fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)
    
    
# #     # x = np.array([144.857,144.524,144.190])
# #     # y = np.array([41.247,41.315,41.383])
    
# #     pygmt.makecpt(
# #         transparency = 80,
# #         cmap="hot",
# #         series=[df["Slip"].min(), df["Slip"].max(),0.5],
# #         continuous=True,
# #         reverse = True
# #     )
    
# #     # fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
# #     fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)
    
    
# #     # fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')

# #     fig.colorbar(frame="x+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')
# #     fig.coast(borders="1/0.5p,black",shorelines=True)
# #     fig_dir = os.path.join(figure_folder,f'{event}.png')
# #     fig.show()
# #     fig.savefig(fig_dir)
    
    
# Gorkha_stn = pd.read_csv('Gorkha_displacement.csv')
# d = np.loadtxt('Gorkha_all_kinematic_model_displacement.txt')
# Gorkha_stn['z'] = d[:,2]
# strike = 293
# strike =90 - strike
# patchsize = 0.37
# df = Gorkha_stn[['lon','lat','z']]
# data = np.array(df)
# region_map = [83.5,87,26.75,29.25]
# fig = pygmt.Figure()
# pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
# fig.basemap(
#     region=region_map,
#     projection="M12c",  # Mercator projection with a width of 12 centimeters
#     frame=["a1", "+tGorkha Uplift"],
# )
# grid_map = pygmt.datasets.load_earth_relief(
#     resolution="30s",
#     region=region_map,
# )
# fig.grdimage(grid=grid_map, cmap="oleron")


# # Plot the downloaded grid with color-coding based on the elevation

# # [[lon, lat, direction, width, height]]

# # data = [[144.857,41.247,256.0,0.65,0.65], [144.524,41.315,256,0.65,0.65], [144.190, 41.383,256,0.65,0.65],[144.013,40.878,256,0.65,0.65],[143.357,41.007,256,0.65,0.65]]
# # fig.plot(data=data, style="j", fill='red', pen="2p,black",transparency =50)


# # x = np.array([144.857,144.524,144.190])
# # y = np.array([41.247,41.315,41.383])

# pygmt.makecpt(
#     transparency = 80,
#     cmap="polar",
#     series=[-2,2,0.5],
#     continuous=True,
#     reverse = False
# )

# # fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
# fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)


# # fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')

# fig.colorbar(frame="x+lUplift(m)",position='JMB+o0c/1c+w5c/0.2c+h')
# fig.coast(borders="1/0.5p,black",shorelines=True)
# fig_dir = os.path.join(figure_folder,'Gorkha_z_displacement.png')
# fig.show()
# fig.savefig(fig_dir)

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




# for event in ['Tohoku','Gorkha','Iquique','Illapel','Pedernales']:
    
#     eq = pd.read_csv(f'{event}_displacement.csv')
#     d = np.loadtxt(f'{event}_all_kinematic_model_displacement.txt')
#     eq['z'] = d[:,2]
#     strike = strike[event]
#     strike =90 - strike
#     patchsize = 0.37
#     df = eq[['lon','lat','z']]
#     data = np.array(df)
#     region_map = boxes[event]
#     fig = pygmt.Figure()
#     pygmt.config(FONT="15p", FONT_ANNOT_PRIMARY="8p",FONT_ANNOT_SECONDARY="8p",MAP_TITLE_OFFSET="2p")
#     fig.basemap(
#         region=region_map,
#         projection="M12c",  # Mercator projection with a width of 12 centimeters
#         frame=["a1", f"+t{event} Uplift"],
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
#         cmap="polar",
#         series=[-max(df['z'].values),max(df['z'].values),0.5],
#         continuous=True,
#         reverse = False
#     )
    
#     # fig.plot(data = data,style=f"j{strike}/0.1/0.1",cmap =True, pen="0.75p,black",transparency = 75)
#     fig.plot(data = data,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.75p,black",transparency = 80)
    
    
#     # fig.colorbar(frame="xa20f10+lSlip(m)",position='JMB+o0c/1c+w5c/0.2c+h')
    
#     fig.colorbar(frame="x+lUplift(m)",position='JMB+o0c/1c+w5c/0.2c+h')
#     fig.coast(borders="1/0.5p,black",shorelines=True)
#     fig_dir = os.path.join(figure_folder,f'{event}_z_displacement.png')
#     fig.show()
#     fig.savefig(fig_dir)
    
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
    





