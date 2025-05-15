# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 09:55:26 2024

@author: joanv
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 09:12:52 2024

@author: joanv
"""


import pygmt 
import pandas as pd
import numpy as np
import os

from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import h5py as h5
import io
import os

working_dir = os.getcwd()

def two_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0)
    
def one_array_formatter(array,shape):
        return np.flip(array.reshape(shape,order='F'),axis=0).flatten()
    
def set_stn(c,factor,patch,control=0):
      
      dpatch = patch/c
      offsetx =  (ncols//factor)*patch
      offsety =  (nrows//factor)*patch
      xstn = np.arange(-control*offsetx + dpatch/2, c*ncols*dpatch + control*offsetx,dpatch)
      ystn = -np.arange(-control*offsety + dpatch/2, c*nrows*dpatch + control*offsety,dpatch)
      ystn = np.flip(ystn)
      return xstn, ystn
  
def proj_ysrc_coords(patch,dip):
      proj_dysrc = -patch*np.cos(dip*np.pi/180) # in meters
      proj_ysrc = np.zeros_like(proj_dysrc)
      for i in range(len(proj_ysrc)):
          proj_ysrc[i] = sum(proj_dysrc[:i]) + (1/2)*proj_dysrc[i] 
      ysrc = np.flip(proj_ysrc)
        
      return ysrc


names = ['Tohoku','Iquique','Illapel','Pedernales','Gorkha']
geoms = [(9,24),(11,12),(10,17),(8,10),(9,18)]
patches = [29,17,18,15,10]
arrow_sizes = [10,5,5,5,5]
nparams = [866,533,682,650,331]
rakes = [90,0,0,0,0]
z_offsets = [2.8,0,0,0,0]
nramps = [0,0,0,0,0]
def build_model(names,geoms,patches,arrow_sizes,nparams,rakes):
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


def get_mu(model_parameters,depth):
        heights = model_parameters['H']
        depths = np.cumsum(heights)
        i = np.argmin(abs(depth - depths))
        if depth > depths[i]:
              i +=1 
        rho = model_parameters['RHO'][i]*1e3 # convert to SI
        vs = model_parameters['VP'][i]*1e3    # convert to SI
        mu = rho*vs**2 
        return mu

models = build_model(names,geoms,patches,arrow_sizes,nparams,rakes)

km = 1000

name = 'Iquique'
model_type = 'kinematic'
nsamples = 100
z_offset = models[name]['z_offset']



nrows,ncols = models[name]['geom'][0],models[name]['geom'][1]
nparam = models[name]['nparam']

shape = (nrows,ncols)
patch = models[name]['patch']
nramp = models[name]['nramp']




    
df = pd.read_csv(f'INPUT/{name}/model/kinematic/all_samples/mean/{name}_mean_kinematic_model.csv')



LON = one_array_formatter(df['lon'].values,shape)
LAT = one_array_formatter(df['lat'].values,shape)

ONSET = np.loadtxt(f'RuptureTime_{name}_all_kinematic_model.txt')


trench = np.loadtxt('southamerica.lonlat')

ONSET = ONSET.reshape(shape)
# df.insert(8, "onset", onset, True)
# df.insert(5,"type",rise_type,True)

# L = np.ones_like(onset)*patchsize*1e3
# W = np.ones_like(onset)*patchsize*1e3

# df.insert(10, "L(m)", L, True)
# df.insert(11,"W(m)",W,True)


Uperp = two_array_formatter(df['U_perp'].values,shape)
Uparallel =  two_array_formatter(df['U_parallel'].values,shape)
SLIP = two_array_formatter(df['Slip'].values,shape)
RAKE = two_array_formatter(np.arctan2(Uparallel,Uperp)*(180/np.pi),shape)
DIP = two_array_formatter(df['dip'].values,shape)
# vertical shift for seafloor
DEPTH = two_array_formatter(df['depth'].values,shape)
# DEPTH = array_formatter(df['depth'].values,shape) 
DURATION = two_array_formatter(df['Tr'].values,shape)

LON = two_array_formatter(df['lon'].values,shape)
LAT = two_array_formatter(df['lat'].values,shape)
STRIKE = two_array_formatter(df['strike'].values,shape)



patch_interp = patch/3

xS_interp = np.arange(patch_interp/2 , ncols*(1.015)*patch,patch_interp)
yS_interp = np.arange(-(nrows*0.97*patch) + patch_interp/2,0 + 0.7*patch_interp,patch_interp)
XS_interp,YS_interp  = np.meshgrid(xS_interp,yS_interp)

xS = np.arange(patch/2 , ncols*patch,patch)
yS = np.arange(-(nrows-1/2)*patch,0,patch)
   # shift accordingly at surface
yS = proj_ysrc_coords(patch,df['dip'].values[:nrows])
XS, YS = np.meshgrid(xS,yS)


DIP = np.flip(df['dip'].values.reshape(nrows,ncols,order='F'),axis=0)
DEPTH = np.flip(df['depth'].values.reshape(nrows,ncols,order='F'),axis=0) 
method = 'linear'
dip_interp = RegularGridInterpolator((yS, xS), DIP,bounds_error=False, fill_value=None,method = method)
strike_interp = RegularGridInterpolator((yS, xS), STRIKE,bounds_error=False, fill_value=None,method = method)

depth_interp = RegularGridInterpolator((yS, xS), DEPTH,bounds_error=False, fill_value=None,method = method)

lat_interp = RegularGridInterpolator((yS, xS), LAT,bounds_error=False, fill_value=None,method = method)
lon_interp = RegularGridInterpolator((yS, xS), LON,bounds_error=False, fill_value=None,method = method)
onset_interp = RegularGridInterpolator((yS, xS), ONSET,bounds_error=False, fill_value=None,method = method)

DIP = dip_interp((YS_interp,XS_interp))
DEPTH = depth_interp((YS_interp,XS_interp))
LON = lon_interp((YS_interp,XS_interp))
LAT = lat_interp((YS_interp,XS_interp))
ONSET = onset_interp((YS_interp,XS_interp))
# pygmt.config(COLOR_BACKGROUND = 'white',
#               MAP_GRID_PEN_PRIMARY = '0.3p,dimgrey',
#               MAP_ANNOT_OFFSET_PRIMARY = '5p',
#               MAP_ANNOT_OFFSET_SECONDARY = '5p', 
#               MAP_ANNOT_OBLIQUE = '30',
#               FONT_ANNOT_PRIMARY = '10p,Helvetica', 
#               FONT_LABEL = '25p',
#               MAP_FRAME_WIDTH = '2p',
#               MAP_FRAME_PEN = '1.2p',
#               MAP_FRAME_TYPE = 'plain',
#               MAP_TICK_LENGTH_PRIMARY = '12p',
#               MAP_LABEL_OFFSET = '5.5p',
#               FORMAT_GEO_MAP = 'F')


boxes = {'Tohoku':[140, 145.5, 35, 42],
         'Gorkha':[84,86.5,27,29],
         'Iquique':[-72.5,-69,-22,-18.0],
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

for event in ['Iquique']:
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
phi = -13.58
strike0 = (phi)*np.pi/180
data={
    "x": vectors[:,0],
    "y": vectors[:,1],
    "east_velocity": df0['U_perp'].values*np.sin(strike0) - df0['U_parallel'].values*np.cos(strike0),
    "north_velocity": df0['U_perp'].values*np.cos(strike0) + df0['U_parallel'].values*np.sin(strike0),
    "east_sigma": np.max(np.column_stack((df0['std_U_perp'].values,df0['std_U_parallel'].values)),axis=1),
    "north_sigma":np.min(np.column_stack((df0['std_U_perp'].values,df0['std_U_parallel'].values)),axis=1),
    "correlation_EN": np.arctan2(df0['U_parallel'].values,df0['U_perp'].values)*(180/np.pi) + 90 -phi
}


def mapper(r,phi=phi):
    if r>1:
        return 180-phi
    else:
        return 90 - phi
    

ratio = df0['std_U_parallel'].values/df0['std_U_perp'].values
data['correlation_EN'] = np.array(list(map(mapper,ratio)))
df = pd.DataFrame(
    data=data
)
fig.basemap(
    region = region,
    map_scale="jBL+o1c/1c+w50k",
    projection = "M20c",
    frame = ["a0.5f0.5",f"+t Mw8.1 {event}"])
fig.coast(shorelines="0.5p,black",land="gray",water="lightskyblue",map_scale="jBL+o1c/1c+w50k")
pygmt.makecpt(
           transparency = 80,
           cmap="hot",
           series=[0,20.5,4],
           continuous=True,
           reverse = True
       )


patch = 17
nrows = 11
has = df0['Hypo_as'].values[0]
hdip = nrows*patch-df0['Hypo_dd'].values[0]

#  equivalent
# hlon = has*np.cos((90-phi)*np.pi/180) + hdip*np.cos((180-phi)*np.pi/180)
# hlat = has*np.sin((90-phi)*np.pi/180) + hdip*np.sin((180-phi)*np.pi/180)

# hlat = hdip*np.cos((phi-90)*np.pi/180) + has*np.cos(phi*np.pi/180)
# hlon = hdip*np.sin((phi-90)*np.pi/180) + has*np.sin(phi*np.pi/180)

hlat = hdip*np.sin(phi*np.pi/180) + has*np.cos(phi*np.pi/180)
hlon = -hdip*np.cos(phi*np.pi/180) + has*np.sin(phi*np.pi/180)

offlon  = patch*np.sin(phi*np.pi/180) + patch*np.cos(phi*np.pi/180)
offlat = patch*np.cos(phi*np.pi/180) + patch*np.sin(phi*np.pi/180)

hlat = df0['lat'].values[nrows-1] + (hlat + offlat/2) * (1/(111)) 
hlon = df0['lon'].values[nrows-1] + (hlon + offlon/2) * (1/(111))

hlat =  -19.621960505766026
hlon = -70.92153656005858 



#hlon, hlat = -70.92164515955,-19.62216108
fig.plot(data=data0,style=f"j{strike}/{patchsize}/{patchsize}",cmap =True, pen="0.2p,black",transparency = 80)
fig.velo(
    data=df,
    region=region,
    pen="0.75p,black",
    uncertaintyfill="lightblue1",
    line=True,
    spec="r0.105/0.67/18",
    projection="M20c",
    vector="0.4c+p1p+e+gblack",
)

fig.contour(
    region=region,
    projection="M20c",
    frame=["a"],
    pen="0.5p",
    # pass the data as 3 1-D data columns
    x= LON.flatten() ,
    y= LAT.flatten(),
    z= ONSET.flatten(),
    # set the contours z values intervals to 10
    levels=10,
    # set the contours annotation intervals to 20
    annotation=20,
)
fig.plot(x = hlon,y = hlat,style="a1.25c", pen="1p,black", fill="blue")
fig.plot(
    x=trench[:,0],  # Longitude in degrees East
    y=trench[:,1],  # Latitude in degrees North
    # Draw a 2-points thick, red, dashed line for the survey line
    pen="3.5p",
    style="f1c/0.25c+r+t+o0.3c+p",
    fill="black"
)

fig.coast(shorelines="0.5p,black")
fig.colorbar(frame="xa4f2+lSlip (m)",position='jBR+o0.6c/2.2c+w3c/0.2c+h',box='+gwhite+p0.8p,black')
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
    #fig.savefig('iquique_ref.pdf')
    fig.savefig('iquique_ref.png',dpi=1200)
    fig.show()

