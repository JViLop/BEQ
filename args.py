

names = ['Tohoku','Illapel','Iquique','Pedernales','Gorkha']
data_type = ['binary','h5','h5','h5','h5']
geom_size = [(9,24),(10,17),(11,12),(8,10),(9,18)]
patch_side = [29,18,17,15,10]
nparams = [866,682,533,331,650]
nrealizations = [1000000,393216,72000,600000,90000]
ramp = [0,0,3,9,0]
model_size = [(866,1000000),(682,393216),(533,72000),(331,600000),(650,90000)]
mw = [9.0,8.3,8.1,7.8,7.8]
nrows_skip = [2,1,1,1,1]
RotAngle = [None,None,None,360-99,None]
rake = [90,90,90,90,107]
header = [{'lon':0,'lat':1,'depth':4,'strike':5,'dip':6},{'lon':0,'lat':1,'depth':4,'strike':5,'dip':6},{'lon':0,'lat':1,'depth':4,'strike':5,'dip':6},{'lon':0,'lat':1,'depth':4,'strike':5,'dip':6},{'lon':2,'lat':1,'depth':3,'strike':4,'dip':5}]
step = [52,52,56,53,12]
def model_args(names,data_type,geom_size,patch_side,nparams,ramp,model_size,mw,nrows_skip,RotAngle,rake,header,step):
    config = dict()
    for i,n in enumerate(names):
        config[n] = {}
        config[n]['data_type'] = data_type[i]
        config[n]['nrows'] = geom_size[i][0]
        config[n]['ncols'] = geom_size[i][1]
        config[n]['patch_length'] = patch_side[i]
        config[n]['nparams'] = nparams[i]
        config[n]['ramp'] = ramp[i]
        config[n]['model_size'] = model_size[i]
        config[n]['mw'] = mw[i]
        config[n]['nrows_skip'] = nrows_skip[i]
        config[n]['RotAngle'] = RotAngle[i]
        config[n]['rake'] = rake[i]
        config[n]['header'] = header[i]
        config[n]['step'] = step[i]

    return config


models =  model_args(names,data_type,geom_size,patch_side,nparams,ramp,model_size,mw,nrows_skip,RotAngle,rake,header,step)
