from args import models
from utils.model_reader import EQ_model

nsamples =  100

for name in list(models.keys()):
    data_type = models[name]['data_type']
    geom_size = (models[name]['nrows'],models[name]['ncols'])
    patch_length = models[name]['patch_length']
    ramp = models[name]['ramp']
    nrows_skip = models[name]['nrows_skip']
    model_size = models[name]['model_size']
    RotAngle = models[name]['RotAngle']
    rake = models[name]['rake']
    header = models[name]['header']
    step = models[name]['step']
    print(name)
    model = EQ_model(name,'kinematic',data_type,geom_size,patch_length,model_file_shape=model_size,nramp=ramp,nrows_skip=nrows_skip,RotAngle=RotAngle,header=header,rake=rake,sampling=True,nsamples=nsamples,step=step)
    model.plot_spatial_corr(0.2)
    model.plot_corr_matrix()
    model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
    model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
    model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
    model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
    model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)




'''

model = EQ_model('Tohoku','kinematic','binary',(9,24),29,nrows_skip=2,model_file_shape=(866,1000000),sampling=True, nsamples=100,step=52)
model.plot_spatial_corr(0.2)
model.plot_corr_matrix()
model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)


model = EQ_model('Illapel','kinematic','h5',(10,17),18,nramp=0,sampling=True, nsamples=100)
model.plot_spatial_corr(0.2)
model.plot_corr_matrix()
model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)

model = EQ_model('Iquique','kinematic','h5',(11,12),17,nramp=3, sampling=True,nsamples=100)
model.plot_spatial_corr(0.2)
model.plot_corr_matrix()
model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)



model = EQ_model('Pedernales','kinematic','h5',(8,10),15,nramp=9,RotAngle=360-99,sampling = True, nsamples=100)

model.plot_spatial_corr(0.2)
model.plot_corr_matrix()
model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)



model = EQ_model('Gorkha','kinematic','h5',(9,18),10,header = {'lon':2,'lat':1,'depth':3,'strike':4,'dip':5},rake=107,sampling=True, nsamples = 100)

model.plot_spatial_corr(0.2)
model.plot_corr_matrix()

model.plot_distribution('mean','Slip',(8,10),'$U (m)$',padding=10)
model.plot_distribution('mean','U_perp',(8,10),'$U_{\perp}(m)$',padding=10)
model.plot_distribution('mean','U_parallel',(8,10),'$U_{||}(m)$',padding=10)
model.plot_distribution('mean','Tr',(8,10),'$Tr(m)$',padding=10)
model.plot_distribution('mean','Vr',(8,10),'$Vr(km/s)$',padding=10)


'''
