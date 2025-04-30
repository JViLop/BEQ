# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 09:28:22 2024

@author: joanv
"""

import sys
import copy 

# adding csi functions to the system path
sys.path.insert(0, '/home/josevilo/csi/csi')


from okadafull import displacement, stress, strain
import pandas as pd
import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


class Ensemble(Observable):
    def __init__(self,
                        parent_dir,
                        name,
                        model_type,
                        patchsize,
                        shape_geom,
                        observable = 'displacement',
                        samples='100',
                        new_station=False,
                        xstation = None,
                        ystation = None,
                        strikeOkada = 90,
                        factor = 4,
                        RotAngle = None,
                        rake=90,
                        mu=30.0*1e9,
                        nu=0.25):
            
            self.RotAngle = RotAngle
            self.rake = rake
            self.mu = mu 
            self.nu = nu
        
            
            super().__init__(
                     parent_dir,
                     name,
                     model_type,
                     patchsize,
                     shape_geom,
                     observable = observable,
                     samples=samples,
                     new_station=new_station,
                     xstation = xstation ,
                     ystation = ystation,
                     strikeOkada = strikeOkada,
                     factor = factor)
            self.Multiple_Okada_displacement()
            
    def bin_data_reader(self):
        self.Np = self.nrows*self.ncols
        file_folder = self.dir_finder(len(self.list_in)-1,tail_dir='bin_data')
        file_dir = os.path.join(file_folder,f'{self.name}_{self.model_type}_n_{self.samples}.dat')
        self.bin_data = np.fromfile(file_dir,'float').reshape((self.nparameters,self.samples))
    
    def Multiple_Okada_displacement(self):
        self.bin_data_reader()
        self.mean_model_reader()
        self.source()
        self.station()
        
        
        h5f_name =  f'{self.name}_{self.observable}_nsamples_{self.samples}.h5'
        h5f_folder = self.make_dirs(tail_dir='data_ensemble')
        h5f_dir = os.path.join(h5f_folder,h5f_name)   
        h5file = h5py.File(h5f_dir,'w')
        
        displacement_dset = h5file.create_dataset(f'{self.observable}',shape=(self.samples,len(self.xstn)*len(self.ystn)*3))
            
        for i in range(self.bin_data.self.samples):
            Np = self.Np
            model = self.bin_data[:,i] 
            model0 = np.copy(model)
            if self.name =='Gorkha':
                model[:Np],model[Np:2*Np] = model0[Np:2*Np],model0[:Np]
            
            if self.name=='Pedernales':
                ss0 = np.zeros(Np)
                ds0 = np.zeros(Np)
                for p in range(Np):
                    RotAngle2 = self.RotAngle*(np.pi/180.)
                    
                    rotation = np.arctan2(np.tan(self.strike[p]*(np.pi/180)) - np.tan(RotAngle2), 
                                        np.cos(self.dip_rad[p])*(1.+np.tan(RotAngle2)*np.tan(self.strike[p]*(np.pi/180))))

                    # If RotAngle within ]90, 270], change rotation
                    if self.RotAngle > 90. and self.RotAngle<=270.:
                        rotation += np.pi

                    rp = model[p]    # rake-perpendicular
                    ar = model[p+Np] # rake-parallel

                    ss0[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                    ds0[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
                Uperp = ss0
                Uparallel = ds0       
            else:
                rake0 = (self.rake - 90)*(np.pi/180)
                Uperp =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0) 
                Uparallel = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)      

            ss = self.array_formatter(Uperp)
            ds = self.array_formatter(Uparallel)
            ts = np.zeros(self.xsrc_flat.shape)
            Displacement = displacement(
                                self.xstn_flat,
                                self.ystn_flat,
                                self.zstn_flat, 
                                self.xsrc_flat, 
                                self.ysrc_flat, 
                                self.zsrc_flat, 
                                self.width, 
                                self.length, 
                                self.strike_rad, 
                                self.dip_rad, 
                                ss, 
                                ds, 
                                ts, nu=self.nu)
            
            displacement_dset[i,:] = Displacement.flatten(order='F')
            
        
        h5file.close()
        
                                      
        
        
def dir_creator(parent_dir,name,model_type,parameter,nsamples):
    
    output_dir = os.path.join(parent_dir,'ouput')
    EQ_dir = os.path.join(output_dir,name)
    model_dir =  os.path.join(EQ_dir,'model')
    model_type_dir =  os.path.join(model_dir,model_type)
    sampling_type_dir = os.path.join(model_type_dir,f'{nsamples}'+'_samples')
    parameter_dir = os.path.join(sampling_type_dir ,parameter)
    parameter_data_dir = os.path.join(parameter_dir, 'data')
    parameter_figures_dir = os.path.join(parameter_dir,'figures')
    
    child_dir = {'data':parameter_data_dir,'figures':parameter_figures_dir}
    
    return child_dir

def dir_finder(parent_dir,name,model_type,nsamples):
    input_dir = os.path.join(parent_dir,'input')
    EQ_dir = os.path.join(input_dir,name)
    model_dir =  os.path.join(EQ_dir,'model')
    model_type_dir =  os.path.join(model_dir,model_type)
    sample_type_dir = os.path.join(model_type_dir,f'{nsamples}'+'_samples')
    bin_data_dir = os.path.join(sample_type_dir,'bin_data')
    bin_input_file_dir = os.path.join(bin_data_dir,f'{name}_{model_type}_n_{nsamples}.dat')
    mean_data_dir = os.path.join(sample_type_dir,'mean')
    mean_input_file_dir = os.path.join(mean_data_dir,f'{name}_mean_{model_type}_model.csv')
    outs = {'sample_dir':sample_type_dir,'bin_data':bin_input_file_dir,'mean_data':mean_input_file_dir}
    return outs

def okada_observable_h5_writer(name,
                                        model_type,
                                        patchsize,
                                        shape_geom_file,
                                        nparameters,
                                        nramp=0,
                                        RotAngle = None,
                                        rake=90,
                                        parameter='stress',
                                        nsamples = 10,
                                        sample_type='mean',
                                        custom_strike=False,
                                        strike_value='Default',
                                        mu=30.0*1e9,
                                        nu=0.25):
    
    working_dir = os.getcwd()
    output_dir = os.path.join(working_dir,'output')
    okada_observables_dir = os.path.join(output_dir,'okada_observables_ensemble')
    samples_dir = os.path.join(okada_observables_dir,f'{nsamples}_samples')
    EQ_dir = os.path.join(samples_dir,name)
    if name =='Tohoku':
        model_type_dir = os.path.join(EQ_dir,model_type)
        file_type_dir = os.path.join(model_type_dir,'binary_data')
    else:
        file_type_dir = os.path.join(EQ_dir,'binary_data')
    h5file_dir = os.path.join(file_type_dir,f'{name}_{parameter}_nsamples_{nsamples}.h5')

    try:
        os.makedirs(file_type_dir)
    except FileExistsError:
        pass 

    mean_file_dir = dir_finder(working_dir,name,model_type,nsamples)['mean_data']    
    bin_file = dir_finder(working_dir,name,model_type,nsamples)['bin_data'] 
    bin_out_dir =   dir_finder(working_dir,name,model_type,nsamples)['sample_dir'] 
    bin_data = np.fromfile(bin_file,'float').reshape((nparameters,nsamples))
    # Accesing data with mean model ; need to include sampled model case
    df = pd.read_csv(mean_file_dir)
    df= df.drop(df.columns[0],axis=1)
    npatches = int(shape_geom_file[0]*shape_geom_file[1])
    nrows = shape_geom_file[0] # before for Tohoku = 9 , Iquique = 11, etc
    ncols = shape_geom_file[1] # before for Tohoku = 24, Iquique = 12, etc
    dip = np.flip(df['dip'].values.reshape(ncols,nrows).transpose(),axis=0).flatten()
    depth = np.flip(df['depth'].values.reshape(ncols,nrows).transpose(),axis=0).flatten()
    strike = np.flip(df['strike'].values.reshape(ncols,nrows).transpose(),axis=0).flatten()
    
    try:
        hypocenter = [df['Hypo_as'][0],-df['Hypo_dd'][0]] # negative down-dip component of hypocenter
    except:
        print('No Hypocenter available')
    #input_file_dir = os.path.join(input_data_directory,'n_10000_sampled_models_kinematic_tohoku.dat')
    #data = np.fromfile(input_file_dir,'float').reshape((866,10000))
    
    #-----------------------------------
    # 1. Dislocation and receiver
    #-----------------------------------
    
    patchsize  = patchsize*1e3  # Tohoku =29, Iquique = 17
    
    
    # DISLOCATION
    x = np.arange((1/2)*patchsize,ncols*patchsize,patchsize)
    y = np.arange(-(nrows-1/2)*patchsize,0,patchsize)
    X,Y = np.meshgrid(x,y)
    xc = X.flatten()
    yc = Y.flatten()
    zc  = depth*1e3
    
    length_c = patchsize*np.ones(xc.shape)
    width_c = patchsize*np.ones(yc.shape)
    dip_c = dip*(np.pi)/180
    
    if custom_strike:
        strike_c = strike_value*np.ones(xc.shape)*(np.pi/180)
    else:    
        strike_c = strike*(np.pi)/180 
    
    # STATION
    nobservables = 3 
    h5file = h5py.File(h5file_dir,'w')
    if parameter=='stress':
        xs = xc
        ys = yc
        zs = -zc
        Nxs,Nys = len(x),len(y)
        stress_dset = h5file.create_dataset(f'{parameter}',shape=(nsamples,Nxs*Nys*nobservables))
        x0geomx_stress = h5file.create_dataset('x0_geometry',shape=(1,Nxs))
        y0geom_stress = h5file.create_dataset('y0_geometry',shape=(1,Nys))
        x0geomx_stress[:],y0geom_stress[:] = x,y

        xgeomx_stress = h5file.create_dataset('x_geometry',shape=(1,Nxs))
        ygeom_stress = h5file.create_dataset('y_geometry',shape=(1,Nys))
        xgeomx_stress[:],ygeom_stress[:] = x,y
    else:
        xs = np.arange(patchsize/2 - (ncols//4)*patchsize, ncols*patchsize+(ncols//4)*patchsize,patchsize)
        ys = np.arange(-(nrows//4)*patchsize - (nrows-1/2)*patchsize, (nrows//4)*patchsize,patchsize)
        #xs = np.arange(-(2*ncols//2 + 1)*patchsize/2,(2*ncols//2 + 1)*ncols*patchsize,patchsize)
        #ys = np.arange(-(nrows + 5/2)*patchsize,(1/2)*nrows*patchsize,patchsize) 
        Nxs,Nys = len(xs),len(ys)
        Nx0,Ny0 = len(x),len(y)
        # Here it where changes go
        x0geomx_disp = h5file.create_dataset('x0_geometry',shape=(1,Nx0))
        y0geom_disp = h5file.create_dataset('y0_geometry',shape=(1,Ny0))
        x0geomx_disp[:],y0geom_disp[:] = x,y

        xgeomx_disp = h5file.create_dataset('x_geometry',shape=(1,Nxs))
        ygeom_disp = h5file.create_dataset('y_geometry',shape=(1,Nys))
        xgeomx_disp[:],ygeom_disp[:] = xs,ys
        displacement_dset = h5file.create_dataset(f'{parameter}',shape=(nsamples,Nxs*Nys*nobservables))
        Xs ,Ys = np.meshgrid(xs,ys)
        xs = Xs.flatten()
        ys = Ys.flatten() 
        zs = np.zeros(xs.shape)   # for station, 0-depth array coincident with  seafloor
    
    ts = np.zeros(xc.shape)
    Np,nramp  = ncols*nrows,nramp
    
    if parameter=='stress':
        S1 = np.zeros((Nxs*Nys,3))
    else:
        S1 = np.zeros((Nxs*Nys,3))
    S2 = np.zeros(S1.shape)
    for i in range(bin_data.shape[1]):
        model = bin_data[:,i] 
        model0 = np.copy(model)
        if name =='Gorkha':
            model[:Np],model[Np:2*Np] = model0[Np:2*Np],model0[:Np]
        if name=='Pedernales':
            ss0 = np.zeros(Np)
            ds0 = np.zeros(Np)
            for p in range(Np):
                RotAngle2 = RotAngle*((np.pi) / 180.)
                
                rotation = np.arctan2(np.tan(strike[p]) - np.tan(RotAngle2), 
                                    np.cos(dip[p])*(1.+np.tan(RotAngle2)*np.tan(strike[p])))

                # If RotAngle within ]90, 270], change rotation
                if RotAngle > 90. and RotAngle<=270.:
                    rotation += np.pi

                rp = model[p]    # rake-perpendicular
                ar = model[p+Np] # rake-parallel

                ss0[p] = ar*np.cos(rotation) - rp*np.sin(rotation)
                ds0[p] = ar*np.sin(rotation) + rp*np.cos(rotation)
            Uperp = ss0
            Uparallel = ds0       
        else:
            rake0 = (rake - 90)*(np.pi/180)
            Uperp =  model[:Np]*np.cos(rake0)  - model[Np:2*Np]*np.sin(rake0) 
            Uparallel = model[:Np]*np.sin(rake0)  + model[Np:2*Np]*np.cos(rake0)      

        ss = np.flip(Uperp.reshape(ncols,nrows).transpose(),axis=0).flatten()
        ds = np.flip(Uparallel.reshape(ncols,nrows).transpose(),axis=0).flatten()
            # Calling Okada module for parameter calculations 
        if parameter=='stress':
            STRESS = stress(xs, ys,zs, xc, yc, zc, width_c, length_c, strike_c, dip_c, ss, ds, ts, nu=0.25)[0]
            tnn = np.zeros(npatches)
            ts1 = np.zeros(npatches)
            ts2 = np.zeros(npatches)
            #compute stress for each patch (total number: ns)
            strike = strike_c     # *(np.pi/180) if  input for Okada is in degrees directly uncomment this
            dip = dip_c           # *(np.pi/180)  ibid

            Sxx = STRESS[:,0] 
            Sxy = STRESS[:,1] 
            Sxz = STRESS[:,2] 
            Syy = STRESS[:,3] 
            Syz = STRESS[:,4]
            Szz = STRESS[:,5]
            
            for k in range(npatches):
                
                sxx = Sxx[k]*1e-6
                sxy = Sxy[k]*1e-6
                sxz = Sxz[k]*1e-6
                syy = Syy[k]*1e-6
                syz = Syz[k]*1e-6
                szz = Szz[k]*1e-6
                
                # assemble stress tensor
                
                tau = [[sxx,sxy,sxz],[sxy,syy,syz],[sxz,syz,szz]]
                
                #tau = np.array([[S[k,0], S[k,1], S[k,2]],[S[k,1], S[k,3], S[k,4]], [S[k,2], S[k,4], S[k,5]]])
                
                # set cosine/sine constants
                CS, SS = np.cos(strike[k]), np.sin(strike[k])
                CD, SD = np.cos(dip[k]), np.sin(dip[k])
            
                # set direction vectors
                NVec = np.array([SD*CS, -SS*SD, CD]) # fault-normal
                Ss = np.array([SS, CS, 0]) # along-strike before it was  -CS
                Sd = np.cross(NVec,Ss) # see Jolivet stressfield.py strikedip2normal() ; manually np.array([-CD*CS, CD*SS, SD]) # along-dip
                
                # compute traction vector
                trn = np.matmul(tau, NVec) 
                tnn[k] = np.dot(trn, NVec) # fault-normal
                ts1[k] = np.dot(trn, Ss) # along-strike 
                ts2[k] = np.dot(trn, Sd) # along-dip

            StressNS = np.array([tnn,ts1,ts2])
            StressNS = StressNS.flatten()
            stress_dset[i,:] = StressNS     
        else:
            DISPLACEMENT= displacement(xs, ys,zs, xc, yc, zc, width_c, length_c, strike_c, dip_c, ss, ds, ts, nu=0.25)
            #STRAIN = strain(xs, ys,zs, xc, yc, zc, width_c, length_c, strike_c, dip_c, ss, ds, ts, nu=0.25)
            Xd,Yd,Zd = DISPLACEMENT[:,0],DISPLACEMENT[:,1],DISPLACEMENT[:,2]
            
            DisplacementNS = np.array([Xd,Yd,Zd])
            DisplacementNS = DisplacementNS.flatten() 
            displacement_dset[i,:] = DisplacementNS 
            S1 += DISPLACEMENT
            S2 += DISPLACEMENT**2
    h5file.close()
