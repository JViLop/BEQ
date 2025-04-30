import sys
import copy 

from .okadafull import displacement, stress, strain
from .observable import Observable
import pandas as pd
import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import matplotlib.patches as mpatches



class Curves():
    def __init__(self,
                parent_dir,
                name,
                model_type,
                step,
                want_step,
                patchsize,
                shape_geom,
                scomp,
                dcomp,
                xsrc,
                ysrc,
                xstn_d,
                ystn_d,
                xstn_s,
                ystn_s,
                samples=100,
                patchID = None):
        
        self.parent_dir = parent_dir
        self.name = name
        self.model_type =model_type 
        self.step = step 
        self.want_step = want_step 
        self.patchsize = patchsize 
        self.shape_geom = shape_geom
        self.scomp = scomp
        self.dcomp = dcomp

        self.xsrc  = xsrc
        self.ysrc = ysrc

        self.xstn_d = xstn_d
        self.ystn_d = ystn_d

        self.xstn_s = xstn_s
        self.ystn_s = ystn_s
    
        self.patchID = patchID
        self.samples = samples
        if self.want_step:
            self.list_mean_slip = [self.parent_dir,self.name,'model',self.model_type,f'step_{str(self.step)}',str(self.samples) + '_samples','mean']
            self.list_out = [self.parent_dir,self.name,'model',self.model_type,f'step_{str(self.step)}','curves_observables',str(self.samples) + '_samples','figures']  
        else:
            self.list_mean_slip = [self.parent_dir,self.name,'model',self.model_type,str(self.samples) + '_samples','mean']
            self.list_out = [self.parent_dir,self.name,'model',self.model_type,'curves_observables',str(self.samples) + '_samples','figures']  
            
        # self.plot_curves_maxpatch()
        #self.plot_all_curves_maxpatch_nonscaled()
        #self.plot_all_curves_no_strike()
        self.save_curves_no_strike()
        self.plot_all_curves_no_strike_combined()
    def list_dirs(self,observable):
        if self.want_step:
            self.list_observable  = [self.parent_dir,self.name,'model',self.model_type,f'step_{str(self.step)}',observable,str(self.samples) + '_samples','data']
            self.list_covariance = [self.parent_dir,self.name,'model',self.model_type,f'step_{str(self.step)}',observable,str(self.samples) + '_samples','covariance']
        else:
            self.list_observable  = [self.parent_dir,self.name,'model',self.model_type,observable,str(self.samples) + '_samples','data']
            self.list_covariance = [self.parent_dir,self.name,'model',self.model_type,observable,str(self.samples) + '_samples','covariance']
        
        return self.list_observable,self.list_covariance

    def dir_manager(self,level,ldir,data_type,tail_dir='data'):
        ldir = ldir
        ldir[-1]=tail_dir
        
        if level ==0:
            path = os.path.join(self.parent_dir,data_type) 
            return path
        else:
            path = os.path.join(self.dir_manager(level-1,ldir,data_type,tail_dir=tail_dir),ldir[level])
            return path

    

    def dir_finder(self,level,ls,folder = 'INPUT',tail_dir='data'):
        return self.dir_manager(level,ls,folder,tail_dir=tail_dir)

    def dir_creator(self,level,tail_dir='data'):
        return self.dir_manager(level,self.list_out,'OUTPUT',tail_dir=tail_dir)
    
    def make_dirs(self,file='figures'):
        folder_dir = self.dir_creator(len(self.list_out)-1,tail_dir=file)
        os.makedirs(folder_dir,exist_ok=True)
        return folder_dir
    
    def plot_dir(self,extra ='',suffix='png'):
        fname = f'{self.name}_{self.samples}_{self.model_type}_observable_curves_{extra}.{suffix}'
        file_dir = os.path.join(self.make_dirs(),fname)
        return file_dir
    


    def array_formatter(self,array):
        return np.flip(array.reshape(self.shape_geom,order='F'),axis=0).flatten()
        
    
    def mean_model_reader(self):
        file_folder = self.dir_finder(len(self.list_mean_slip)-1,self.list_mean_slip,tail_dir='mean')
        file_dir = os.path.join(file_folder,f'{self.name}_mean_{self.model_type}_model.csv')
        df = pd.read_csv(file_dir)
        df = df.drop(df.columns[0],axis=1)
    
        self.Slip = self.array_formatter(df['Slip'].values)
        self.std_Slip = self.array_formatter(df['std_U'].values)
    def find_slip(self,x,y):
         

        n0cols,n0rows = len(self.xsrc),len(self.ysrc)
        
        max_slip_patch = np.argmin(abs(self.Slip-max(self.Slip)))

        
        col0_max_slip,row0_max_slip = max_slip_patch%n0cols, max_slip_patch//n0cols
        
        
        # displacement is defined on customized geometry overlapping prescribed surface 

        
        ncols,nrows = len(x),len(y)
        extra_cols = (ncols-n0cols)//2
        extra_rows = (nrows-n0rows)//2
        
        self.ncol_target = col0_max_slip + extra_cols 
        self.nrow_target = row0_max_slip + extra_rows
        self.max_slip_patch_in_stn = self.nrow_target*ncols + self.ncol_target 
        

        return self.ncol_target,self.nrow_target, self.max_slip_patch_in_stn
            

    def prepare_files(self):
        ls_stress, ls_cov_stress = self.list_dirs('Stress change')
        ls_displacement, ls_cov_displacement = self.list_dirs('displacement')
        
        cov_stress_folder  = self.dir_finder(len(ls_cov_stress)-1,ls_cov_stress,folder = 'OUTPUT',tail_dir= ls_cov_stress[-1])
        cov_stress_fname = f'covariance_{self.name}_Stress change_nsamples_{self.samples}.txt'
        cov_stress_dir = os.path.join(cov_stress_folder,cov_stress_fname)
        
        cov_displacement_folder  = self.dir_finder(len(ls_cov_displacement)-1,ls_cov_displacement,folder = 'OUTPUT',tail_dir = ls_cov_displacement[-1])
        cov_displacement_fname = f'covariance_{self.name}_displacement_nsamples_{self.samples}.txt'
        cov_displacement_dir = os.path.join(cov_displacement_folder,cov_displacement_fname)

        stress_folder  = self.dir_finder(len(ls_cov_stress)-1,ls_cov_stress,folder = 'OUTPUT')
        stress_fname =  f'{self.name}_{self.samples}_{self.model_type}_model_Stress change.txt'
        stress_dir = os.path.join(stress_folder,stress_fname)
        
        displacement_folder  = self.dir_finder(len(ls_cov_displacement)-1,ls_cov_displacement,folder = 'OUTPUT')
        displacement_fname = f'{self.name}_{self.samples}_{self.model_type}_model_displacement.txt'
        displacement_dir = os.path.join(displacement_folder,displacement_fname)
            
        cov_stress = np.loadtxt(cov_stress_dir)    
        cov_displacement = np.loadtxt(cov_displacement_dir)
        self.std_stress = np.sqrt(np.diagonal(cov_stress).reshape(cov_stress.shape[0]//3,3,order='F'))
        self.std_displacement = np.sqrt(np.diagonal(cov_displacement).reshape(cov_displacement.shape[0]//3,3,order='F'))
        self.stress = np.loadtxt(stress_dir) 
        self.displacement = np.loadtxt(displacement_dir)  

    def prepare_input(self):
        self.prepare_files()
        self.mean_model_reader()
        scomp = self.scomp
        dcomp = self.dcomp

        ncol_target_d,nrow_target_d, max_slip_patch_in_d = self.find_slip(self.xstn_d,self.ystn_d)
        ncol_target_s,nrow_target_s, max_slip_patch_in_s = self.find_slip(self.xstn_s,self.ystn_s)
        ncol_src_target,nrow_src_target, max_slip_patch_in_src = self.find_slip(self.xsrc,self.ysrc)
        
        nrows_src, ncols_src = len(self.ysrc),len(self.xsrc)
        nrows_d,ncols_d = len(self.ystn_d),len(self.xstn_d)
        nrows_s,ncols_s = len(self.ystn_s),len(self.xstn_s)

        self.target_stress = self.stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]
        self.target_displacement = self.displacement[:,dcomp].reshape(nrows_d,ncols_d)[:,ncol_target_d]
        self.max_stress = max(abs(max(self.target_stress)),abs(min(self.target_stress)))
        self.max_displacement = max(abs(max(self.target_displacement)),abs(min(self.target_displacement)))
        
        self.std_up_target_stress = self.target_stress + self.std_stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]
        self.std_down_target_stress = self.target_stress - self.std_stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]

        self.std_up_target_displacement = self.target_displacement + self.std_displacement[:,dcomp].reshape(nrows_d,ncols_d)[:,ncol_target_d]
        self.std_down_target_displacement = self.target_displacement - self.std_displacement[:,dcomp].reshape(nrows_d,ncols_d)[:,ncol_target_d]

        self.mean_slip_target = self.Slip.reshape(nrows_src,ncols_src)[:,ncol_src_target] 
        self.std_up_slip = self.mean_slip_target + self.std_Slip.reshape(nrows_src,ncols_src)[:,ncol_src_target]
        self.std_down_slip = self.mean_slip_target - self.std_Slip.reshape(nrows_src,ncols_src)[:,ncol_src_target]
        self.error_slip =   self.std_Slip.reshape(nrows_src,ncols_src)[:,ncol_src_target]

        self.max_slip = max(self.mean_slip_target)
    
    def plot_curves_maxpatch(self):
        self.prepare_input()
        scomp = self.scomp
        dcomp = self.dcomp
        from pylab import figure,setp
        label_stress = ['normal','along-strike','along-dip']
        label_displacement = ['x','y','z']

        supertitle = f"{self.name} Curves for {label_stress[scomp]} stress and {label_displacement[dcomp]} displacement from {self.samples} samples"

       
        
        fig = figure(dpi=500)
        

        # yprops = dict(rotation=90,
        #             horizontalalignment='right',
        #             verticalalignment='center')
        axprops = dict()
        
        ax1 =fig.add_axes([0.1, 0.60, 0.8, 0.25], **axprops)
        
        ax1.plot(self.ystn_d*1e-3,self.target_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='red')
        ax1.fill_between(self.ystn_d*1e-3,self.std_down_target_displacement,self.std_up_target_displacement,facecolor='pink',alpha=0.5)
        ax1.set_title(supertitle,fontsize=10)
        ax1.axvline(x = 0,ls='--')
        ax1.set_ylabel('Displacement (m)',fontsize=8)
        ax1.tick_params(labelsize=8)
        axprops['sharex'] = ax1
        # force x axes to remain in register, even with toolbar navigation
        ax2 = fig.add_axes([0.1, 0.35, 0.8, 0.25], **axprops)

        
        ax2.plot(self.ystn_s*1e-3,self.target_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='blue')
        ax2.fill_between(self.ystn_s*1e-3,self.std_down_target_stress,self.std_up_target_stress,facecolor='skyblue',alpha=0.5)
        ax2.set_ylabel('Stress Change (Pa)',fontsize=8)
        ax2.axvline(x = 0,ls='--')
        ax2.tick_params(labelsize=8)

        ax3 = fig.add_axes([0.1, 0.1, 0.8, 0.25], **axprops)
        ax3.errorbar(self.ysrc*1e-3,self.mean_slip_target,yerr = self.error_slip,drawstyle='steps-mid',capthick = 0.25,capsize=1.5,marker='.',linewidth=0.75,ms=3, color='black')
        ax3.axvline(x = 0,ls='--')
        ax3.set_ylabel('Slip (m)',fontsize=8)
        ax3.annotate('Trench',(0,0),(2,0))
        ax3.set_xlabel('Trench-perpendicular distance from trench (km)',fontsize=8)
        ax3.tick_params(labelsize=8)
        fig.align_ylabels(fig.get_axes())
        for ax in ax1, ax2:
            setp(ax.get_xticklabels(), visible=False)
        
        
        # fig, ax = plt.subplots(figsize = (7,5),dpi=600)
        # ax.plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='red',label='Displacement')
        # ax.fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='pink',alpha=0.5)
        
        
        
        
        
        # ax.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='blue',label='Stress Change')
        # ax.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='skyblue',alpha=0.5)
        
        

       
    
        # ax.errorbar(self.ysrc*1e-3,self.mean_slip_target/self.max_slip,yerr = self.error_slip/self.max_slip,drawstyle='steps-mid',capthick = 0.25,capsize=2,marker='.',linewidth=0.75,ms=3, color='black',label='Mean Slip')
        
        # # ax.fill_between(self.ysrc*1e-3,self.std_down_slip/self.max_slip,self.std_up_slip/self.max_slip,facecolor='grey',alpha=0.5)

        # ax.axvline(x = 0,ls='--')
        # ax.annotate('Trench',(0,-1),(2,-1))
        # # missing mean slip std ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        
        # ax.set_xlabel('Trench-perpendicular distance from trench (km)')
        # ax.set_ylabel('Normalized')
        # ax.legend()

        # fig, axes = plt.subplots(2,1, figsize = (8,10))
        # line1, =  axes[0].plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Displacement')
        # axes[0].fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='green',alpha=0.5)
        # axes[0].set_xlabel('trench-perpendicular distance from trench')
        # axes[0].set_ylabel('$bar{d}$')
        
        # ax2 = axes[0].twinx() 
        # line2, = ax2.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Stress change')
        # ax2.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='red',alpha=0.5)
        # ax2.set_ylabel('$bar{\sigma} $')
        # # axes[0].legend(loc='upper left')
        # ax2.legend(handles=[line1,line2])
       
    
        # axes[1].plot(self.ysrc*1e-3,self.mean_slip_target,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='mean slip')
        
        # # missing mean slip ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        # axes[1].legend()
        # axes[1].set_xlabel('trench-perpendicular distance from trench')
        # axes[1].set_ylabel('Slip(m)')
        # axes[1].set_title(f'{self.name} Slip along trench-perpendicular ')    
        # fig.suptitle(supertitle,x=0.5,y=1.00,fontweight='bold',fontsize=8)
        # plt.tight_layout()
        fig.savefig(self.plot_dir(extra=f's_{label_stress[scomp]}_d_{label_displacement[dcomp]}'))    
        
        # plt.close()


    def strs(self,scomp):
        self.prepare_files()
        self.mean_model_reader()

        ncol_target_s,nrow_target_s, max_slip_patch_in_s = self.find_slip(self.xstn_s,self.ystn_s)
        ncol_src_target,nrow_src_target, max_slip_patch_in_src = self.find_slip(self.xsrc,self.ysrc)
        
        nrows_src, ncols_src = len(self.ysrc),len(self.xsrc)
        nrows_s,ncols_s = len(self.ystn_s),len(self.xstn_s)

        s = self.stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]
        error_s = self.std_stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]

        out = {'s':s,
                'err_s':error_s
                }
        return out 

    def strs_norm(self,scomp):
        s = self.strs(scomp)['s']
        ds = self.strs(scomp)['err_s']

        max_stress = max(abs(max(s)),abs(min(s)))
        arg_max_stress = np.argmin(abs(abs(s) - max_stress))
    

        s_norm = s/max_stress
        error_s_norm = abs(s_norm)*(np.sqrt((ds/s)**2 + (ds[arg_max_stress]/max_stress)**2))


        out = {'s':s_norm,
                'err_s':error_s_norm
                }
        return out 

    """
    
    def s_norm(self,scomp):
        self.prepare_files()
        self.mean_model_reader()

        ncol_target_s,nrow_target_s, max_slip_patch_in_s = self.find_slip(self.xstn_s,self.ystn_s)
        ncol_src_target,nrow_src_target, max_slip_patch_in_src = self.find_slip(self.xsrc,self.ysrc)
        
        nrows_src, ncols_src = len(self.ysrc),len(self.xsrc)
        nrows_s,ncols_s = len(self.ystn_s),len(self.xstn_s)

        target_stress = self.stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]
        max_stress = max(abs(max(target_stress)),abs(min(target_stress)))
        arg_max_stress = np.argmin(abs(target_stress - max_stress))
    
        
        st_stress = self.std_stress[:,scomp].reshape(nrows_s,ncols_s)[:,ncol_target_s]
        std_up_target_stress, std_down_target_stress= target_stress + st_stress, target_stress - st_stress
    

        stress_norm = target_stress/max_stress
        err_stress = abs(stress_norm)*(np.sqrt((st_stress/target_stress)**2 + (st_stress[arg_max_stress]/max_stress)**2))
        stress_up = stress_norm + err_stress
        stress_down = stress_norm - err_stress

        out = {'s':stress_norm,
                's_up':stress_up,
                's_down': stress_down,
                'err_s':err_stress
                }
        return out 


    """


    def dis(self,dcomp):
        self.prepare_files()
        self.mean_model_reader()
        ncol_target_d,nrow_target_d, max_slip_patch_in_d = self.find_slip(self.xstn_d,self.ystn_d)
        # ncol_src_target,nrow_src_target, max_slip_patch_in_src = self.find_slip(self.xsrc,self.ysrc)
        
        nrows_src, ncols_src = len(self.ysrc),len(self.xsrc)
        nrows_d,ncols_d = len(self.ystn_d),len(self.xstn_d)

        d = self.displacement[:,dcomp].reshape(nrows_d,ncols_d)[:,ncol_target_d]
        error_d =  self.std_displacement[:,dcomp].reshape(nrows_d,ncols_d)[:,ncol_target_d]
        
        out = {'d':d,
                'err_d':error_d
                }

        return out


    def dis_norm(self,dcomp):
        d = self.dis(dcomp)['d']
        dd = self.dis(dcomp)['err_d']

        max_displacement = max(abs(abs(max(d)),abs(min(d))))

        arg_max_displacement = np.argmin(abs(abs(d) - max_displacement))

        d_norm = d/max_displacement
        error_d_norm  = abs(d_norm)*(np.sqrt((dd/d)**2 + (dd[arg_max_displacement]/max_displacement)**2))
        
        out = {'d':d_norm,
                'err_d':error_d_norm
                }

        return out

    def m(self):
        self.prepare_files()
        self.mean_model_reader()
        ncol_src_target,nrow_src_target, max_slip_patch_in_src = self.find_slip(self.xsrc,self.ysrc)
        
        nrows_src, ncols_src = len(self.ysrc),len(self.xsrc)

        mean_slip = self.Slip.reshape(nrows_src,ncols_src)[:,ncol_src_target] 
        error_slip =   self.std_Slip.reshape(nrows_src,ncols_src)[:,ncol_src_target]
        out = {'m':mean_slip,
                'err_m':error_slip
                }
        return out

    
    
    def m_norm(self):
        m = self.m()['m']
        dm = self.m()['err_m']

        max_slip = max(m)
        arg_mix_slip = np.argmin(abs(m - max_slip))

        slip_norm = m/max_slip
        error_slip_norm = slip_norm*(np.sqrt((dm/m)**2+(dm[arg_mix_slip]/max_slip)**2))
        
        out = {'m':slip_norm,
                'err_m':error_slip_norm
                }
        return out

    def plot_all_curves_maxpatch_nonscaled(self):
        
        from pylab import figure,setp
        label_stress = ['normal','along-strike','along-dip']
        label_displacement = ['along-strike','trench-normal','vertical']

        # supertitle = f"{self.name} Curves for {label_stress[scomp]} stress and {label_displacement[dcomp]} displacement from {self.samples} samples"

        color_stress = ['deepskyblue','blue','darkviolet']    
        color_displacement = ['darkred','red','deeppink']  
        lns = ['solid','dashed','dotted']  
        markers = ['o','s','v']
        fig = figure(figsize =(8,6),dpi=900)
        

        # yprops = dict(rotation=90,
        #             horizontalalignment='right',
        #             verticalalignment='center')
        axprops = dict()

        if self.name == 'Gorkha':
            offset = 75 #km
            ystn_d = self.ystn_d - offset*1e3
            ystn_s = self.ystn_s - offset*1e3
            ysrc = self.ysrc - offset*1e3
        else:
            offset = 0 #km
            ystn_d = self.ystn_d - offset
            ystn_s = self.ystn_s - offset
            ysrc = self.ysrc - offset



        ax1 = fig.add_axes([0.1, 0.60, 0.8, 0.25], **axprops)
        for i in range(0,3):
            # ax1.plot(self.ystn_d*1e-3,self.d_norm(i)['d'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_displacement[i],label=label_displacement[i])
            # ax1.fill_between(self.ystn_d*1e-3,self.d_norm(i)['d_down'],self.d_norm(i)['d_up'],facecolor=color_displacement[i],alpha=0.2)
            
            ax1.errorbar(ystn_d*1e-3,self.dis(i)['d'],yerr = self.dis(i)['err_d'],capthick = 0.15,capsize=1.25,linestyle = lns[i],marker=markers[i],markeredgecolor='black',markeredgewidth = 0.2,linewidth=0.8,ms=1, color=color_displacement[i],label=label_displacement[i])
            ax1.set_title(f'{self.name} Slip Model and Observables',fontweight = 'bold',fontsize=13.5)
            if self.name !='Gorkha':
                ax1.axvline(x = 0,ls='--',linewidth=0.6)
            
            ax1.axhline(y = 0,ls='-',linewidth=0.4,color='black')
            ax1.set_ylabel('Displacement (m)',fontsize=9)
            ax1.tick_params(labelsize=8)
            ax1.legend(fontsize=6,loc='upper right',handlelength= 1.5,markerscale = 1.5)
            axprops['sharex'] = ax1
            # force x axes to remain in register, even with toolbar navigation
        ax2 = fig.add_axes([0.1, 0.35, 0.8, 0.25], **axprops)
        for i in range(0,3):    
            # ax2.plot(self.ystn_s*1e-3,self.s_norm(i)['s'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_stress[i],label=label_stress[i])
            # ax2.fill_between(self.ystn_s*1e-3,self.s_norm(i)['s_down'],self.s_norm(i)['s_up'],facecolor=color_stress[i],alpha=0.2)
            ax2.errorbar(ystn_s*1e-3,self.strs(i)['s'],yerr = self.strs(i)['err_s'],capthick = 0.15,capsize=1.25,linestyle = lns[i],marker=markers[i],markeredgecolor='black',markeredgewidth = 0.2,linewidth=0.8,ms=1, color=color_stress[i],label=label_stress[i])
            ax2.set_ylabel('Stress Change (MPa)',fontsize=9)
            if self.name !='Gorkha':
                ax2.axvline(x = 0,ls='--',linewidth=0.6)
            ax2.axhline(y = 0,ls='-',linewidth=0.4,color='black')
            ax2.legend(fontsize=6,loc='upper right',handlelength= 1.5,markerscale = 1.5 )
            ax2.tick_params(labelsize=9)

        ax3 = fig.add_axes([0.1, 0.1, 0.8, 0.25], **axprops)
        # ax3.errorbar(self.ysrc*1e-3,self.m()['m'],yerr = self.m()['err_m'],drawstyle='steps-mid',capthick = 0.2,capsize=1.25,marker='.',linewidth=0.6,ms=2, color='black')
        y0 = ysrc*1e-3
        dysrcleft = np.abs(y0[0] - y0[1])/2
        dysrcright = np.abs(y0[-1] -y0[-2])/2 
        ysrcleft = np.array([y0[0]-dysrcleft,y0[0]])
        ysrcright = np.array([y0[-1],y0[-1] + dysrcright])
        mleft = np.array([self.m()['m'][0],self.m()['m'][0]])
        mright = np.array([self.m()['m'][-1],self.m()['m'][-1]])


        # Create the step function
        
        ax3.plot(ysrcleft, mleft,linewidth=0.8,color='black')
        ax3.plot(ysrcright, mright,linewidth=0.8,color='black')
        ax3.fill_between(ysrcleft, mleft-self.m()['err_m'][0],mleft+self.m()['err_m'][0],step='pre',facecolor='grey',alpha=0.75)
        ax3.fill_between(ysrcright, mright-self.m()['err_m'][-1],mright+self.m()['err_m'][-1],step='post',facecolor='grey',alpha=0.75)
        ax3.plot(ysrc*1e-3,self.m()['m'],drawstyle='steps-mid',marker = '.',linewidth=0.8,ms=2, color='black')
        ax3.fill_between(ysrc*1e-3,self.m()['m'] + self.m()['err_m'] ,self.m()['m'] - self.m()['err_m'],step='mid',facecolor='grey',alpha=0.75)
        if self.name !='Gorkha':
            ax3.axvline(x = 0,ls='--',linewidth=0.6)
        ax3.axhline(y = 0,ls='-',linewidth=0.4,color='black')
        ax3.set_ylabel('Slip (m) ',fontsize=9)
        if self.name =='Gorkha':
            ax3.text(7.5-offset, 3, "Trench",
            ha="center", va="center", rotation=0, size=8,
            bbox=dict(boxstyle="rarrow,pad=0.3",
                      fc="white", ec="black", lw=0.6))
            # ax3.annotate((0,0.5),(2,0.5),ha="center", va="center", rotation=0, size=4,
            # bbox=dict(boxstyle="rarrow,pad=0.3",fc="blue", ec="black", lw=2))
        else:
            ax3.annotate('Trench',(0,0.5),(2,0.5),fontsize=9,rotation = 90)
        ax3.set_xlabel('Distance from trench (km)',fontsize=9)
        ax3.tick_params(labelsize=8)
        fig.align_ylabels(fig.get_axes())
        for ax in ax1, ax2:
            setp(ax.get_xticklabels(), visible=False)
        
        
        # fig, ax = plt.subplots(figsize = (7,5),dpi=600)
        # ax.plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='red',label='Displacement')
        # ax.fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='pink',alpha=0.5)
        
        
        
        
        
        # ax.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='blue',label='Stress Change')
        # ax.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='skyblue',alpha=0.5)
        
        

       
    
        # ax.errorbar(self.ysrc*1e-3,self.mean_slip_target/self.max_slip,yerr = self.error_slip/self.max_slip,drawstyle='steps-mid',capthick = 0.25,capsize=2,marker='.',linewidth=0.75,ms=3, color='black',label='Mean Slip')
        
        # # ax.fill_between(self.ysrc*1e-3,self.std_down_slip/self.max_slip,self.std_up_slip/self.max_slip,facecolor='grey',alpha=0.5)

        # ax.axvline(x = 0,ls='--')
        # ax.annotate('Trench',(0,-1),(2,-1))
        # # missing mean slip std ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        
        # ax.set_xlabel('Trench-perpendicular distance from trench (km)')
        # ax.set_ylabel('Normalized')
        # ax.legend()

        # fig, axes = plt.subplots(2,1, figsize = (8,10))
        # line1, =  axes[0].plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Displacement')
        # axes[0].fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='green',alpha=0.5)
        # axes[0].set_xlabel('trench-perpendicular distance from trench')
        # axes[0].set_ylabel('$bar{d}$')
        
        # ax2 = axes[0].twinx() 
        # line2, = ax2.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Stress change')
        # ax2.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='red',alpha=0.5)
        # ax2.set_ylabel('$bar{\sigma} $')
        # # axes[0].legend(loc='upper left')
        # ax2.legend(handles=[line1,line2])
       
    
        # axes[1].plot(self.ysrc*1e-3,self.mean_slip_target,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='mean slip')
        
        # # missing mean slip ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        # axes[1].legend()
        # axes[1].set_xlabel('trench-perpendicular distance from trench')
        # axes[1].set_ylabel('Slip(m)')
        # axes[1].set_title(f'{self.name} Slip along trench-perpendicular ')    
        # fig.suptitle(supertitle,x=0.5,y=1.00,fontweight='bold',fontsize=8)
        # plt.tight_layout()
        fig.savefig(self.plot_dir(extra=f's_d_all'))  

    def save_curves_no_strike(self):
        if self.name == 'Gorkha':
            offset = 75 #km
            ystn_d = self.ystn_d - offset*1e3
            ystn_s = self.ystn_s - offset*1e3
            ysrc = self.ysrc - offset*1e3
        else:
            offset = 0 #km
            ystn_d = self.ystn_d - offset
            ystn_s = self.ystn_s - offset
            ysrc = self.ysrc - offset

        
        model = {'ysrc':ysrc,'m':self.m()['m'],'m':self.m()['err_m']} 
        displacements = {'ystn_d':ystn_d,'tn':self.disp(1)['d'],'err_tn':self.disp(1)['err_d'],'v':self.disp(2)['d'],'err_v':self.disp(2)['err_d']} 
        stresses = {'ytn_s':ystn_s,'n':self.strs(0)['s'],'err_n':self.strs(0)['err_s'],'ad':self.strs(2)['s'],'err_ad':self.strs(2)['err_s']}
        df_m = pd.DataFrame(model)
        file_m = f'{self.name}_{self.samples}_{self.model_type}_curve_m.csv'
        file_dir = os.path.join(self.make_dirs(file=str(self.samples) + '_samples'),file_m)
        df_m.to_csv(file_dir)
        
        df_disp = pd.DataFrame(displacements)
        file_disp = f'{self.name}_{self.samples}_{self.model_type}_curve_displacement.csv'
        file_dir = os.path.join(self.make_dirs(file=str(self.samples) + '_samples'),file_disp)
        df_disp.to_csv(file_dir)
        
        df_stress = pd.DataFrame(stresses)
        file_stress = f'{self.name}_{self.samples}_{self.model_type}_curve_stress.csv'
        file_dir = os.path.join(self.make_dirs(file=str(self.samples) + '_samples'),file_stress)
        df_stress.to_csv(file_dir)
    
    def plot_all_curves_no_strike(self):
        
        from pylab import figure,setp
        label_stress = ['normal','along-dip']
        label_displacement = ['trench-normal','vertical']

        # supertitle = f"{self.name} Curves for {label_stress[scomp]} stress and {label_displacement[dcomp]} displacement from {self.samples} samples"

        color_stress = ['deepskyblue','blue']    
        color_displacement = ['red','deeppink']  
        lns = ['solid','dashed']  
        markers = ['o','s']
        fig = figure(figsize =(8,6),dpi=900)
        

        # yprops = dict(rotation=90,
        #             horizontalalignment='right',
        #             verticalalignment='center')
        axprops = dict()

        if self.name == 'Gorkha':
            offset = 75 #km
            ystn_d = self.ystn_d - offset*1e3
            ystn_s = self.ystn_s - offset*1e3
            ysrc = self.ysrc - offset*1e3
        else:
            offset = 0 #km
            ystn_d = self.ystn_d - offset
            ystn_s = self.ystn_s - offset
            ysrc = self.ysrc - offset



        ax1 = fig.add_axes([0.1, 0.60, 0.8, 0.25], **axprops)
        for n,(k,i) in enumerate(zip([0,2],[1,2])):
            # ax1.plot(self.ystn_d*1e-3,self.d_norm(i)['d'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_displacement[i],label=label_displacement[i])
            # ax1.fill_between(self.ystn_d*1e-3,self.d_norm(i)['d_down'],self.d_norm(i)['d_up'],facecolor=color_displacement[i],alpha=0.2)
            
            ax1.errorbar(ystn_d*1e-3,self.dis(i)['d'],yerr = self.dis(i)['err_d'],capthick = 0.15,capsize=1.25,linestyle = lns[n],marker=markers[n],markeredgecolor='black',markeredgewidth = 0.2,linewidth=1.25,ms=2, color=color_displacement[n],label=label_displacement[n])
            ax1.set_title(f'{self.name} Slip Model and Observables',fontweight = 'bold',fontsize=13.5)
            if self.name !='Gorkha':
                ax1.axvline(x = 0,ls='--',linewidth=0.6)
            
            ax1.axhline(y = 0,ls='-',linewidth=0.4,color='black')
            ax1.set_ylabel('Displacement (m)',fontsize=9)
            ax1.tick_params(labelsize=8)
            ax1.legend(fontsize=7.5,loc='upper right',handlelength= 1.5,markerscale = 1.0)
            axprops['sharex'] = ax1
            # force x axes to remain in register, even with toolbar navigation
        ax2 = fig.add_axes([0.1, 0.35, 0.8, 0.25], **axprops)
        for n,(k,i) in enumerate(zip([0,2],[1,2])):   
            # ax2.plot(self.ystn_s*1e-3,self.s_norm(i)['s'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_stress[i],label=label_stress[i])
            # ax2.fill_between(self.ystn_s*1e-3,self.s_norm(i)['s_down'],self.s_norm(i)['s_up'],facecolor=color_stress[i],alpha=0.2)
            ax2.errorbar(ystn_s*1e-3,self.strs(k)['s'],yerr = self.strs(k)['err_s'],capthick = 0.15,capsize=1.25,linestyle = lns[n],marker=markers[n],markeredgecolor='black',markeredgewidth = 0.2,linewidth=1.25,ms=2, color=color_stress[n],label=label_stress[n])
            ax2.set_ylabel('Stress Change (MPa)',fontsize=9)
            if self.name !='Gorkha':
                ax2.axvline(x = 0,ls='--',linewidth=0.6)
            ax2.axhline(y = 0,ls='-',linewidth=0.4,color='black')
            ax2.legend(fontsize=7.5,loc='upper right',handlelength= 1.5,markerscale = 1.0 )
            ax2.tick_params(labelsize=9)

        ax3 = fig.add_axes([0.1, 0.1, 0.8, 0.25], **axprops)
        # ax3.errorbar(self.ysrc*1e-3,self.m()['m'],yerr = self.m()['err_m'],drawstyle='steps-mid',capthick = 0.2,capsize=1.25,marker='.',linewidth=0.6,ms=2, color='black')
        y0 = ysrc*1e-3
        dysrcleft = np.abs(y0[0] - y0[1])/2
        dysrcright = np.abs(y0[-1] -y0[-2])/2 
        ysrcleft = np.array([y0[0]-dysrcleft,y0[0]])
        ysrcright = np.array([y0[-1],y0[-1] + dysrcright])
        mleft = np.array([self.m()['m'][0],self.m()['m'][0]])
        mright = np.array([self.m()['m'][-1],self.m()['m'][-1]])


        # Create the step function
        
        ax3.plot(ysrcleft, mleft,linewidth=0.8,color='black')
        ax3.plot(ysrcright, mright,linewidth=0.8,color='black')
        ax3.fill_between(ysrcleft, mleft-self.m()['err_m'][0],mleft+self.m()['err_m'][0],step='pre',facecolor='grey',alpha=0.75)
        ax3.fill_between(ysrcright, mright-self.m()['err_m'][-1],mright+self.m()['err_m'][-1],step='post',facecolor='grey',alpha=0.75)
        ax3.plot(ysrc*1e-3,self.m()['m'],drawstyle='steps-mid',marker = '.',linewidth=1.25,ms=2.5, color='black')
        ax3.fill_between(ysrc*1e-3,self.m()['m'] + self.m()['err_m'] ,self.m()['m'] - self.m()['err_m'],step='mid',facecolor='grey',alpha=0.75)
        if self.name !='Gorkha':
            ax3.axvline(x = 0,ls='--',linewidth=0.6)
        ax3.axhline(y = 0,ls='-',linewidth=0.4,color='black')
        ax3.set_ylabel('Slip (m) ',fontsize=9)
        if self.name =='Gorkha':
            ax3.text(7.5-offset, 3, "Trench",
            ha="center", va="center", rotation=0, size=8,
            bbox=dict(boxstyle="rarrow,pad=0.3",
                      fc="white", ec="black", lw=0.6))
            # ax3.annotate((0,0.5),(2,0.5),ha="center", va="center", rotation=0, size=4,
            # bbox=dict(boxstyle="rarrow,pad=0.3",fc="blue", ec="black", lw=2))
        else:
            ax3.annotate('Trench',(0,0.5),(2,0.5),fontsize=9,rotation = 90)
        ax3.set_xlabel('Distance from trench (km)',fontsize=9)
        ax3.tick_params(labelsize=8)
        fig.align_ylabels(fig.get_axes())
        for ax in ax1, ax2:
            setp(ax.get_xticklabels(), visible=False)
        
        
        # fig, ax = plt.subplots(figsize = (7,5),dpi=600)
        # ax.plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='red',label='Displacement')
        # ax.fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='pink',alpha=0.5)
        
        
        
        
        
        # ax.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='blue',label='Stress Change')
        # ax.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='skyblue',alpha=0.5)
        
        

       
    
        # ax.errorbar(self.ysrc*1e-3,self.mean_slip_target/self.max_slip,yerr = self.error_slip/self.max_slip,drawstyle='steps-mid',capthick = 0.25,capsize=2,marker='.',linewidth=0.75,ms=3, color='black',label='Mean Slip')
        
        # # ax.fill_between(self.ysrc*1e-3,self.std_down_slip/self.max_slip,self.std_up_slip/self.max_slip,facecolor='grey',alpha=0.5)

        # ax.axvline(x = 0,ls='--')
        # ax.annotate('Trench',(0,-1),(2,-1))
        # # missing mean slip std ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        
        # ax.set_xlabel('Trench-perpendicular distance from trench (km)')
        # ax.set_ylabel('Normalized')
        # ax.legend()

        # fig, axes = plt.subplots(2,1, figsize = (8,10))
        # line1, =  axes[0].plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Displacement')
        # axes[0].fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='green',alpha=0.5)
        # axes[0].set_xlabel('trench-perpendicular distance from trench')
        # axes[0].set_ylabel('$bar{d}$')
        
        # ax2 = axes[0].twinx() 
        # line2, = ax2.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Stress change')
        # ax2.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='red',alpha=0.5)
        # ax2.set_ylabel('$bar{\sigma} $')
        # # axes[0].legend(loc='upper left')
        # ax2.legend(handles=[line1,line2])
       
    
        # axes[1].plot(self.ysrc*1e-3,self.mean_slip_target,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='mean slip')
        
        # # missing mean slip ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        # axes[1].legend()
        # axes[1].set_xlabel('trench-perpendicular distance from trench')
        # axes[1].set_ylabel('Slip(m)')
        # axes[1].set_title(f'{self.name} Slip along trench-perpendicular ')    
        # fig.suptitle(supertitle,x=0.5,y=1.00,fontweight='bold',fontsize=8)
        # plt.tight_layout()
        fig.savefig(self.plot_dir(extra=f's_d_all')) 

    def plot_all_curves_no_strike_combined(self):
        
        from pylab import figure,setp
        label_stress = ['normal','along-dip']
        label_displacement = ['trench-normal','vertical']

        # supertitle = f"{self.name} Curves for {label_stress[scomp]} stress and {label_displacement[dcomp]} displacement from {self.samples} samples"

        color_stress = ['deepskyblue','blue']    
        color_displacement = ['red','deeppink']  
        lns = ['solid','solid']  
        markers = ['o','s']
        if self.name=='Pedernales':
            x = 6.6
        else:
            x = 6
        fig = figure(figsize =(x,6),dpi=900)
        

        # yprops = dict(rotation=90,
        #             horizontalalignment='right',
        #             verticalalignment='center')
        axprops = dict()

        if self.name == 'Gorkha':
            offset = 75 #km
            ystn_d = self.ystn_d - offset*1e3
            ystn_s = self.ystn_s - offset*1e3
            ysrc = self.ysrc - offset*1e3
        else:
            offset = 0 #km
            ystn_d = self.ystn_d - offset
            ystn_s = self.ystn_s - offset
            ysrc = self.ysrc - offset



        ax1 = fig.add_axes([0.1, 0.4, 0.8, 0.3], **axprops)
        for n,(k,i) in enumerate(zip([0,2],[1,2])):
            # ax1.plot(self.ystn_d*1e-3,self.d_norm(i)['d'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_displacement[i],label=label_displacement[i])
            # ax1.fill_between(self.ystn_d*1e-3,self.d_norm(i)['d_down'],self.d_norm(i)['d_up'],facecolor=color_displacement[i],alpha=0.2)
            
            ax1.errorbar(ystn_d*1e-3,self.dis(i)['d'],yerr = self.dis(i)['err_d'],capthick = 0.2,capsize=2,linestyle = lns[n],marker=markers[n],markeredgecolor='black',markeredgewidth = 0.2,linewidth=1.25,ms=2, color=color_displacement[n],label=label_displacement[n])
            #ax1.set_title(f'{self.name} Slip Model and Observables',fontweight = 'bold',fontsize=13.5)
            if self.name !='Gorkha':
                ax1.axvline(x = 0,ls='--',linewidth=0.6)
            
            ax1.axhline(y = 0,ls='--',linewidth=0.4,color='red')
            ax1.set_ylabel('Displacement (m)',fontsize=9)
            ax1.tick_params(axis='y',labelsize=9,labelcolor='red')
            ax1.legend(fontsize=7.5,loc='upper left',handlelength= 1.5,markerscale = 1.0)

            axprops['sharex'] = ax1
            #ax1.tick_params(axis='x',labelsize=0)

            # force x axes to remain in register, even with toolbar navigation
        ax2 = ax1.twinx()
        for n,(k,i) in enumerate(zip([0,2],[1,2])):   
            # ax2.plot(self.ystn_s*1e-3,self.s_norm(i)['s'],'o', ls='-',linewidth=0.75,ms=0.75, color=color_stress[i],label=label_stress[i])
            # ax2.fill_between(self.ystn_s*1e-3,self.s_norm(i)['s_down'],self.s_norm(i)['s_up'],facecolor=color_stress[i],alpha=0.2)
            ax2.errorbar(ystn_s*1e-3,self.strs(k)['s'],yerr = self.strs(k)['err_s'],capthick = 0.2,capsize=2,linestyle = lns[n],marker=markers[n],markeredgecolor='black',markeredgewidth = 0.2,linewidth=1.25,ms=2, color=color_stress[n],label=label_stress[n])
            ax2.set_ylabel('Stress Change (MPa)',fontsize=9)
            if self.name !='Gorkha':
                ax2.axvline(x = 0,ls='--',linewidth=0.6)
            ax2.axhline(y = 0,ls='--',linewidth=0.4,color='blue')
            ax2.legend(fontsize=7.5,loc='upper right',handlelength= 1.5,markerscale = 1.0 )
            ax2.tick_params(axis='y',labelsize=9,labelcolor='blue')
            #ax2.tick_params(axis='x',labelsize=0)

            #ax2.set_xticks([])
         #**axprops
        ax3 = fig.add_axes([0.1, 0.1, 0.8, 0.3],**axprops)
        # ax3.errorbar(self.ysrc*1e-3,self.m()['m'],yerr = self.m()['err_m'],drawstyle='steps-mid',capthick = 0.2,capsize=1.25,marker='.',linewidth=0.6,ms=2, color='black')
        y0 = ysrc*1e-3
        dysrcleft = np.abs(y0[0] - y0[1])/2
        dysrcright = np.abs(y0[-1] -y0[-2])/2 
        ysrcleft = np.array([y0[0]-dysrcleft,y0[0]])
        ysrcright = np.array([y0[-1],y0[-1] + dysrcright])
        mleft = np.array([self.m()['m'][0],self.m()['m'][0]])
        mright = np.array([self.m()['m'][-1],self.m()['m'][-1]])


        # Create the step function
        ax3.text(0.025,0.825,f'{self.name}',fontweight='bold',transform=ax3.transAxes)
        ax3.plot(ysrcleft, mleft,linewidth=0.8,color='black')
        ax3.plot(ysrcright, mright,linewidth=0.8,color='black')
        ax3.fill_between(ysrcleft, mleft-self.m()['err_m'][0],mleft+self.m()['err_m'][0],step='pre',facecolor='grey',alpha=0.75)
        ax3.fill_between(ysrcright, mright-self.m()['err_m'][-1],mright+self.m()['err_m'][-1],step='post',facecolor='grey',alpha=0.75)
        ax3.plot(ysrc*1e-3,self.m()['m'],drawstyle='steps-mid',marker = '.',linewidth=1.25,ms=2.5, color='black')
        ax3.fill_between(ysrc*1e-3,self.m()['m'] + self.m()['err_m'] ,self.m()['m'] - self.m()['err_m'],step='mid',facecolor='grey',alpha=0.75)
        if self.name !='Gorkha':
            ax3.axvline(x = 0,ls='--',linewidth=0.6)
        ax3.axhline(y = 0,ls='--',linewidth=0.4,color='black')
        ax3.set_ylabel('Slip (m) ',fontsize=9)
        if self.name =='Gorkha':
            ax3.text(7.5-offset, 3, "Trench",
            ha="center", va="center", rotation=0, size=8,
            bbox=dict(boxstyle="rarrow,pad=0.3",
                      fc="white", ec="black", lw=0.6))
            # ax3.annotate((0,0.5),(2,0.5),ha="center", va="center", rotation=0, size=4,
            # bbox=dict(boxstyle="rarrow,pad=0.3",fc="blue", ec="black", lw=2))
        else:
            ax3.annotate('Trench',(0,0.5),(2,0.5),fontsize=9,rotation = 90)
        ax3.set_xlabel('Distance from trench (km)',fontsize=9)

        ax3.tick_params(labelsize=9)
        fig.align_ylabels(fig.get_axes())
        for ax in ax1, ax2:
            setp(ax.get_xticklabels(), visible=False)
        
        
        # fig, ax = plt.subplots(figsize = (7,5),dpi=600)
        # ax.plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='red',label='Displacement')
        # ax.fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='pink',alpha=0.5)
        
        
        
        
        
        # ax.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='blue',label='Stress Change')
        # ax.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='skyblue',alpha=0.5)
        
        

       
    
        # ax.errorbar(self.ysrc*1e-3,self.mean_slip_target/self.max_slip,yerr = self.error_slip/self.max_slip,drawstyle='steps-mid',capthick = 0.25,capsize=2,marker='.',linewidth=0.75,ms=3, color='black',label='Mean Slip')
        
        # # ax.fill_between(self.ysrc*1e-3,self.std_down_slip/self.max_slip,self.std_up_slip/self.max_slip,facecolor='grey',alpha=0.5)

        # ax.axvline(x = 0,ls='--')
        # ax.annotate('Trench',(0,-1),(2,-1))
        # # missing mean slip std ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        
        # ax.set_xlabel('Trench-perpendicular distance from trench (km)')
        # ax.set_ylabel('Normalized')
        # ax.legend()

        # fig, axes = plt.subplots(2,1, figsize = (8,10))
        # line1, =  axes[0].plot(self.ystn_d*1e-3,self.target_displacement/self.max_displacement,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Displacement')
        # axes[0].fill_between(self.ystn_d*1e-3,self.std_down_target_displacement/self.max_displacement,self.std_up_target_displacement/self.max_displacement,facecolor='green',alpha=0.5)
        # axes[0].set_xlabel('trench-perpendicular distance from trench')
        # axes[0].set_ylabel('$bar{d}$')
        
        # ax2 = axes[0].twinx() 
        # line2, = ax2.plot(self.ystn_s*1e-3,self.target_stress/self.max_stress,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='Stress change')
        # ax2.fill_between(self.ystn_s*1e-3,self.std_down_target_stress/self.max_stress,self.std_up_target_stress/self.max_stress,facecolor='red',alpha=0.5)
        # ax2.set_ylabel('$bar{\sigma} $')
        # # axes[0].legend(loc='upper left')
        # ax2.legend(handles=[line1,line2])
       
    
        # axes[1].plot(self.ysrc*1e-3,self.mean_slip_target,'o', ls='-',linewidth=0.75,ms=0.75, color='black',label='mean slip')
        
        # # missing mean slip ###
        # # axes[1].fill_between(y,std_down_target_obs,std_up_target_obs,facecolor='b',alpha=0.5)
        # # axes[i].plot(y,std_up_target_obs,label='mean $ + \sigma$')
        # # axes[i].plot(y,std_down_target_obs,label='mean $ - \sigma$')
        # axes[1].legend()
        # axes[1].set_xlabel('trench-perpendicular distance from trench')
        # axes[1].set_ylabel('Slip(m)')
        # axes[1].set_title(f'{self.name} Slip along trench-perpendicular ')    
        # fig.suptitle(supertitle,x=0.5,y=1.00,fontweight='bold',fontsize=8)
        # plt.tight_layout()
        fig.savefig(self.plot_dir(extra=f's_d_all_combined')) 
      
      
        






