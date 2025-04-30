

import sys
import os  
# adding csi functions to the system path


from displacement import Displacement
from curves import Curves
from ensemble import Observable,Ensemble_Displacement,Ensemble_Stress
from stress_change import Stress
# from csi.ensemble_01 import Ensemble_Stress
pwd = os.getcwd()

# ### Computing condensed curves ###



### Tohoku ###


# Tohoku = Observable(pwd,'Tohoku','kinematic',29,(9,24),observable='Rupture_Time',samples='all', new_station=None,xstation=None,ystation=None,strikeOkada =90,factor=None)
# Tohoku.plot_rupture(smax=35,yloc= 0.70,hspace=-0.72)
# Tohoku.save_T0()


### Gorkha ###

Gorkha = Observable(pwd,'Gorkha','kinematic',10,(9,18),observable='RuptureTime',samples=5000, new_station=None,xstation=None,ystation=None,strikeOkada =90,factor=None)
Gorkha.plot_rupture(dt=3,yloc=0.93,smax=15)
Gorkha.save_T0()
Gorkha = Displacement(pwd,'Gorkha','kinematic',10,(9,18),larrow=1,scale=0.05,GF_on=False,samples= 5000)
Gorkha = Stress(pwd,'Gorkha','kinematic',10,(9,18),new_station=True,xstation=Gorkha.xsrc,ystation=Gorkha.ysrc,samples= 5000)
Gorkha = Ensemble_Displacement(pwd,'Gorkha','kinematic',10,(9,18),nparameters = 650,rake=107,samples = 5000) 
Gorkha = Ensemble_Stress(pwd,'Gorkha','kinematic',10,(9,18),nparameters = 650,rake=107,samples =5000,new_station=True,xstation=Gorkha.xsrc,ystation=Gorkha.ysrc) 
Gorkha_displacement = Displacement(pwd,'Gorkha','kinematic',10,(9,18),larrow=1,scale=0.05,GF_on=False,samples=5000)
Curves(pwd,'Gorkha','kinematic',10,(9,18),1,2,Gorkha_displacement.xsrc,Gorkha_displacement.ysrc,Gorkha_displacement.xstn,Gorkha_displacement.ystn,Gorkha_displacement.xsrc,Gorkha_displacement.ysrc,samples =5000)


# ### Iquique ###


# Iquique = Observable(pwd,'Iquique','kinematic',17,(11,12),observable='RuptureTime',samples='all', new_station=None,xstation=None,ystation=None,strikeOkada =90,factor=None)
# Iquique.plot_rupture(dt=8,yloc=0.93,smax=15)
# Iquique.save_T0()



# ### Illapel ###

# Illapel = Observable(pwd,'Illapel','kinematic',18,(10,17),observable='RuptureTime',samples='all', new_station=None,xstation=None,ystation=None,strikeOkada =90,factor=None)
# Illapel.plot_rupture(dt=8,yloc=0.93,smax=15)
# Illapel.save_T0()


# ### Pedernales ###

# Pedernales = Observable(pwd,'Pedernales','kinematic',15,(8,10),observable='RuptureTime',samples='all',new_station=None,xstation=None,ystation=None,strikeOkada =90,factor=None)
# Pedernales.plot_rupture(dt=5,yloc=0.93,smax=10)
# Pedernales.save_T0()


# # Displacement #

# Gorkha = Displacement(pwd,'Gorkha','kinematic',10,(9,18),larrow=1,scale=0.05,GF_on=False)
# # Gorkha = Ensemble_Stress(pwd,'Gorkha','kinematic',10,(9,18),nparameters = 650,rake=107,new_station=True,xstation=Gorkha.xsrc,ystation=Gorkha.ysrc) 

# # Ensemble of realizations #
# # Note:
# # default `samples`=100 
# # Gorkha requires `rake` = 107
# # Pedernales requires `RotAngle` = 360 -99

# Gorkha = Ensemble_Displacement(pwd,'Gorkha','kinematic',10,(9,18),nparameters = 650,rake=107,samples = 'all') 
# #Tohoku = Ensemble_Displacement(pwd,'Tohoku','static',29,(9,24),nparameters = 432)
# #Tohoku = Ensemble_Displacement(pwd,'Tohoku','kinematic',29,(9,24),nparameters = 866)
# #Iquique = Ensemble_Displacement(pwd,'Iquique','kinematic',17,(11,12),nparameters = 533)
# #Illapel = Ensemble_Displacement(pwd,'Illapel','kinematic',18,(10,17),nparameters = 682)
# # Pedernales = Ensemble_Displacement(pwd,'Pedernales','kinematic',15,(8,10),nparameters = 331,RotAngle=360-99)

# #Gorkha = Ensemble_Stress(pwd,'Gorkha','kinematic',10,(9,18),nparameters = 650,rake=107,new_station=True,xstation=Gorkha.xsrc,ystation=Gorkha.ysrc) 
# # Tohoku = Ensemble_Stress(pwd,'Tohoku','static',29,(9,24),nparameters = 432,new_station=True,xstation=Tohoku.xsrc,ystation=Tohoku.ysrc)
# # Tohoku = Ensemble_Stress(pwd,'Tohoku','kinematic',29,(9,24),nparameters = 866,samples = 'all',new_station=True,xstation=Tohoku.xsrc,ystation=Tohoku.ysrc)
# # Iquique = Ensemble_Stress(pwd,'Iquique','kinematic',17,(11,12),nparameters = 533,new_station=True,xstation=Iquique.xsrc,ystation=Iquique.ysrc)
# # Illapel = Ensemble_Stress(pwd,'Illapel','kinematic',18,(10,17),nparameters = 682,new_station=True,xstation=Illapel.xsrc,ystation=Illapel.ysrc)
# # Pedernales = Ensemble_Stress(pwd,'Pedernales','kinematic',15,(8,10),nparameters = 331,RotAngle=360-99,new_station=True,xstation=Pedernales.xsrc,ystation=Pedernales.ysrc)


# # Displacement #

# # Gorkha = Displacement(pwd,'Gorkha','kinematic',10,(9,18),larrow=1,scale=0.08,GF_on=False)
# # Tohoku = Displacement(pwd,'Tohoku','static',29,(9,24),GF_on=False)
# # Tohoku = Displacement(pwd,'Tohoku','kinematic',29,(9,24),GF_on=False)
# # Iquique = Displacement(pwd,'Iquique','kinematic',17,(11,12),larrow=1,scale=0.05,GF_on=False)
# # Illapel = Displacement(pwd,'Illapel','kinematic',18,(10,17),larrow=2,scale=0.075,GF_on=False)
# # Pedernales = Displacement(pwd,'Pedernales','kinematic',15,(8,10),larrow=1,scale=0.075,GF_on=False)

# # # Stress #
# # Note:
# # To initialize Stress require Gorkha to get initialized first to set new station_coordintes

# # Gorkha = Stress(pwd,'Gorkha','kinematic',10,(9,18),new_station=True,xstation=Gorkha.xsrc,ystation=Gorkha.ysrc)
# # Tohoku = Stress(pwd,'Tohoku','static',29,(9,24),new_station=True,xstation=Tohoku.xsrc,ystation=Tohoku.ysrc)
# # Tohoku = Stress(pwd,'Tohoku','kinematic',29,(9,24),new_station=True,xstation=Tohoku.xsrc,ystation=Tohoku.ysrc)
# # Iquique = Stress(pwd,'Iquique','kinematic',17,(11,12),new_station=True,xstation=Iquique.xsrc,ystation=Iquique.ysrc)
# # Illapel = Stress(pwd,'Illapel','kinematic',18,(10,17),new_station=True,xstation=Illapel.xsrc,ystation=Illapel.ysrc)
# # Pedernales = Stress(pwd,'Pedernales','kinematic',15,(8,10),new_station=True,xstation=Pedernales.xsrc,ystation=Pedernales.ysrc)


# # # Testing `curves.py` : summary curves 
# # Gorkha_displacement = Displacement(pwd,'Gorkha','kinematic',10,(9,18),larrow=1,scale=0.05,GF_on=False,samples=100)
# # # # Gorkha_stress = Stress(pwd,'Gorkha','kinematic',10,(9,18),new_station=True,xstation=Gorkha_displacement.xsrc,ystation=Gorkha_displacement.ysrc,samples= 100)
# # Curves(pwd,'Gorkha','kinematic',10,(9,18),1,2,Gorkha_displacement.xsrc,Gorkha_displacement.ysrc,Gorkha_displacement.xstn,Gorkha_displacement.ystn,Gorkha_displacement.xsrc,Gorkha_displacement.ysrc)

# # Iquique_displacement = Displacement(pwd,'Iquique','kinematic',17,(11,12),larrow=1,scale=0.05,GF_on=False,samples = 100)
# # # # Iquique_stress = Stress(pwd,'Iquique','kinematic',17,(11,12),new_station=True,xstation=Iquique_displacement.xsrc,ystation=Iquique_displacement.ysrc,samples=100)
# # Curves(pwd,'Iquique','kinematic',17,(11,12),1,2,Iquique_displacement.xsrc,Iquique_displacement.ysrc,Iquique_displacement.xstn,Iquique_displacement.ystn,Iquique_displacement.xsrc,Iquique_displacement.ysrc)

# # Illapel_displacement = Displacement(pwd,'Illapel','kinematic',18,(10,17),larrow=2,scale=0.075,GF_on=False,samples = 100)
# # # # Illapel_stress = Stress(pwd,'Illapel','kinematic',18,(10,17),new_station=True,xstation=Illapel_displacement.xsrc,ystation=Illapel_displacement.ysrc,samples = 100)
# # Curves(pwd,'Illapel','kinematic',18,(10,17),1,2,Illapel_displacement.xsrc,Illapel_displacement.ysrc,Illapel_displacement.xstn,Illapel_displacement.ystn,Illapel_displacement.xsrc,Illapel_displacement.ysrc)

# # Tohoku_displacement = Displacement(pwd,'Tohoku','static',29,(9,24),GF_on=False,samples= 100)
# # # # Tohoku_stress = Stress(pwd,'Tohoku','static',29,(9,24),new_station=True,xstation=Tohoku_displacement.xsrc,ystation=Tohoku_displacement.ysrc,samples = 100)
# # Curves(pwd,'Tohoku','static',29,(9,24),1,2,Tohoku_displacement.xsrc,Tohoku_displacement.ysrc,Tohoku_displacement.xstn,Tohoku_displacement.ystn,Tohoku_displacement.xsrc,Tohoku_displacement.ysrc)


# # Tohoku_displacement = Displacement(pwd,'Tohoku','kinematic',29,(9,24),GF_on=False,samples=100)
# # # # Tohoku_stress = Stress(pwd,'Tohoku','kinematic',29,(9,24),new_station=True,xstation=Tohoku_displacement.xsrc,ystation=Tohoku_displacement.ysrc,samples = 100)
# # Curves(pwd,'Tohoku','kinematic',29,(9,24),1,2,Tohoku_displacement.xsrc,Tohoku_displacement.ysrc,Tohoku_displacement.xstn,Tohoku_displacement.ystn,Tohoku_displacement.xsrc,Tohoku_displacement.ysrc)


# # Pedernales_displacement = Displacement(pwd,'Pedernales','kinematic',15,(8,10),larrow=1,scale=0.075,GF_on=False,samples=100)
# # # # Pedernales_stress = Stress(pwd,'Pedernales','kinematic',15,(8,10),new_station=True,xstation=Pedernales_displacement.xsrc,ystation=Pedernales_displacement.ysrc,samples = 100)
# # Curves(pwd,'Pedernales','kinematic',15,(8,10),1,2,Pedernales_displacement.xsrc,Pedernales_displacement.ysrc,Pedernales_displacement.xstn,Pedernales_displacement.ystn,Pedernales_displacement.xsrc,Pedernales_displacement.ysrc)
