# -*- coding: utf-8 -*-
"""
Main program of DSE-23
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""




from inputs.concept_1 import *
from modules.initialsizing_fuselage import *
from modules.initialsizing_planform import *
from modules.initialsizing_loading import *


conf_1 = [90,4000]                                                              # [pax,range]
conf_2 = [120,2000]                                                             # [pax,range]
conf_3 = [120,4000]                                                             # [pax,range]

N_pax = [90,120,120]
R = [4000,2000,4000]


CD0, CD0_TO, CD0_land=dragcoefficient(Cfe,Swet_S)

loadingdiagram=plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,c,f,sigma, TOP, CD0,100,7100,100)
