# -*- coding: utf-8 -*-
"""
Main program of DSE-23
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""
from inputs.concept_1 import *
from inputs.concept_2 import *
from inputs.concept_3 import *


concept = 1

if concept == 1:
    from inputs.concept_1 import *
if concept == 2:
    from inputs.concept_2 import *
if concept == 3:
    from inputs.concept_3 import *



from modules.initialsizing_loading import *
#CD0, CD0_TO, CD0_land=dragcoefficient(Cfe,Swet_S)
#loadingdiagram=plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,c,f,sigma, TOP, CD0,100,7100,100)




print('Done!')