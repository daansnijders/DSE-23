# -*- coding: utf-8 -*-
"""
Main program of DSE-23
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""




from inputs.concept_1 import *
from inputs.constants import *
from inputs.performance_inputs import *

from modules.Stability.main_stability import *



'first module of the iteration'
iterate_loading=plot_loadingdiagram(Sland,Cl_TO,Cl_clean,Cl_land,V_climb,c,f,sigma, TOP, CD0,1000,7000,100)


'update the values from iterate_loading into concept 1 MANUALLY'
'Update the fuel fractions into concept 1 '
'Update the OEW from class 2'

'do not forget to delete update Canard size in concept 1 and delete the landing gear and empennage sizing'


'perform class2 weight and cg estimation'

'run aero, stability and control, performance and propulsion and structures'

'perform class 2 weigth estimation again'
