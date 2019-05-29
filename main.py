# -*- coding: utf-8 -*-
"""
Main program of DSE-23
Created on Fri May  3 09:45:17 2019

@author: Lisa
"""

from modules.class2_struct_defs import *
from modules.class2_powerplant_defs import *
from modules.class2_fixedequipment_defs import *
from modules.class2 import *

from inputs.concept_1 import *
from inputs.constants import *
from inputs.performance_inputs import *


#update this 
concept = 1

#if concept == 1:
#    from inputs.concept_1 import *
#if concept == 2:
#    from inputs.concept_2 import *
#if concept == 3:
#    from inputs.concept_3 import *



config1     = Class2_weight(N_pax[0],MTOW[0], loadfactor[0],V_dive[0],M_fuel[0], T_req[0], l_f[0],d_f_inner[0],d_f_outer[0],l_cabin[0],l_h[0],S[0], b[0], S_v[0],S_h[0],Cr_t[0],lambda_2_rad[0],lambda_h_2_rad[0], lambda_v_2_rad[0], S_fus[0])     
config2     = Class2_weight(N_pax[1],MTOW[1], loadfactor[1],V_dive[1],M_fuel[1], T_req[1], l_f[1],d_f_inner[1],d_f_outer[1],l_cabin[1],l_h[1],S[1], b[1], S_v[1],S_h[1],Cr_t[1],lambda_2_rad[1],lambda_h_2_rad[1], lambda_v_2_rad[1], S_fus[1])
config3     = Class2_weight(N_pax[2],MTOW[2], loadfactor[2],V_dive[2],M_fuel[2], T_req[2], l_f[2],d_f_inner[2],d_f_outer[2],l_cabin[2],l_h[2],S[2], b[2], S_v[2],S_h[2],Cr_t[2],lambda_2_rad[2],lambda_h_2_rad[2], lambda_v_2_rad[2], S_fus[2])

#get a while loop which results in a percentage difference less than 1 percent
#get structural weight
config1.M_structural=config1.structural_mass()
config2.M_structural=config2.structural_mass()
config3.M_structural=config3.structural_mass()
#get powerplant mass
config1.M_powerplant=config1.powerplant_mass()
config2.M_powerplant=config2.powerplant_mass()
config3.M_powerplant=config3.powerplant_mass()
#get fixed equipment mass
config1.M_fixedeq=config1.fixed_equipment_mass()
config2.M_fixedeq=config2.fixed_equipment_mass()
config3.M_fixedeq=config3.fixed_equipment_mass()
#get OEW total
config1.class2_OEW=config1.OEW(config1.M_structural,config1.M_powerplant,config1.M_fixedeq)
config2.class2_OEW=config2.OEW(config2.M_structural,config2.M_powerplant,config2.M_fixedeq)
config3.class2_OEW=config3.OEW(config3.M_structural,config3.M_powerplant,config3.M_fixedeq)

#update MTOW with class 1 and perform iteration
config1.MTOW_class1=get_MTOW_class1(config1.class2_OEW)
config2.MTOW_class1=get_MTOW_class1(config2.class2_OEW)
config3.MTOW_class1=get_MTOW_class1(config3.class2_OEW)