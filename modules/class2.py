# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:06:58 2019

@author: Lisa
"""
from modules.class2_struct_defs import *

class Class2_weight:
    
    def __init__(self, m_takeoff, m_empty, loadfactor ):
        self.MTOW=m_takeoff
        self.OEW=m_empty
        self.n_ult=1.5*loadfactor
    def structural_mass(self):
        M_wing=get_wing_mass(self.MTOW)
        return M_wing
#        M_fuselage=get_fuselage_mass()
#        M_nacelle= get_nacelle_mass()
#        M_empennage=get_empennage_mass()
#        M_landinggear=get_landinggear_mass()
#        
#        M_structure=M_wing +M_fuselage+M_nacelle+M_empennage+M_landinggear
#        
#        
    
#    def powerplant_mass():
#        M_engine=get_engine_mass()
#        M_airinduction=get_airinductuin_mass()
#        M_fuelsystem=get_fuelsystem_mass()
#        M_propulsionsystem=get_propulsionsystem_mass()
#         
        
    def fixed_equipment_mass(self):
       M_fixed_equipment=int(self.n_ult*1000)
       return fixed_equipment
        
MTOW=80000
OEW=40000
n= 2.3
config1= Class2_weight(MTOW, OEW, n )      
 
struct_1=config1.structural_mass()
