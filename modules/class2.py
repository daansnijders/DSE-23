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
        M_fuselage=get_fuselage_mass()
        M_nacelle= get_nacelle_mass()
        M_empennage=get_empennage_mass()
        M_landinggear=get_landinggear_mass()
        
        M_structure=M_wing +M_fuselage+M_nacelle+M_empennage+M_landinggear
        
        
    
#    def powerplant_mass():
#        M_engine=get_engine_mass()
#        M_airinduction=get_airinduction_mass()
#        M_fuelsystem=get_fuelsystem_mass()
#        M_propulsionsystem=get_propulsionsystem_mass()
#         
        
    def fixed_equipment_mass(self):
        M_fc=get_flightcontrolsystem_mass(M_TO):
       M_fixedequipment=(self.n_ult*1000)
       return M_fixedequipment
   
   def get_flightcontrolsystem_mass(W_TO):
    M_fc = K_fc*M_TO*kg_to_lbs**(2/3)
    return M_fc
def get_hydraulic_pneumatic_mass(W_TO):
    return 0.009*W_TO*kg_to_lbs
def get_electricalsystem_mass(d_f_inner, l_cabin):
    V_pax =  0.25*pi*d_f_inner**2*l_cabin*m_to_ft**2
    M_els=10.8*V_pax**0.7*(1-0.018*V_pax**0.35)
    return M_els 
def get_avioncis_mass(W_TO):
    return 120+20*n_engines + 0.006*W_TO*kg_to_lbs
def get_environmentsystem_mass(l_cabin):
    return 6.75*l_cabin**1.28*m_to_ft

def get_oxygensystem_mass(N_pax):
    return 40+2.4*N_pax

def get_apu_mass(W_TO):
    return 0.0085*W_TO*kg_to_lbs

def get_furnish_mass(W_TO, M_fuel):
    return 0.211*(W_TO*kg_to_lbs - M_fuel*kg_to_lbs)**0.91

def get_cargohandling_mass(S_ff):
    return 3*S_ff

def get_operationitems_mass():
    return 0

def get_flighttestinstrumentation_mass():
    return 0

def get_paint_mass(W_TO):
    return 0.0045*W_TO

def get_fixedequipment_mass(M_fc,M_hydr,M_els,M_avion,M_environ,M_oxygen,M_apu,M_furnish,M_cargohand,M_operation,M_flightest,M_paint):
    M_fixedequipment=M_fc+M_hydr+M_els+M_avion+M_environ+M_oxygen+M_apu+M_furnish+M_cargohand+M_operation+M_flightest+M_paint
    
    return  M_fixedequipment # should be the total mass of the aformentioned 
    
    
        
MTOW=80000
OEW=40000
n= 2.3
config1= Class2_weight(MTOW, OEW, n )      
 
struct_1=config1.structural_mass()
