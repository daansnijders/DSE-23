# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:21:58 2019

@author: Lisa
"""
from math import *
from inputs.constants import *
                                                                       # [lbs] take-off weight

def get_flightcontrolsystem_mass(M_TO):
    M_fc = K_fc*(M_TO*kg_to_lbs)**(2/3)
    return M_fc

def get_hydraulic_pneumatic_mass(M_TO):
    return 0.009*M_TO*kg_to_lbs

def get_electricalsystem_mass(d_f_inner, l_cabin):
    V_pax =  0.25*pi*(d_f_inner*m_to_ft)**2*l_cabin*m_to_ft
    M_els=10.8*V_pax**0.7*(1-0.018*V_pax**0.35)
    return M_els 

def get_avioncis_mass(M_TO):
    return 120+20*n_engines + 0.006*M_TO*kg_to_lbs

def get_environmentsystem_mass(l_cabin):
    return 6.75*(l_cabin*m_to_ft)**1.28

def get_oxygensystem_mass(N_pax):
    return 40+2.4*N_pax

def get_apu_mass(M_TO):
    return 0.0085*M_TO*kg_to_lbs

def get_furnish_mass(M_TO, M_fuel):
    return 0.211*(M_TO*kg_to_lbs - M_fuel*kg_to_lbs)**0.91

def get_cargohandling_mass(l_cabin):
    uspace = 3.68 * 0.75 / 3.18 *m_to_ft
    S_ff = uspace * l_cabin *m_to_ft
    return 3*S_ff

def get_operationitems_mass():
    return 0.

def get_flighttestinstrumentation_mass():
    return 0.

def get_paint_mass(M_TO):
    return 0.0045*M_TO*kg_to_lbs

def get_fixedequipment_mass(M_fc,M_hydr,M_els,M_avion,M_environ,M_oxygen,M_apu,M_furnish,M_cargohand,M_operation,M_flightest,M_paint):
    M_fixedequipment=M_fc+M_hydr+M_els+M_avion+M_environ+M_oxygen+M_apu+M_furnish+M_cargohand+M_operation+M_flightest+M_paint
    
    return  M_fixedequipment # should be the total mass of the aformentioned