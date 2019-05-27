# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:21:58 2019

@author: Lisa
"""
from math import *
from inputs.constants import *

W_TO=0                                                                          # [lbs] take-off weight

def get_flightcontrolsystem_mass(W_TO):
    M_fc = K_fc*W_TO*kg_to_lbs**(2/3)
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

def get_cargohandling_mass(l_cabin):
    uspace = 3.68 * 0.75 / 3.18
    S_ff = uspace * l_cabin
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