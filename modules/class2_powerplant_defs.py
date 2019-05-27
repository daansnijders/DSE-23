# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:20:12 2019

@author: Lisa
"""
from inputs.constants import *

                                                      
                                                    
def get_engine_mass():                                                          #[lbs] Total engine mass (6.1, page 83 torenbeek V)
    M_engine_total = M_engine * n_engines
    return M_engine_total

def get_airinduction_mass(): 
    M_airinduction=0                                                   # [lbs] equal to 0 (6.2, page 83 torenbeek V)
    return M_airinduction

def get_fuelsystem_mass(M_fuel,K_fsp):                                                      # [lbs] mass of the fuel system (6.4, page 83 torenbeek V)
    M_fuelsystem = 80*(n_engines+n_fueltanks-1)+15*n_fueltanks**0.5*(M_fuel*kg_to_lbs/K_fsp)**0.333
    return M_fuelsystem

def get_propulsionsystem_mass(l_f,b):                                           # [lbs] mass of propulsion system (6.5, page 83 torenbeek V)
    M_ec = 88.46*((m_to_ft*l_f + m_to_ft*b)*n_engines/100.)**0.294            # [lbs] engine control        
    M_ess = 9.33 * (get_engine_mass()/1000)**1.078                           # [lbs] engine startup system
    M_propsystem=M_ec+M_ess
    return M_propsystem

def get_totalpowerplant_mass(M_engine_total,M_airinduction,M_fuelsystem, M_propsystem):
    M_powerplant=M_engine_total+M_airinduction+M_fuelsystem+M_propsystem
    return M_powerplant 