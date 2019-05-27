# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:06:58 2019

@author: Lisa
"""
from modules.class2_struct_defs import *
from modules.class2_powerplant_defs import *
from modules.class2_fixedequipment_defs import *

from inputs.constants import *
from inputs.performance_inputs import *

class Class2_weight:
    def __init__(self, m_takeoff, m_empty, loadfactor,m_fuel, l_f, N_pax, S_v,S_h,Cr_t,b,S,lambdas,  ):
        self.MTOW=m_takeoff
        self.OEW=m_empty
        self.n_ult=1.5*loadfactor
        self.M_fuel
        self.l_f
        self.
    def structural_mass(self):
        M_wing          =get_wing_mass(M_MZF,b,S,Cr_t,lambda_2_rad,n_ult)
        M_fuselage      =get_fuselage_mass(V_dive, l_h, w_fus, h_fus, S_fus)
        M_nacelle       = get_nacelle_mass(T_req_TO)
        
        M_horizontaltail_mass   =get_horizontaltail_mass(K_h,S_h,V_dive,lambda_h_2_rad)
        M_verticaltail_mass     =get_verticaltail_mass(K_v,S_v,V_dive,lambda_v_2_rad)
        
        M_landinggear_nose      =get_landinggear_mass(K_gr,Ag_nose,Bg_nose,Cg_nose,Dg_nose,M_TO)
        M_landinggear_main      =get_landinggear_mass(K_gr,Ag_main,Bg_main,Cg_main,Dg_main,M_TO)
        M_landinggear           =M_landinggear_nose+M_landinggear_main
        
        M_structure =get_structural_mass(M_wing,M_fuselage,M_nacelle,M_horizontaltail,M_verticaltail,M_landinggear)
        
        
    
    def powerplant_mass():
        M_engine_total          = get_engine_mass(M_engine)
        M_airinduction          = get_airinduction_mass()
        M_fuelsystem            = get_fuelsystem_mass(M_fuel,K_fsp)
        M_propulsionsystem      = get_propulsionsystem_mass(l_f,b)
        
        M_powerplant            = get_totalpowerplant_mass(M_engine_total,M_airinduction,M_fuelsystem, M_propsystem)
         
        
    def fixed_equipment_mass(self):
        M_fc         = get_flightcontrolsystem_mass(M_TO)
        M_hydr       = get_hydraulic_pneumatic_mass(M_TO)
        M_els        = get_electricalsystem_mass(d_f_inner, l_cabin)
        M_avion      = get_avioncis_mass(M_TO)
        M_environ    = get_environmentsystem_mass(l_cabin)
        M_oxygen     = get_oxygensystem_mass(N_pax)
        M_apu        = get_apu_mass(M_TO)
        M_furnish    = get_furnish_mass(M_TO, M_fuel)
        M_cargohand  = get_cargohandling_mass(S_ff)
        M_operations = get_operationitems_mass()
        M_flighttest = get_flighttestinstrumentation_mass()
        M_paint      = get_paint_mass(M_TO)
        
        M_fixedequipment= get_fixedequipment_mass(M_fc,M_hydr,M_els,M_avion,M_environ,M_oxygen,M_apu,M_furnish,M_cargohand,M_operation,M_flighttest,M_paint)
       return M_fixedequipment
   
    def OEW(self):
       return M_structure+M_powerplant+M_fixedequipment
   


    
#test          
MTOW=80000
OEW=40000
n= 2.3
config1= Class2_weight(MTOW, OEW, n )      
struct_1=config1.structural_mass()
