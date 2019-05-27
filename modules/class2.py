# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:06:58 2019

@author: Lisa
"""
from modules.class2_struct_defs import *
from modules.class2_powerplant_defs import *
from modules.class2_fixedequipment_defs import *

from inputs.concept1 import *
from inputs.constants import *
from inputs.performance_inputs import *

class Class2_weight:
    def __init__(self,N_pax m_takeoff, loadfactor,M_fuel, T_req,l_f,d_f_inner,l_cabin,S, b, S_v,S_h,Cr_t,lambda_2_rad,lambda_h_2_rad, lambda_v_2_rad, S_fus):
        self.MTOW=m_takeoff
        self.n_ult=1.5*loadfactor
        self.M_fuel=M_fuel
        self.T_req_TO=T_req
        self.l_f=l_f
        self.S=S
        self.b=b
        self.S_v=S_v
        self.S_h=S_h
        self.N_pax=N_pax
        self.S_fus=S_fus
        self.lambda_2_rad=lambda_2_rad
        self.lambda_h_2_rad=lambda_h_2_rad
        self.lambda_v_2_rad=lambda_v_2_rad
        self.l_h=l_h
        self.w_fus=w_fus
        self.h_fus=h_fus
        self.d_f_inner=d_f_inner
        self.l_cabin=l_cabin
        
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
        return M_structure
        
    
    def powerplant_mass():
        M_engine_total          = get_engine_mass()
        M_airinduction          = get_airinduction_mass()
        M_fuelsystem            = get_fuelsystem_mass(M_fuel,K_fsp)
        M_propulsionsystem      = get_propulsionsystem_mass(l_f,b)
        
        M_powerplant            = get_totalpowerplant_mass(M_engine_total,M_airinduction,M_fuelsystem, M_propsystem)
        
        return M_powerplant
        
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

config1     = Class2_weight(N_pax[0], m_takeoff[0], loadfactor[0],M_fuel[0], T_req,l_f[0],d_f_inner[0],l_cabin[0],S[0], b[0], S_v[0],S_h[0],Cr_t[0],lambda_2_rad[0],lambda_h_2_rad[0], lambda_v_2_rad[0], S_fus[0] )     
config2     = Class2_weight(N_pax[1],m_takeoff[1], loadfactor[1],M_fuel[1], T_req, l_f[1],d_f_inner[1],l_cabin[1],S[1], b[1], S_v[1],S_h[1],Cr_t[1],lambda_2_rad[1],lambda_h_2_rad[1], lambda_v_2_rad[1], S_fus[1])
config3     = Class2_weight(N_pax [2],m_takeoff[2], loadfactor[2],M_fuel[2], T_req, l_f[2],d_f_inner[2],l_cabin[2],S[2], b[2], S_v[2],S_h[2],Cr_t[2],lambda_2_rad[2],lambda_h_2_rad[2], lambda_v_2_rad[2], S_fus[2])


struct_1=config1.structural_mass()
struct_2=config2.structural_mass()
struct_3=config.structural_mass()

power_1=config1.powerplant_mass()
power_2=config2.powerplant_mass()
power_3=config3.powerplant_mass()

fixedeq_1=config1.fixed_equipment_mass()
fixedeq_2=config2.fixed_equipment_mass()
fixedeq_3=config3.fixed_equipment_mass()

OEW_1=config1.OEW()
OEW_2=config2.OEW()
OEW_3=config3.OEW()

