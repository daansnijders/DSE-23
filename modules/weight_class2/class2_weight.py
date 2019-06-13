# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:06:58 2019

@author: Lisa
"""
from modules.weight_class2.class2_struct_defs import *
from modules.weight_class2.class2_powerplant_defs import *
from modules.weight_class2.class2_fixedequipment_defs import *

from inputs.concept_1 import *
from inputs.constants import *
from inputs.performance_inputs import *

class Class2_weight:
    def __init__(self,config, N_pax, MTOW,M_carried_canard_MZF,M_MZF, n_max,V_dive, M_fuel ,T_req,l_f,d_f_inner,d_f_outer,l_cabin,l_h,S,S_c, b,b_c, S_v,S_h,Cr_t,Cr_t_c,lambda_2_rad,lambda_h_2_rad, lambda_v_2_rad,lambda_c_2_rad, S_fus):
        self.config         = config
        self.M_TO           = MTOW                                              # [kg] 
        self.n_ult          = 1.5*n_max                                         # [-]
        self.V_dive         = V_dive                                            # []?
        self.M_fuel         = M_fuel                                            # []?
        self.T_req_TO       = T_req                                             # [N]?
        self.l_f            = l_f                                               # [m]
        self.l_h            = l_h                                               # []?
        self.S              = S                                                 # [m^2]
        self.b              = b                                                 # [m]
        self.Cr_t           = Cr_t                                              # [m]
        self.S_v            = S_v                                               # [m^2]
        self.S_h            = S_h                                               # [m^2]
        self.S_fus          = S_fus                                             # [m^2]?
        self.N_pax          = N_pax                                             # [-]
        self.lambda_2_rad   = lambda_2_rad                                      # [rad]
        self.lambda_h_2_rad = lambda_h_2_rad                                    # [rad]
        self.lambda_v_2_rad = lambda_v_2_rad                                    # [rad]
        self.d_f_inner      = d_f_inner                                         # [m]
        self.d_f_outer      = d_f_outer                                         # [m]
        self.l_cabin        = l_cabin                                           # [m]
        self.b_c            = b_c                                                # [m]
        self.S_c            = S_c                                                # [m^2]
        self.Cr_t_c         = Cr_t_c                                            # [m]
        self.lambda_c_2_rad = lambda_c_2_rad                                    # [rad]
        self.w_fus          = self.d_f_outer/2                                  # [m]
        self.h_fus          = self.d_f_outer/2                                  # [m]
        self.M_MZF          = M_MZF                                             # [kg]?
        self.M_carried_canard_MZF= M_carried_canard_MZF
        
    def structural_mass(self):
        self.M_wing          =get_wing_mass(self.M_MZF,self.b,self.S,self.Cr_t,self.lambda_2_rad,self.n_ult)*0.95*lbs_to_kg
        self.M_fuselage      =get_fuselage_mass(self.V_dive, self.l_h, self.w_fus, self.h_fus, self.S_fus)*lbs_to_kg
        self.M_nacelle       = get_nacelle_mass(self.T_req_TO)*lbs_to_kg   #TOTAL NACELLE
        if self.config==1:
            self.M_canard    =0
        else:
            self.M_canard    =get_wing_mass(self.M_carried_canard_MZF,self.b_c,self.S_c,self.Cr_t_c,self.lambda_c_2_rad,self.n_ult)*lbs_to_kg
        
        self.M_horizontaltail   =get_horizontaltail_mass(K_h,self.S_h,self.V_dive,self.lambda_h_2_rad)*lbs_to_kg
        self.M_verticaltail     =get_verticaltail_mass(K_v,self.S_v,self.V_dive,self.lambda_v_2_rad)*lbs_to_kg
        #add canard for the configguration 2 and 3 
        
        self.M_landinggear_nose      =get_landinggear_mass(K_gr,Ag_nose,Bg_nose,Cg_nose,Dg_nose,max(MTOW))*lbs_to_kg
        self.M_landinggear_main      =get_landinggear_mass(K_gr,Ag_main,Bg_main,Cg_main,Dg_main,max(MTOW))*lbs_to_kg
        self.M_landinggear           =(self.M_landinggear_nose+self.M_landinggear_main)
        
        M_structure =get_structural_mass(self.M_wing,self.M_fuselage,self.M_nacelle,self.M_horizontaltail,self.M_verticaltail,self.M_landinggear, self.M_canard)
        
        return  M_structure # M_wing*lbs_to_kg, M_canard*lbs_to_kg #,M_wing* lbs_to_kg,M_fuselage* lbs_to_kg,M_nacelle* lbs_to_kg,M_horizontaltail* lbs_to_kg,M_verticaltail* lbs_to_kg,M_landinggear* lbs_to_kg,
    
        
    
    def powerplant_mass(self):
        self.M_engines_total          = get_engine_mass()*lbs_to_kg
        self.M_airinduction          = get_airinduction_mass()* lbs_to_kg
        self.M_fuelsystem            = get_fuelsystem_mass(self.M_fuel,K_fsp)* lbs_to_kg
        self.M_propulsionsystem      = get_propulsionsystem_mass(self.l_f,self.b)* lbs_to_kg
        
        M_powerplant            = get_totalpowerplant_mass(self.M_engines_total,self.M_airinduction,self.M_fuelsystem, self.M_propulsionsystem)
        
        return M_powerplant# M_engine_total*lbs_to_kg, M_airinduction*lbs_to_kg,M_fuelsystem*lbs_to_kg, M_propulsionsystem*lbs_to_kg
    
    def get_wing_group_mass(self):
        M_wing_group=self.M_wing+self.M_engines_total+self.M_landinggear_main+self.M_nacelle
        return M_wing_group    
        
    def get_fuselage_group_mass(self):
        M_fuselage_group=self.M_canard+self.M_horizontaltail+self.M_verticaltail+self.M_fuselage+self.M_landinggear_nose
        return M_fuselage_group
        
    def fixed_equipment_mass(self):
        M_fc         = get_flightcontrolsystem_mass(self.M_TO)* lbs_to_kg
        M_hydr       = get_hydraulic_pneumatic_mass(self.M_TO)* lbs_to_kg
        M_els        = get_electricalsystem_mass(self.d_f_inner, self.l_cabin)* lbs_to_kg
        M_avion      = get_avioncis_mass(self.M_TO)* lbs_to_kg
        M_environ    = get_environmentsystem_mass(self.l_cabin)* lbs_to_kg
        M_oxygen     = get_oxygensystem_mass(self.N_pax)* lbs_to_kg
        M_apu        = get_apu_mass(self.M_TO)* lbs_to_kg
        M_furnish    = get_furnish_mass(self.M_TO, self.M_fuel)* lbs_to_kg
        M_cargohand  = get_cargohandling_mass(self.l_cabin)* lbs_to_kg
        M_operation  = get_operationitems_mass()* lbs_to_kg
        M_flighttest = get_flighttestinstrumentation_mass()* lbs_to_kg
        M_paint      = get_paint_mass(self.M_TO)* lbs_to_kg
        
        M_fixedequipment= get_fixedequipment_mass(M_fc,M_hydr,M_els,M_avion,M_environ,M_oxygen,M_apu,M_furnish,M_cargohand,M_operation,M_flighttest,M_paint)

        #return  M_fixedequipment* lbs_to_kg
        #return M_fc* lbs_to_kg, M_hydr* lbs_to_kg, M_els* lbs_to_kg, M_avion* lbs_to_kg, M_environ* lbs_to_kg, M_oxygen* lbs_to_kg, M_apu* lbs_to_kg, M_cargohand* lbs_to_kg, M_operation* lbs_to_kg, M_flighttest* lbs_to_kg, M_paint* lbs_to_kg
        #return M_fc* lbs_to_kg+ M_hydr* lbs_to_kg+ M_els* lbs_to_kg+ M_avion* lbs_to_kg+ M_environ* lbs_to_kg+ M_oxygen* lbs_to_kg+ M_apu* lbs_to_kg+ M_cargohand* lbs_to_kg+ M_operation* lbs_to_kg+ M_flighttest* lbs_to_kg+ M_paint* lbs_to_kg,M_fixedequipment* lbs_to_kg 

        return  M_fixedequipment
        #return M_fc, M_hydr, M_els, M_avion, M_environ, M_oxygen, M_apu, M_furnish, M_cargohand, M_operation, M_flighttest, M_paint, M_fixedequipment


   
    def OEW(self,M_structure,M_powerplant,M_fixedequipment):
      return M_structure+M_powerplant+M_fixedequipment
   
    

         
def get_difference_iteration_MTOW(MTOW_old,MTOW_new):
    return (MTOW_new-MTOW_old)/MTOW_old *100

