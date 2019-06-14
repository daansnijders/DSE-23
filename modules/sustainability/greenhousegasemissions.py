# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:00:22 2019

@author: Lisa
"""
from inputs.constants import *
#engine inputs from easa database
#EINOx from emission database, PW1525G, 16PW110  [g/kg]

#EI_NOx_TO=28.1 
#EI_NOx_cruise=21.2
#EI_NOx_app=11.1
#EI_NOx_idle=6.1

#Dp/Foo, the grams emmited during a 
#standard/reference LTO cycle in[g/kN] reference data 
Dp_Foo_NOx_av= 37.6
Dp_Foo_NOx_sigma= 0             #standard deviation
Dp_Foo_NOx_range=[36.5,38.3]    #range
Dp_Foo_NOx_character=39.8       #charachteristic value


#from emmision database
NOx_n_test=3
NOx_n_engines_tested=3

NOx_total_mass=4125 #[g], for the reference take-off and landing cycle (trivial)
fuel_burnt_ref= 299   #[kg]

pressure_ratio=38.67   #ratio of the end of the combustor and the beginning of the combustor (at ISA and take-off conditions)
bypass_ratio=11.05      # other sources claim to be 12?


#check with the requirements 
#co2
CO_2_reduction=30.8             #in percentage
E175_E1_1500=1727 *gallons_to_l *rho_fuel/1000            #fuel used in E1 for 1500 nm
E190_E1_1500=2133 *gallons_to_l *rho_fuel/1000
fuel_target_1500=E190_E1_1500*(1-CO_2_reduction/100) #for 1500 nm fuel consumption
#NOx
Dp_Foo_NOx_caep6=-1.04+2*pressure_ratio
Dp_Foo_NOx_caep8=-9.88+2*pressure_ratio

def get_EI_NOX_fuelflow(fuel_flow):
    EI_NOx=29.171*fuel_flow+3.8626
    return EI_NOx



class greehousegas_emissions():
    def __init__(self,performance,configuration):
        self.performance=performance
        self.configuration=configuration-1
        
    def get_CO2_per_passenger_per_km(self):
        return self.performance.fuel_mass_nominal/N_pax[self.configuration]/R[self.configuration]
    
    def get_CO2_emissions(self):
        return 3.15*self.performance.fuel_mass_nominal
    
#    def get_fuel_flow_per_phase(fuel_flow_TO, fuel_flow_climb, fuel_flow_cruise,fuel_flow_descend, fuel_flow_landing):
#        fuel_flow_list=[fuel_flow_TO, fuel_flow_climb, fuel_flow_cruise,fuel_flow_descend, fuel_flow_landing]
#        return fuel_flow_list
#    def get_total_fuel_per_phase(M_fuel_TO, M_fuel_climb, M_fuel_cruise,M_fuel_descend, M_fuel_landing):
#        M_fuel_list=[M_fuel_TO, M_fuel_climb, M_fuel_cruise,M_fuel_descend, M_fuel_landing]
#        return M_fuel_list
    
def get_NOx_emissions_total(self):  #[g of NOx]
    return fuel_flow*time_in_flight_phase*EI_NOx 

def get_tot_NOx(self):
    EI_NOx=[get_EI_NOX_fuelflow(fuel_flow) for fuel_flow in fuel_flow_list]
    M_NOx=[get_NOx_emissions_total(M_fuel_list[i],EI_NOx[i]) for i in range(len(EI_NOx))]
    return M_NOx
    

def get_Dp_Foo_NOx_specific(NOx_total,T):               #[g/kN]
    return NOx_total/T

def get_fuel_flow_per_phase(fuel_flow_TO, fuel_flow_climb, fuel_flow_cruise,fuel_flow_descend, fuel_flow_landing):
    fuel_flow_list=[fuel_flow_TO, fuel_flow_climb, fuel_flow_cruise,fuel_flow_descend, fuel_flow_landing]
    return fuel_flow_list
def get_total_fuel_per_phase(M_fuel_TO, M_fuel_climb, M_fuel_cruise,M_fuel_descend, M_fuel_landing):
    M_fuel_list=[M_fuel_TO, M_fuel_climb, M_fuel_cruise,M_fuel_descend, M_fuel_landing]
    return M_fuel_list
#
def get_tot_NOx(fuel_flow_list,M_fuel_list):
    EI_NOx=[get_EI_NOX_fuelflow(fuel_flow) for fuel_flow in fuel_flow_list]
    M_NOx=[get_NOx_emissions_total(M_fuel_list[i],EI_NOx[i]) for i in range(len(EI_NOx))]
    return M_NOx
    
    
    
def get_CO2_emissions(M_fuel_burnt):
    return 3.15*M_fuel_burnt

def get_NOx_reduction_CAEP(Dp_Foo_NOx_caep,Dp_Foo_flight):  #should come to 45%
    return (Dp_Foo_NOx_caep-Dp_Foo_flight)/Dp_Foo_NOx_caep*100

def get_requirement_CO2_check(M_fuel_burnt):
    return (M_fuel_burnt-fuel_target_1500)/fuel_target_1500*100

def get_CO2_per_passenger_per_km(M_fuel_burnt,N_pax,R):
    return M_fuel_burnt/N_pax/R


