# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:00:22 2019

@author: Lisa
"""
from inputs.constants import *
import inputs.concept_1 as c1
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

NOx_total_mass=4125 #[g], for the reference take-off and landing cycle 
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





class greenhousegas_emissions(object):
    def __init__(self,performance,configuration):
        self.performance=performance
        self.configuration=configuration-1
        
    def get_CO2_per_passenger_per_km(self):
        
        self.fuel_per_pax=self.performance.fuel_mass_nominal/c1.N_pax[self.configuration]/c1.R[self.configuration]*1000
        self.CO2_reduction=(1-self.fuel_per_pax/0.021)*100
        return self.fuel_per_pax, self.CO2_reduction
    
    def get_CO2_emissions(self):
        return 3.15*self.performance.fuel_mass_nominal
    
    def get_NOx_mass(self):
        self.fuel_flow_list=[self.performance.fuel_flow_take_off, self.performance.fuel_flow_climb, self.performance.fuel_flow_cruise_breguet,self.performance.fuel_flow_descent, self.performance.fuel_flow_landing] #self.performance.fuel_flow_descent
        self.M_fuel_list=[self.performance.fuel_mass_take_off,self.performance.fuel_mass_climb,self.performance.fuel_mass_cruise_breguet,self.performance.fuel_mass_descent, self.performance.fuel_mass_landing] #self.performance.fuel_mass_take_off
        self.EI_NOx= [29.171*self.fuel_flow_list[i]+3.8626 for i in range(len(self.fuel_flow_list))]
        self.M_NOx=[self.EI_NOx[i]*self.M_fuel_list[i] for i in range(len(self.fuel_flow_list))]
        return sum(self.M_NOx)



#def get_Dp_Foo_NOx_specific(NOx_total,T):               #[g/kN]
#    return NOx_total/T
#    
#
#def get_NOx_reduction_CAEP(Dp_Foo_NOx_caep,Dp_Foo_flight):  #should come to 45%
#    return (Dp_Foo_NOx_caep-Dp_Foo_flight)/Dp_Foo_NOx_caep*100
#




