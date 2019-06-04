# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 15:00:22 2019

@author: Lisa
"""
#engine inputs from easa database
<<<<<<< HEAD
#EINOx from emission database, [g/kg]

EI_NOx_TO=26.0 
EI_NOx_cruise=20.5
EI_NOx_app=10.4
EI_NOx_idle=6.0

#Dp/Foo, the grams emmited during a 
#standard/reference LTO cycle in[g/kN]
Dp_Foo_NOx_av= 36
Dp_Foo_NOx_sigma= 0             #standard deviation
Dp_Foo_NOx_range=[35.6,36.3]    #range
Dp_Foo_NOx_character=41.7       #charachteristic value
=======
#EINOx from emission database, PW1525G [g/kg]

EI_NOx_TO=28.1 
EI_NOx_cruise=21.2
EI_NOx_app=11.1
EI_NOx_idle=6.1

#Dp/Foo, the grams emmited during a 
#standard/reference LTO cycle in[g/kN]
Dp_Foo_NOx_av= 37.6
Dp_Foo_NOx_sigma= 0             #standard deviation
Dp_Foo_NOx_range=[36.5,38.3]    #range
Dp_Foo_NOx_character=39.8       #charachteristic value
>>>>>>> 92eeb8a5ba1d71a2d16437666985a387fa5f1aca


#from emmision database
NOx_n_test=3
<<<<<<< HEAD
NOx_n_engines_tested=1

NOx_total_mass=3909 #[g], for the reference take-off and landing cycle (trivial)
=======
NOx_n_engines_tested=3

NOx_total_mass=4125 #[g], for the reference take-off and landing cycle (trivial)
fuel_burnt_ref= 299   #[kg]

>>>>>>> 92eeb8a5ba1d71a2d16437666985a387fa5f1aca
pressure_ratio=38.67   #ratio of the end of the combustor and the beginning of the combustor (at ISA and take-off conditions)
bypass_ratio=11.05
T_max=108.53            #[kN]


<<<<<<< HEAD


#check with the requirements

CO_2_reduction=30.8
=======
#fuelflow [kg/s]
fuel_flow_TO=0.79
fuel_flow_cruise=0.65
fuel_flow_app=0.23
fuel_flow_idle=0.08

#power setting in percentage 
power_TO=100
power_climb=85
power_app=30
power_idle=7


#check with the requirements 
CO_2_reduction=30.8             #in percentage
E175_E1_1500=1727               #fuel used in E1 for 1500 nm
E190_E1_1500=2133

>>>>>>> 92eeb8a5ba1d71a2d16437666985a387fa5f1aca

Dp_Foo_NOx_caep6=-1.04+2*pressure_ratio
Dp_Foo_NOx_caep8=-9.88+2*pressure_ratio

def get_NOx_emissions_total(fuel_flow,time_in_flight_phase,EI_NOx):  #[g of NOx]
    return fuel_flow*time_in_flight_phase*EI_NOx 

def get_Dp_Foo_NOx_specific(NOx_total,T):
    return NOx_total/T

def get_CO2_emissions(M_fuel_burnt):
    return 3.1*M_fuel_burnt

def get_NOx_reduction(Dp_Foo_NOx_caep,Dp_Foo_flight):
    return (Dp_Foo_NOx_caep-Dp_Foo_flight)/Dp_Foo_NOx_caep*100

werk=get_NOx_reduction(Dp_Foo_NOx_caep6,Dp_Foo_NOx_av)

