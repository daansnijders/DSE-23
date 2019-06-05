# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:27:30 2019

@author: Lisa
"""

from math import * 

#conversion
inch_to_m=0.0254                                                                # [m/in] inches to metres
inchsq_to_msq=0.00064516                                                        # [m^2/in^2] inches squared to metres squared
ft_to_m=0.3048                                                                  # [m/ft] feet to metres
kts_to_ms = 0.514444444                                                         # [m/s/kts] knots to metres per second
ms_to_kts=1/(kts_to_ms)
m_to_ft = 1/(ft_to_m)
kg_to_lbs = 2.20462262 
lbs_to_kg=1/(kg_to_lbs)
nm_to_m=1852
gallons_to_l=3.78541178


#FLIGHT PARAMETERS parameters
H_ft=37000                                                                      # [ft] cruising altitude in feet
H_m=H_ft*ft_to_m                                                                # [m] cruising altitude in metres
M_cruise=0.75                                                                   # [-] cruising mach number
M_x=0.935                                                                       # [-] technology factor

#air constants
gamma=1.4                                                                       # [-] heat capacity ratio
R=287.05                                                                        # [J/kg/K] specific gas constant


T_0=288.15                                                                      # [K] ISA temperature at sea level
p_0=101325                                                                      # [Pa] ISA pressure at sea level
rho_0=1.225                                                                     # [kg/m^3] ISA density at sea level
g=9.80665                                                                       # [m/s^2] gravitational acceleration

mu_37 = 0.0000143226                                                            # [Pa*s] Dynamic Viscosity at 37000ft
mu_sl = 0.00001789                                                              # [Pa*s] Dynamic Viscosity at sea level

if H_m<11000:
    T=T_0-0.0065*(H_m)                                                          # [K] temperature at cruising altitude
    p=p_0*(T/T_0)**(-g/(R*-0.0065))                                             # [Pa] pressure at cruising altitude
    rho= rho_0*(T/T_0)**(-g/(R*-0.0065)-1)                                      # [kg/m^3] density at cruising altitude

elif H_m<20000:
    T1=216.65
    a1=0
    p1=22631.7

    T=T1+a1*(H_m)                                                               # [K] temperature at cruising altitude
    p=p1*e**(-g/(R*T)*(H_m-11000))                                              # [Pa] pressure at cruising altitude
    rho= rho_0*(T1/T_0)**(-g/(R*-0.0065)-1)*p/p1                                # [kg/m^3] density at cruising altitude
    
#get cruise velocities and speed of sound
a=(gamma*R*T)**0.5                                                              # [m/s] speed of sound
V_cruise=M_cruise*a                                                             # [m/s] true airspeed speed cruise



#FUSELAGE
#data on the interior
N_sa = 5                                                                        # [-] Number of seats abreast
N_aisle = 1                                                                     # [-] Number of aisles
aisle_width = 0.51                                                              # [m] width of aisles
s_clearance = 0.02                                                              # [m] clearance between seat and fuselage
seat_pitch = 32*inch_to_m                                                       # [m] seat pitch length
seat_width = 20*inch_to_m                                                       # [m] seat width length
armrest = 2*inch_to_m                                                           # [m] Armrest width
Xfirst = 5                                                                      # [m] Distance between the c.g. of the first seat and the nose

Xfirst = 5.                                                                     # [m] Distance between the c.g. of the first seat and the nose
Xcargo1 = 5.
Xcargo2 = 20.



#average masses on passengers and payload
M_carry_on = 6.1                                                                # [kg] average carry-on luggage weight
M_check_in  = 16.7                                                              # [kg] average check-in luggage weight
M_pax = 83.8 - M_carry_on                                                       # [kg] average passenger weight
rho_lugg = 170                                                                  # [kg/m3] average luggage density
rho_cargo = 160                                                                 # [kg/m3] average cargo density

#constant dimensions in the aircraft
h_max = 0.9                                                                     # [m] distance between the center of the fuselage and the floor
h_floor = 0.2                                                                   # [m] floor heigth
l_cockpit = 4                                                                   # [m] cockpit length

#engine properties
M_engine = 4800                                                                 # [lbs] mass of a single engine
n_engines = 2                                                                   # [-] number of engines used
n_fueltanks  = 4                                                                # [-] number of integral fuel tanks  
d_eng = 2.006                                                                   # [m] diameter of the engine
d_fan = 2.006                                                                   # [m] diameter of the fan
l_eng = 3.184                                                                   # [m] length of the engine
rho_fuel=810                                                                    # [kg/m*3] wikipedia? at 15 degrees
i_e = 0                                                                         # [deg] incidence angle of the engine
# Undercarriage
N_mw = 4                                                                        # [-] number of wheels mlg
N_nw = 2                                                                        # [-] number of wheels nlg
N_struts = 2                                                                    # [-] number of struts used
stroke = 0.3                                                                    # [m] shock absorber stroke

LCN = 45                                                                        # [-] load classification number    
 
#angles clearance
theta = 15                                                                      # [deg] scrape angle
beta = 17                                                                       # [deg] tip-back angle
phi = 5                                                                         # [deg] tip clearance angle
psi = 55                                                                        # [deg] overturn angle

#constants used for Class-II                     
K_fsp = 6.55                                                                    # [lbs/gal] see page 91 from torenbeek V
K_fc = 0.64     
K_h=1.1                                                                             # for movable incidence stabilizers
K_v= 1                                                                              # is for fuselage mounted H tail
K_f=1.08                                                                             # pressurized fuselage

#landing gear retreactable 
K_gr=1                                                                                  #depends on the low wing or high wing
Ag_main=40                                                                      #torenbeek landinggear main and nose (both retractable)
Bg_main=0.16
Cg_main=0.019
Dg_main=1.5E-5

Ag_nose=20
Bg_nose=0.1
Cg_nose=0
Dg_nose=2E-6      

# Structures constants
loc_front_spar = 0.25                                                           # [-] percentage of the chord for front spar location
loc_rear_spar = 0.75                                                            # [-] percentage of the chord for rear spar location
                  
# class 1 estimations                                        
x_c_wcg = 0.205
x_c_oewcg = 0.225
x_f_cg_l_f = 0.435

# V-n loads diagram
V_dive = [276.62626727779644, 276.62626727779644, 276.62626727779644]           # [m/s]all the same with current inputs, but might change, so keep as list
n_max = [2.661047580279859, 2.661047580279859, 2.661047580279859]               # all the same with current inputs, but might change, so keep as list
