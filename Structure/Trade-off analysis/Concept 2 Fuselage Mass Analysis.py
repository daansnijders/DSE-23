# -*- coding: utf-8 -*-
"""
Created on Tue May 14 12:01:24 2019

@author: thong

Analysis of the mass contribution of the modules to the fuselage

Assumptions:
1) Cabin shape is a semi cylinder
2) The total overhead luggage mass at module uniformly distributed along the length
3) The lateral position of the luggage weight is 1/3 of cabin width from each 
    side of the cabin
4) Bending moment is calculated by splitting the module frame by half and clamped at the bottom
    (NEEDS RECHECK)
5) Baggage weight is point load distributed equally across the length of module
6) Material used is Aluminium 7075T6. Yield tensile strength used because not
    primary load structure
    
Description:
1) Calculates the maximum bending moment at the base of the module due to the 
    baggage weight.
2) Calculate the thickness of the module based on the bending stress caused by the baggage
3) Calculate the mass of the module by the product of the density and cross-section 
    and length of module.
4) Levelled, symmetric, stationary flight

Input:
1) The mass of luggage per passenger
2) Length of modules (based on the cabin length and the rows of passenger to be configured)
3) The lateral position of the baggage position
    
Output:
1) Thickness of module structure
2) Mass of module structure
"""
import numpy as np
"""
Parameters
Previous update time: - 
Current update time: 13/5/2019, 12:11PM
-------------------------------------------------------------------------------
"""

""" Total mass of overhead luggage at modules [kg]"""
m_luggage = 6.1*30

""" Length of modules [m]"""
l_modules = 25.92-19.44

""" Width of modules [m]"""
w_modules = 3.446
r_modules = w_modules/2

""" Lateral location of point location from the center [m] """
y = w_modules/6

""" Gravitational Acceleration [m/s^2] """
g = 9.81

""" Safety Factor """
safety_factor = 1.5

""" Load Factor """
load_factor = 2.5

""" Tensile Yield Strength [MPa] """
strength = 503
""" Density [kg/m^3]"""
density = 2801

"""
Analysis
-------------------------------------------------------------------------------
"""
angle = np.arcsin(y/r_modules)

# Vertical Reaction Force [N] distributed across the module length
R_y = m_luggage*g
a = w_modules/2-y
c = 2*y
# Reaction Bending Moment [Nm]
R_M = R_y*a*(a+c)/(2*w_modules)


# Thickness required [mm]
t = ((6 * R_M *safety_factor*load_factor)/(l_modules*10**6*strength))**0.5*1000

# Mass of module
mass = (w_modules + np.pi*w_modules/2)*t/1000*l_modules*density

print("Module Mass = ",mass, " kg")
print("Module Thickness = " ,t , "mm ")



