# -*- coding: utf-8 -*-
"""
Created on Thu May 16 15:36:07 2019

@author: thong
"""
"""
Calculate the pylon and structure mass of additional fuel tank

Hi Daan!
Can you help me get 
1) the maximum pylon mass of each concept by inputting the maximum mass of fuel
to be carried in the external fuel tank.
2) The structure mass of each concept and configuration by inputting the mass of 
fuel to be carried in the external fuel tank. 

"""

""" Mass of fuel in external fuel tank [kg] 
Concept 1
Concept 2
Concept 3
"""
fuel_mass1 = 4000
fuel_mass2 = 4000
fuel_mass3 = 4000

""" Mass of reference pylon [kg] """
pylon_mass_ref = 61.235

""" Mass of reference fuel tank structure [kg] """
structure_mass_ref = 297.557

""" Mass of reference fuel in external fuel tank [kg] """
fuel_mass_ref = 4131.773
def calc_struc_mass(fuel_mass):
    structure_mass = fuel_mass/fuel_mass_ref*(structure_mass_ref)
    pylon_mass = fuel_mass/fuel_mass_ref*(pylon_mass_ref)
    return [pylon_mass,structure_mass]

pylon_mass1 = calc_struc_mass(fuel_mass1)[0]
pylon_mass2 = calc_struc_mass(fuel_mass2)[0]
pylon_mass3 = calc_struc_mass(fuel_mass3)[0]

structure_mass1 = calc_struc_mass(fuel_mass1)[1]
structure_mass2 = calc_struc_mass(fuel_mass2)[1]
structure_mass3 = calc_struc_mass(fuel_mass3)[1]

print("Fuel Tank Structure Mass Concept 1 = ", structure_mass1, " kg")
print("Fuel Tank Structure Mass Concept 2 = ", structure_mass2, " kg")
print("Fuel Tank Structure Mass Concept 3 = ", structure_mass3, " kg")

print("Fuel Tank Pylon Mass Concept 1 = ", pylon_mass1, " kg")
print("Fuel Tank Pylon Mass Concept 2 = ", pylon_mass2, " kg")
print("Fuel Tank Pylon Mass Concept 3 = ", pylon_mass3, " kg")