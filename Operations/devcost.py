# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:38:04 2019

@author: ahmadmahmoud
"""
############################IMPORT MASSES##########################
M_wing = 6359
M_verticaltail = 83
M_horizontaltail = 151
M_canard = 1000
M_empennage = M_canard + M_horizontaltail + M_verticaltail
M_fuselage = 9100
M_landinggear = 2700
M_engines_total = 4354
M_airinduction = 0
M_fuelsystem = 0
M_propulsionsystem = 0
M_fc = 0
M_hydr = 0
M_els = 0
M_avion = 0
M_environ = 0
M_oxygen = 0
M_apu = 6481 
M_systems = M_airinduction + M_fuelsystem + M_propulsionsystem + M_fc + M_hydr + M_els + M_avion + M_environ + M_oxygen + M_apu
M_furnish = 0
M_cargohand = 0
M_operation = 9000
M_payloads = M_furnish + M_cargohand + M_operation

########################PRICES PER KG##############################
wing_per_kg = 39090.12
empennage_per_kg = 114984.16
fuselage_per_kg = 70752.87
landinggear_per_kg = 5509.35
installedengines_per_kg = 19160.35
systems_per_kg = 75633.90
payloads_per_kg = 23728.33

############################TOTAL COSTS############################
wing_cost = wing_per_kg * M_wing
empennage_cost = empennage_per_kg * M_empennage
fuselage_cost = fuselage_per_kg * M_fuselage
landinggear_cost = landinggear_per_kg * M_landinggear
installedengines_cost = installedengines_per_kg * M_engines_total
systems_cost = systems_per_kg * M_systems
payloads_cost = payloads_per_kg * M_payloads

##########################ENG COSTS################################
wing_eng = wing_cost * 0.4
empennage_eng = empennage_cost * 0.4
fuselage_eng = fuselage_cost * 0.4
landinggear_eng = landinggear_cost * 0.4
installedengines_eng = installedengines_cost * 0.4
systems_eng = systems_cost * 0.4
payloads_eng = payloads_cost * 0.4

#########################ME COSTS##################################
wing_me = wing_cost * 0.1
empennage_me = empennage_cost * 0.1
fuselage_me = fuselage_cost * 0.1
landinggear_me = landinggear_cost * 0.1
installedengines_me = installedengines_cost * 0.1
systems_me = systems_cost * 0.1
payloads_me = payloads_cost * 0.1

###########################TOOL DESIGN COSTS#######################
wing_tooldes = wing_cost * 0.105
empennage_tooldes = empennage_cost * 0.105
fuselage_tooldes = fuselage_cost * 0.105
landinggear_tooldes = landinggear_cost * 0.105
installedengines_tooldes = installedengines_cost * 0.105
systems_tooldes = systems_cost * 0.105
payloads_tooldes = payloads_cost * 0.105

###########################TOOL FAB COSTS##########################
wing_toolfab = wing_cost * 0.348
empennage_toolfab = empennage_cost * 0.348
fuselage_toolfab = fuselage_cost * 0.348
landinggear_toolfab = landinggear_cost * 0.348
installedengines_toolfab = installedengines_cost * 0.348
systems_toolfab = systems_cost * 0.348
payloads_toolfab = payloads_cost * 0.348

###########################SUPPORT COSTS###########################
wing_support = wing_cost * 0.047
empennage_support = empennage_cost * 0.047
fuselage_support = fuselage_cost * 0.047
landinggear_support = landinggear_cost * 0.047
installedengines_support = installedengines_cost * 0.047
systems_support = systems_cost * 0.047
payloads_support = payloads_cost * 0.047















