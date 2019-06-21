# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 13:38:04 2019

@author: ahmadmahmoud
"""
############################IMPORT MASSES##########################
M_wing = 5614.98#
M_verticaltail = 58.97#
M_horizontaltail = 75.73#
M_canard = 626.18#
M_empennage = M_canard + M_horizontaltail + M_verticaltail
M_fuselage = 9604.23#
M_landinggear_nose = 364.07# 
M_landinggear_main = 1967.13#
M_landinggear = M_landinggear_nose + M_landinggear_main
M_engines_total = 4380#
M_airinduction = 0#
M_fuelsystem = 345.10#
M_propulsionsystem = 111.19#
M_fc = 748.20#
M_hydr = 534.15#
M_els = 1583.13#
M_avion = 428.68#
M_environ = 903.40#
M_oxygen = 148.77#
M_apu = 504.48#
M_systems = M_airinduction + M_fuelsystem + M_propulsionsystem + M_fc + M_hydr + M_els + M_avion + M_environ + M_oxygen + M_apu
M_furnish = 3790.83#
M_cargohand = 329.514#
M_operation = 0.0#
M_payloads = M_furnish + M_cargohand + M_operation


########################PRICES PER KG##############################
wing_per_kg = 39090.12 * 1.419106771007
empennage_per_kg = 114984.16 * 1.419106771007
fuselage_per_kg = 70752.87 * 1.419106771007
landinggear_per_kg = 5509.35 * 1.419106771007
installedengines_per_kg = 19160.35 * 1.419106771007
systems_per_kg = 75633.90 * 1.419106771007
payloads_per_kg = 23728.33 * 1.419106771007

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

#==============================================================================
#                               TOTAL COSTS
#==============================================================================
totalcost = wing_cost + empennage_cost+fuselage_cost+landinggear_cost+installedengines_cost +systems_cost +payloads_cost












