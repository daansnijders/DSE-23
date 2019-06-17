# -*- coding: utf-8 -*-
"""
Created on Tue May 14 17:20:01 2019

@author: Stijn
"""
# Loading all variables
import inputs.concept_1 as c1
#from modules.main_class2 import 
import output.class2_itegration as massclass2
import inputs.constants as const
import output.connection_departments as detsiz
import csv 

concept=1
# Opening file
safety_check = input('Really want to write and start a new iteration? If so, write "yes" ')
if safety_check == 'yes':
    
    if concept==1:
        output_file = open('output_concept1.csv' ,  'w')
    if concept==2:
        output_file = open('output_concept2.csv' ,  'w')    
    if concept==3:
        output_file = open('output_concept3.csv' ,  'w')
    
    # Selecting variables to save
    output_file.write('MISSION PROFILE' + '\n')
    output_file.write('N_pax = ' + str(c1.N_pax) + '\n')
    output_file.write('R = ' + str(c1.R) + '\n')
   
    output_file.write('WING LOADING/ THRUST LOADING' + '\n')
    output_file.write('T_W = ' + str(c1.T_W) + '\n')
    output_file.write('W_S = ' + str(c1.W_S) + '\n')
    output_file.write('MASSES' + '\n')
    output_file.write('OEW = ' + str(c1.OEW) + '\n')
    output_file.write('MTOW = ' + str(c1.MTOW) + '\n')
    output_file.write('M_ff = ' + str(c1.M_ff) + '\n')
    output_file.write('M_fuel = ' + str(c1.M_fuel) + '\n')
    output_file.write('T_req = ' + str(c1.T_req) + '\n')
    
    output_file.write('FUSELAGE PARAMETERS' + '\n')
    output_file.write('l_cabin = ' + str(c1.l_cabin) + '\n')
    output_file.write('d_f_inner = ' + str(c1.d_f_inner) + '\n')
    output_file.write('d_f_outer = ' + str(c1.d_f_outer) + '\n')
#    output_file.write('l_nose = ' + str(l_nose) + '\n')
 #   output_file.write('l_tailcone = ' + str(l_tailcone) + '\n')
  #  output_file.write('l_tail = ' + str(l_tail) + '\n')
    output_file.write('l_f = ' + str(c1.l_f) + '\n')
#    output_file.write('d_f_outer = ' + str(d_f_outer) + '\n')
#    output_file.write('V_os = ' + str(V_os) + '\n')
#    output_file.write('V_cc = ' + str(V_cc) + '\n')
#    output_file.write('V_carry_on = ' + str(V_carry_on) + '\n')
#    output_file.write('V_check_in = ' + str(V_check_in) + '\n')
    output_file.write('V_cargo_available = ' + str(c1.V_cargo_available) + '\n')
    
    output_file.write('PROPULSION' + '\n')
    output_file.write('T_req = ' + str(c1.T_req) + '\n')
#    output_file.write('fuel_cruise = ' + str(fuel_cruise) + '\n')
    output_file.write('d_fan = ' + str(d_fan) + '\n')
    output_file.write('d_nacel = ' + str(c1.d_nacel) + '\n')
    output_file.write('l_eng = ' + str(l_eng) + '\n')
    output_file.write('l_nacel = ' + str(c1.l_nacel) + '\n')
#    output_file.write('y_eng = ' + str(y_eng) + '\n')
#    output_file.write('d_eng = ' + str(d_eng) + '\n')
#    output_file.write('z_eng = ' + str(z_eng) + '\n')
    
    
    output_file.write('WING PARAMETERS' + '\n')
    output_file.write('A = ' + str(c1.A) + '\n')
    output_file.write('e = ' + str(c1.e) + '\n')
    output_file.write('S = ' + str(c1.S) + '\n')
    output_file.write('b = ' + str(c1.b) + '\n')
    output_file.write('lambda_4_rad = ' + str(c1.lambda_4_rad) + '\n')
    output_file.write('lambda_2_rad = ' + str(c1.lambda_2_rad) + '\n')
    output_file.write('lambda_le_rad = ' + str(c1.lambda_le_rad) + '\n')
    output_file.write('taper_ratio = ' + str(c1.taper_ratio) + '\n')
    output_file.write('Cr = ' + str(c1.Cr) + '\n')
    output_file.write('Ct = ' + str(c1.Ct) + '\n')
    output_file.write('t_c = ' + str(c1.t_c) + '\n')
    output_file.write('MAC = ' + str(c1.MAC) + '\n')
    output_file.write('y_MAC = ' + str(c1.y_MAC) + '\n')
    output_file.write('dihedral_rad = ' + str(c1.dihedral_rad ) + '\n')
    
    output_file.write('HORIZONTAL TAIL PARAMETERS' + '\n')
#    output_file.write('V_h =' + str(V_h) + '\n')
    output_file.write('A_h  =' + str(detsiz.A_h) + '\n')
    output_file.write('taper_ratio_h =' + str(detsiz.taper_ratio_h) + '\n')
    output_file.write('lambda_h_le_rad =' + str(detsiz.lambda_h_le_rad) + '\n')
#    output_file.write('x_le_h =' + str(x_le_h) + '\n')
    output_file.write('S_h=' + str(detsiz.S_h) + '\n')
    output_file.write('b_h =' + str(detsiz.b_h) + '\n')
    output_file.write('Cr_h =' + str(detsiz.Cr_h) + '\n')
    output_file.write('Ct_h =' + str(detsiz.Ct_h) + '\n')
    
    output_file.write('VERTICAL TAIL PARAMETERS' + '\n')
#    output_file.write('V_v =' + str(V_v) + '\n')
    output_file.write('A_v =' + str(detsiz.A_v) + '\n')
    output_file.write('lambda_v_le_rad =' + str(detsiz.lambda_v_le_rad) + '\n')
#   output_file.write('x_le_v =' + str(x_le_v) + '\n')
    output_file.write('S_v  =' + str(detsiz.S_v) + '\n')
    output_file.write('b_v =' + str(detsiz.b_v) + '\n')
    output_file.write('Cr_v =' + str(detsiz.Cr_v) + '\n')
    output_file.write('Ct_v =' + str(detsiz.Ct_v) + '\n')
    
    
    output_file.write('CANARD PARAMETERS CONFIG 2' + '\n')
#    output_file.write('V_v =' + str(V_v) + '\n')
    output_file.write('A_c2 =' + str(detsiz.A_c2) + '\n')
    output_file.write('lambda_c_le_rad2 =' + str(detsiz.lambda_c_le_rad2) + '\n')
#    output_file.write('x_le_v =' + str(x_le_v) + '\n')
    output_file.write('S_c2  =' + str(detsiz.S_c2) + '\n')
    output_file.write('b_c2 =' + str(detsiz.b_c2) + '\n')
    output_file.write('Cr_c2 =' + str(detsiz.Cr_c2) + '\n')
    output_file.write('Ct_c2 =' + str(detsiz.Ct_c2) + '\n')
    
    output_file.write('CANARD PARAMETERS CONFIG 3' + '\n')
#    output_file.write('V_v =' + str(V_v) + '\n')
    output_file.write('A_c3 =' + str(detsiz.A_c3) + '\n')
    output_file.write('lambda_c_le_rad3 =' + str(detsiz.lambda_c_le_rad3) + '\n')
#    output_file.write('x_le_v =' + str(x_le_v) + '\n')
    output_file.write('S_c3  =' + str(detsiz.S_c3) + '\n')
    output_file.write('b_c3 =' + str(detsiz.b_c3) + '\n')
    output_file.write('Cr_c3 =' + str(detsiz.Cr_c3) + '\n')
    output_file.write('Ct_c3 =' + str(detsiz.Ct_c3) + '\n')
    
    

    output_file.write('UNDERCARRIAGE PARAMETERS' + '\n')
    output_file.write('theta = ' + str(const.theta) + '\n')
    output_file.write('beta = ' + str(const.beta) + '\n')
    output_file.write('phi = ' + str(const.phi) + '\n')
    output_file.write('psi = ' + str(const.psi) + '\n')
    output_file.write('P_mw= ' + str(P_mw) + '\n')
    output_file.write('P_nw = ' + str(P_nw) + '\n')
    output_file.write('x_mlg = ' + str(x_mlg) + '\n')
    output_file.write('y_mlg = ' + str(y_mlg) + '\n')
    output_file.write('z_mlg = ' + str(z_mlg) + '\n')
    output_file.write('l_m = ' + str(l_m) + '\n')
    output_file.write('l_n = ' + str(l_n) + '\n')
    output_file.write('x_nlg = ' + str(x_nlg) + '\n')
    output_file.write('x_nlg = ' + str(y_mlg) + '\n')
    output_file.write('z_nlg = ' + str(z_nlg) + '\n')
    
    
    output_file.write('AERODYNAMICS' + '\n')

    
    output_file.write('CLASS MASS 2 ESTIMATION' + '\n')
    output_file.write('config 1 M struct' + str(massclass2.config1_M_structural)+ '\n')
    output_file.write('config 2 M struct' + str(massclass2.config2_M_structural)+ '\n')
    output_file.write('config 3 M struct' + str(massclass2.config3_M_structural)+ '\n')
    output_file.write('config 1 M pwerplant' + str(massclass2.config1_M_powerplant)+ '\n')
    output_file.write('config 2 M pwerplant' + str(massclass2.config2_M_powerplant)+ '\n')
    output_file.write('config 3 M pwerplant' + str(massclass2.config3_M_powerplant)+ '\n')
    output_file.write('config 1 M fixedeq' + str(massclass2.config1_M_fixedeq)+ '\n')
    output_file.write('config 2 M fixedeq' + str(massclass2.config2_M_fixedeq)+ '\n')
    output_file.write('config 3 M fixedeq' + str(massclass2.config3_M_fixedeq)+ '\n')
    output_file.write('OEW config 1 CLASS 2' + str(massclass2.config1_class2_OEW)+ '\n')
    output_file.write('OEW config 2 CLASS 2' + str(massclass2.config2_class2_OEW)+ '\n')
    output_file.write('OEW config 3 CLASS 2' + str(massclass2.config3_class2_OEW)+ '\n')
    
    
    output_file.write('CLASS CG 2 ESTIMATION' + '\n')
    output_file.write('config 1 X CG ' + str(massclass2.config1_cg_x)+ '\n')
    output_file.write('config 2 X CG ' + str(massclass2.config2_cg_x)+ '\n')
    output_file.write('config 3 X CG ' + str(massclass2.config3_cg_x)+ '\n')
    output_file.write('config 1 y CG ' + str(massclass2.config1_cg_y)+ '\n')
    output_file.write('config 2 y CG ' + str(massclass2.config2_cg_y)+ '\n')
    output_file.write('config 3 y CG ' + str(massclass2.config3_cg_y)+ '\n')
    output_file.write('config 1 z CG ' + str(massclass2.config1_cg_z)+ '\n')
    output_file.write('config 2 z CG ' + str(massclass2.config2_cg_z)+ '\n')
    output_file.write('config 3 z CG ' + str(massclass2.config3_cg_z)+ '\n')
      
    # Closing file
    output_file.close()
    print('Saved!')
    

    # Opening files
    if concept==1:
        exec(open("./output_concept1.csv").read())
    
    if concept==2:
        exec(open("./output_concept2.csv").read())
        
    if concept==3:
        exec(open("./output_concept3.csv").read())





