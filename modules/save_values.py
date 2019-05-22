# -*- coding: utf-8 -*-
"""
Created on Tue May 14 17:20:01 2019

@author: Stijn
"""
# Loading all variables
from main import *
import csv 


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
    output_file.write('N_pax = ' + str(N_pax) + '\n')
    output_file.write('R = ' + str(R) + '\n')
    output_file.write('MISSION PROFILE' + '\n')
    output_file.write('WING LOADING/ THRUST LOADING' + '\n')
    output_file.write('T_W = ' + str(T_W) + '\n')
    output_file.write('W_S = ' + str(W_S) + '\n')
    output_file.write('MASSES' + '\n')
    output_file.write('OEW = ' + str(OEW) + '\n')
    output_file.write('MTOW = ' + str(MTOW) + '\n')
    output_file.write('M_ff = ' + str(M_ff) + '\n')
    output_file.write('M_fuel = ' + str(M_fuel) + '\n')
    output_file.write('T_req = ' + str(T_req) + '\n')
    
    output_file.write('FUSELAGE PARAMETERS' + '\n')
    output_file.write('l_cabin = ' + str(l_cabin) + '\n')
    output_file.write('d_f_inner = ' + str(d_f_inner) + '\n')
    output_file.write('d_f_outer = ' + str(d_f_outer) + '\n')
    output_file.write('l_nose = ' + str(l_nose) + '\n')
    output_file.write('l_tailcone = ' + str(l_tailcone) + '\n')
    output_file.write('l_tail = ' + str(l_tail) + '\n')
    output_file.write('l_f = ' + str(l_f) + '\n')
    output_file.write('d_f_outer = ' + str(d_f_outer) + '\n')
    output_file.write('V_os = ' + str(V_os) + '\n')
    output_file.write('V_cc = ' + str(V_cc) + '\n')
    output_file.write('V_carry_on = ' + str(V_carry_on) + '\n')
    output_file.write('V_check_in = ' + str(V_check_in) + '\n')
    output_file.write('V_cargo_available = ' + str(V_cargo_available) + '\n')
    
    output_file.write('PROPULSION' + '\n')
    output_file.write('T_req = ' + str(T_req) + '\n')
    output_file.write('fuel_cruise = ' + str(fuel_cruise) + '\n')
    output_file.write('d_fan = ' + str(d_fan) + '\n')
    output_file.write('d_nacel = ' + str(d_nacel) + '\n')
    output_file.write('l_eng = ' + str(l_eng) + '\n')
    output_file.write('l_nacel = ' + str(l_nacel) + '\n')
    output_file.write('y_eng = ' + str(y_eng) + '\n')
    output_file.write('d_eng = ' + str(d_eng) + '\n')
    output_file.write('z_eng = ' + str(z_eng) + '\n')
    
    
    output_file.write('WING PARAMETERS' + '\n')
    output_file.write('A = ' + str(A) + '\n')
    output_file.write('e = ' + str(e) + '\n')
    output_file.write('S = ' + str(S) + '\n')
    output_file.write('b = ' + str(b) + '\n')
    output_file.write('lambda_4_rad = ' + str(lambda_4_rad) + '\n')
    output_file.write('lambda_2_rad = ' + str(lambda_2_rad) + '\n')
    output_file.write('lambda_le_rad = ' + str(lambda_le_rad) + '\n')
    output_file.write('taper_ratio = ' + str(taper_ratio) + '\n')
    output_file.write('Cr = ' + str(Cr) + '\n')
    output_file.write('Ct = ' + str(Ct) + '\n')
    output_file.write('t_c = ' + str(t_c) + '\n')
    output_file.write('MAC = ' + str(MAC) + '\n')
    output_file.write('y_MAC = ' + str(y_MAC) + '\n')
    output_file.write('dihedral_rad = ' + str(dihedral_rad ) + '\n')
    
    output_file.write('HORIZONTAL TAIL PARAMETERS' + '\n')
    output_file.write('V_h =' + str(V_h) + '\n')
    output_file.write('A_h  =' + str(A_h) + '\n')
    output_file.write('taper_ratio_h =' + str(taper_ratio_h) + '\n')
    output_file.write('lambda_h_le =' + str(lambda_h_le) + '\n')
    output_file.write('x_le_h =' + str(x_le_h) + '\n')
    output_file.write('S_h=' + str(S_h) + '\n')
    output_file.write('b_h =' + str(b_h) + '\n')
    output_file.write('Cr_h =' + str(Cr_h) + '\n')
    output_file.write('Ct_h =' + str(Ct_h) + '\n')
    
    output_file.write('VERTICAL TAIL PARAMETERS' + '\n')
    output_file.write('V_v =' + str(V_v) + '\n')
    output_file.write('A_v =' + str(A_v) + '\n')
    output_file.write('lambda_v_le =' + str(lambda_v_le) + '\n')
    output_file.write('x_le_v =' + str(x_le_v) + '\n')
    output_file.write('S_v  =' + str(S_v) + '\n')
    output_file.write('b_v =' + str(b_h) + '\n')
    output_file.write('Cr_v =' + str(Cr_v) + '\n')
    output_file.write('Ct_v =' + str(Ct_v) + '\n')
    
    output_file.write('CG LOCATIONS' + '\n')
    output_file.write('x_cg =' + str(x_cg) + '\n')
    output_file.write('y_cg = ' + str(y_cg) + '\n')
    output_file.write('z_cg = ' + str(z_cg) + '\n')
    output_file.write('cg_loc='+ str(cg_loc)+'\n')
    
    output_file.write('AERODYNAMIC PARAMETERS' + '\n')
    output_file.write('Cl_des = '+str(Cl_des) + '\n')
    output_file.write('CLdes = '+str(CLdes) + '\n')
    output_file.write('CL_alpha = '+str(CL_alpha) + '\n')
    output_file.write('CLmax = '+str(CLmax) + '\n')
    output_file.write('CLmaxto = '+str(CLmaxto) + '\n')
    output_file.write('CD0 = '+str(CD0) + '\n')
    output_file.write('CDcruise = '+str(CDcruise) + '\n')
    output_file.write('L_over_D = '+str(LoverD) + '\n')
    
    output_file.write('UNDERCARRIAGE PARAMETERS' + '\n')
    output_file.write('theta = ' + str(theta) + '\n')
    output_file.write('beta = ' + str(beta) + '\n')
    output_file.write('phi = ' + str(phi) + '\n')
    output_file.write('psi = ' + str(psi) + '\n')
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
    # Closing file
    output_file.close()
    print('Saved!')
    

# Opening files
exec(open("./output.csv").read())







