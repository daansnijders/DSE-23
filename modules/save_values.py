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
    
    output_file = open('output.csv' ,  'w')
    
    
    # Selecting variables to save
    output_file.write('N_pax= ' + str(N_pax) + '\n')
    output_file.write('R = ' + str(R) + '\n')
    output_file.write('T_W ' + str(T_W) + '\n')
    output_file.write('W_S ' + str(W_S) + '\n')
    output_file.write('OEW = ' + str(OEW) + '\n')
    output_file.write('MTOW = ' + str(MTOW) + '\n')
    output_file.write('M_ff ' + str(M_ff) + '\n')
    output_file.write('M_fuel ' + str(M_fuel) + '\n')
    output_file.write('T_req ' + str(T_req) + '\n')
    output_file.write('l_cabin' + str(l_cabin) + '\n')
    output_file.write('d_f_inner' + str(d_f_inner) + '\n')
    output_file.write('d_f_inner' + str(d_f_outer) + '\n')
    output_file.write('A' + str(A[0]) + '\n')
    output_file.write('e' + str(e[0]) + '\n')
    output_file.write('S' + str(S) + '\n')
    output_file.write('b' + str(b) + '\n')
    output_file.write('x_cg' + str(x_cg) + '\n')
    output_file.write('y_cg' + str(y_cg) + '\n')
    output_file.write('z_cg' + str(z_cg) + '\n')
    output_file.write('Cl_des'+str(Cl_des) + '\n')
    output_file.write('CD0'+str(CD0) + '\n')
    output_file.write('CDcruise'+str(CDcruise) + '\n')
    
    # Closing file
    output_file.close()
    print('Saved!')
    

# Opening files
exec(open("./output.csv").read())