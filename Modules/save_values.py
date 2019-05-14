# -*- coding: utf-8 -*-
"""
Created on Tue May 14 17:20:01 2019

@author: Stijn
"""
# Loading all variables
from main import *

# Opening file
safety_check = input('Really want to write and start a new iteration? If so, write "yes" ')
if safety_check == 'yes':
    
    output_file = open('output.csv' ,  'w')
    
    # Selecting variables to save
    output_file.write('OEW = ' + str(OEW) + '\n')
    output_file.write('MTOW = ' + str(MTOW) + '\n')
    
    # Closing file
    output_file.close()
    print('Saved!')
    
# Opening files
exec(open("./output.csv").read())