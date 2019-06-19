# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:43:30 2019

@author: Lisa
"""

'CLASS 2 WEIGHT ESTIMATION'
import inputs.concept_1 as c1
import modules.weight_class2.class2_weight as c2m
import modules.CG.class2_CG  as c2cg
#import Output.connection_departments as detailedsizing

from inputs.constants import n_max,V_dive
import csv 
import pandas as pd 


output_file = open('output_detailedsizing.dat' ,  'r')
lines= output_file.readlines()
values=[]
labels=[]

#datastored = pd.read_csv('output_detailedsizing.csv' ,  'w'))


for s in lines:
    words=s.split('=')
    
    label=words[0].strip()
    value=words[-1].strip()
    labels.append(label)
    values.append(float(value))

output_file.close()
    
A_h=values[0]  
taper_ratio_h=values[1]
lambda_h_2_rad =values[2]
S_h=values[3]
b_h =values[4]
Cr_h =values[5]
Ct_h =values[6]
A_v =values[7]
lambda_v_2_rad =values[8]
S_v  =values[9]
b_v =values[10]
Cr_v =values[11]
Ct_v =values[12]
A_c2 =values[13]
lambda_c_2_rad2 =values[14]
S_c2  =values[15]
b_c2 =values[16]
Cr_c2 =values[17]
Ct_c2 =values[18]
Cr_t_c2=values[19]
A_c3 =values[20]
lambda_c_2_rad3 =values[21]
S_c3  =values[22]
b_c3 =values[23]
Cr_c3 =values[24]
Ct_c3 =values[25]
Cr_t_c3=values[26]
l_h=values[27]
x_le_MAC1=values[28]
x_le_MAC2=values[29]
x_le_MAC3=values[30]
x_le_MAC=[x_le_MAC1,x_le_MAC2,x_le_MAC3]





'CLASS 2 WEIGHT EXECUTION AND CG LOCATION ESTIMATION'
#configuration 1 
config1_class2     = c2m.Class2_weight(1,c1.N_pax[0],c1.MTOW[0],max(c1.MTOW),c1.M_carried_canard_MZF[0],min(c1.M_MZF), n_max[0],V_dive[0],c1.M_fuel[0], max(c1.T_req), c1.l_f[0],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[0], l_h, c1.S, 0., c1.b, 0., S_v,S_h,c1.Cr_t,0,c1.lambda_2_rad,lambda_h_2_rad, lambda_v_2_rad,0, c1.S_fus[0])     


config1_M_structural            =config1_class2.structural_mass()
config1_M_powerplant            =config1_class2.powerplant_mass()
config1_M_fixedeq               =config1_class2.fixed_equipment_mass()

config1_M_winggroup             =config1_class2.get_wing_group_mass()
config1_M_fuselagegroup         =config1_class2.get_fuselage_group_mass()


config1_class2_OEW              =config1_class2.OEW(config1_M_structural,config1_M_powerplant,config1_M_fixedeq)
#cg locations
config1_cg = c2cg.get_cg(x_le_MAC,config1_class2)   
config1_cg_x=config1_cg.calc_x_cg()
config1_cg_y=config1_cg.calc_y_cg()
config1_cg_z=config1_cg.calc_z_cg()


#configuration 3
#masses
config3_class2     = c2m.Class2_weight(3,c1.N_pax[2],c1.MTOW[2],max(c1.MTOW),c1.M_carried_canard_MZF[2],min(c1.M_MZF), n_max[2],V_dive[2],c1.M_fuel[2], max(c1.T_req), c1.l_f[2],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[2], l_h, c1.S,S_c3, c1.b, b_c3, S_v, S_h,c1.Cr_t, Cr_t_c3, c1.lambda_2_rad,lambda_h_2_rad, lambda_v_2_rad,lambda_c_2_rad3, c1.S_fus[2])

config3_M_structural            =config3_class2.structural_mass()
config3_M_powerplant            =config3_class2.powerplant_mass()
config3_M_fixedeq               =config3_class2.fixed_equipment_mass()

config3_M_winggroup             =config3_class2.get_wing_group_mass()
config3_M_fuselagegroup         =config3_class2.get_fuselage_group_mass()

config3_class2_OEW              =config3_class2.OEW(config3_M_structural,config3_M_powerplant,config3_M_fixedeq)

#cg locations
config3_cg = c2cg.get_cg(x_le_MAC,config3_class2) 
config3_cg_x=config3_cg.calc_x_cg()
config3_cg_y=config3_cg.calc_y_cg()
config3_cg_z=config3_cg.calc_z_cg()


#configuration 2
config2_class2     = c2m.Class2_weight(2,c1.N_pax[1],c1.MTOW[1],max(c1.MTOW),c1.M_carried_canard_MZF[1],min(c1.M_MZF), n_max[1],V_dive[1],c1.M_fuel[1], max(c1.T_req), c1.l_f[1],c1.d_f_inner,c1.d_f_outer,c1.l_cabin[1], l_h, c1.S, S_c2, c1.b, b_c2, S_v,S_h,c1.Cr_t,Cr_t_c2,c1.lambda_2_rad,lambda_h_2_rad, lambda_v_2_rad,lambda_c_2_rad2, c1.S_fus[1])
config2_M_structural            =config2_class2.structural_mass()
config2_M_powerplant            =config2_class2.powerplant_mass()
config2_M_powerplant            =config3_M_powerplant
config2_M_fixedeq               =config3_M_fixedeq                                  #needed as the fixed equipment needs to be the same for both
   
config2_M_winggroup             =config2_class2.get_wing_group_mass()
config2_M_fuselagegroup         =config2_class2.get_fuselage_group_mass()

                                          
config2_class2_OEW              =config2_class2.OEW(config2_M_structural,config2_M_powerplant,config2_M_fixedeq)
#cg locations
config2_cg = c2cg.get_cg(x_le_MAC,config2_class2)   
config2_cg_x=config2_cg.calc_x_cg()
config2_cg_y=config2_cg.calc_y_cg()
config2_cg_z=config2_cg.calc_z_cg()



output_file = open('output_class2iteration.dat' ,  'w')


output_file.write('config 1 M struct=' + str(config1_M_structural)+ '\n')
output_file.write('config 2 M struct=' + str(config2_M_structural)+ '\n')
output_file.write('config 3 M struct=' + str(config3_M_structural)+ '\n')
output_file.write('config 1 M pwerplant=' + str(config1_M_powerplant)+ '\n')
output_file.write('config 2 M pwerplant=' + str(config2_M_powerplant)+ '\n')
output_file.write('config 3 M pwerplant=' + str(config3_M_powerplant)+ '\n')
output_file.write('config 1 M fixedeq=' + str(config1_M_fixedeq)+ '\n')
output_file.write('config 2 M fixedeq=' + str(config2_M_fixedeq)+ '\n')
output_file.write('config 3 M fixedeq=' + str(config3_M_fixedeq)+ '\n')
output_file.write('OEW config 1 CLASS 2=' + str(config1_class2_OEW)+ '\n')
output_file.write('OEW config 2 CLASS 2=' + str(config2_class2_OEW)+ '\n')
output_file.write('OEW config 3 CLASS 2=' + str(config3_class2_OEW)+ '\n')



output_file.write('config 1 X CG =' + str(config1_cg_x)+ '\n')
output_file.write('config 2 X CG =' + str(config2_cg_x)+ '\n')
output_file.write('config 3 X CG =' + str(config3_cg_x)+ '\n')
output_file.write('config 1 y CG =' + str(config1_cg_y)+ '\n')
output_file.write('config 2 y CG =' + str(config2_cg_y)+ '\n')
output_file.write('config 3 y CG =' + str(config3_cg_y)+ '\n')
output_file.write('config 1 z CG =' + str(config1_cg_z)+ '\n')
output_file.write('config 2 z CG =' + str(config2_cg_z)+ '\n')
output_file.write('config 3 z CG =' + str(config3_cg_z)+ '\n')

output_file.close()