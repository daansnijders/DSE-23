# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:17:24 2019

@author: Lisa
"""



'READ AND IMPORT ALL THE OUTPUTS '

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

MAC_c2=values[31]
MAC_c3 =values[32]
l_c2=values[33]
l_c3=values[34]

x_mlg1 =values[35]
x_mlg2 =values[36]
x_mlg3 =values[37]
x_mlg=[x_mlg1,x_mlg2,x_mlg3]

l_m1 =values[38]
l_m2 =values[39]
l_m3 =values[40]
l_n1 =values[41]
l_n2 =values[42]
l_n3 =values[43]

L_strut_mlg =values[44]
L_strut_nlg =values[45]
D_strut_mlg =values[46]
D_strut_nlg =values[47]
CL_clean_max=values[48]
CL_flaps_max=values[49]
CL_TO=values[50]
CL_land=values[51]
CL_cruise1=values[52]
CL_cruise2=values[53]
CL_cruise3=values[54]
Cd_cruise1=values[55]
Cd_cruise2=values[56]
Cd_cruise3=values[57]


