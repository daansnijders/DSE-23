# -*- coding: utf-8 -*-
"""
Created on Tue May 14 13:42:00 2019

@author: Stijn
"""


def get_z_cg(d_f_outer):
    return [d_f_outer[i]/2 for i in range(3)]


def get_y_cg():
    return [0, 0, 0]


def get_x_cg(l_f, MTOW, MAC, concept_3 = False):
    x_c_wcg = 0.205
    x_c_oewcg = 0.225
    l_engine = 3.184
    x_f_cg_l_f = 0.435
    
    # Mass wing group
    M_wing = [0.092 * MTOW[i] for i in range(3)]
    M_eng = [0.083 * MTOW[i] for i in range(3)]
    M_wing_group = [M_wing[i] + M_eng[i]  for i in range(3)]
    
    
    # Fuselage
    x_cg_fuselage = [x_f_cg_l_f * l_f[i] for i in range(3)]
    M_fuselage = [0.234 * MTOW[i] for i in range(3)]
    
    # Tail
    x_cg_tail = [0.9 * l_f[i] for i in range(3)]
    M_tail = [0.024 * MTOW[i] for i in range(3)]
    
    # Fuselage group
    M_fuselage_group = [M_fuselage[i] + M_tail[i] for i in range(3)]
    x_cg_fuselage_group = [(M_tail[i]*x_cg_tail[i] \
                            +M_fuselage[i]*x_cg_fuselage[i])/ M_fuselage_group[i]  for i in range(3)]
    
    x_le_MAC_org = [x_f_cg_l_f * l_f[i] + MAC[i] * (x_c_wcg*(M_wing_group[i]/M_fuselage_group[i])\
                    -x_c_oewcg * (1+M_wing_group[i]/M_fuselage_group[i])) for i in range(3)]
    
    x_le_MAC = [x_le_MAC_org[0], x_le_MAC_org[0] \
                + (l_f[1] - l_f[0]), x_le_MAC_org[0] + (l_f[2] - l_f[0])]
    if concept_3:
        x_le_MAC = [max(x_le_MAC_org) for i in range(3)]

    x_cg_wing = [x_le_MAC[i] + x_c_wcg * MAC[i] for i in range(3)]
    x_cg_eng = [x_le_MAC[i]-0.6*l_engine for i in range(3)]
    
    x_cg_wing_group = [(M_wing[i]*x_cg_wing[i] \
                        +x_cg_eng[i]*M_eng[i])/M_wing_group[i] for i in range(3)]

    x_cg = [(M_wing_group[i]*x_cg_wing_group[i] \
             + M_fuselage_group[i]*x_cg_fuselage_group[i])/(M_wing_group[i] \
             + M_fuselage_group[i]) for i in range(3)]
    
    return x_cg