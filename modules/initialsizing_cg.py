# -*- coding: utf-8 -*-
"""
Created on Tue May 14 13:42:00 2019

@author: Stijn
"""
import inputs.constants as const

def get_z_cg(d_f_outer):
    return d_f_outer/2


def get_y_cg():
    return [0, 0, 0]

def get_mass_winggroup(MTOW):
    M_wing = [0.092 * MTOW[i] for i in range(3)]
    M_eng = [0.083 * MTOW[i] for i in range(3)]
    M_wing_group = [M_wing[i] + M_eng[i]  for i in range(3)]
    return M_wing, M_eng, M_wing_group

def get_mass_fuselage(MTOW,l_f):
    M_fuselage = [0.234 * MTOW[i] for i in range(3)]
    x_cg_fuselage = [const.x_f_cg_l_f * l_f[i] for i in range(3)]
    return M_fuselage, x_cg_fuselage


def get_mass_tail(MTOW,l_f):
    M_tail = [0.024 * MTOW[i] for i in range(3)]
    x_cg_tail = [0.9 * l_f[i] for i in range(3)]
    return M_tail,x_cg_tail

def get_mass_fuselagegroup(M_fuselage,M_tail,x_cg_fuselage,x_cg_tail):
    M_fuselage_group = [M_fuselage[i] + M_tail[i] for i in range(3)]
    x_cg_fuselage_group = [(M_tail[i]*x_cg_tail[i] \
                            +M_fuselage[i]*x_cg_fuselage[i])/ M_fuselage_group[i]  for i in range(3)]
    return M_fuselage_group, x_cg_fuselage_group


def get_x_le_MAC(l_f,MAC,M_wing_group, M_fuselage_group):
    x_le_MAC_org = [const.x_f_cg_l_f * l_f[i] + MAC * (const.x_c_wcg*(M_wing_group[i]/M_fuselage_group[i])\
                    -const.x_c_oewcg * (1+M_wing_group[i]/M_fuselage_group[i])) for i in range(3)]
    
    x_le_MAC = [x_le_MAC_org[0], x_le_MAC_org[0] \
                + (l_f[1] - l_f[0]), x_le_MAC_org[0] + (l_f[2] - l_f[0])]

    return x_le_MAC

def get_cg_winggroup(x_le_MAC, MAC,M_wing, M_eng, M_wing_group ):
    x_cg_wing = [x_le_MAC[i] +const. x_c_wcg * MAC for i in range(3)]
    x_cg_eng = [x_le_MAC[i]-0.6*const.l_eng for i in range(3)]
    
    x_cg_wing_group = [(M_wing[i]*x_cg_wing[i] \
                        +x_cg_eng[i]*M_eng[i])/M_wing_group[i] for i in range(3)]
        
    return x_cg_wing,x_cg_eng,x_cg_wing_group
    
def get_x_cg(M_wing_group, M_fuselage_group,x_cg_wing_group, x_cg_fuselage_group):
    x_cg = [(M_wing_group[i]*x_cg_wing_group[i] \
             + M_fuselage_group[i]*x_cg_fuselage_group[i])/(M_wing_group[i] \
             + M_fuselage_group[i]) for i in range(3)]
    
    return x_cg