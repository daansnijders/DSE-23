# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:30:42 2019

@author: Lisa
"""
#
import inputs.constants as const
from math import *

def get_wing_mass(M_MZF,b,S,Cr_t,lambda_2_rad,n_ult):
    return 0.0017*M_MZF*const.kg_to_lbs*(b*const.m_to_ft/cos(lambda_2_rad))**0.75*(1+(6.3*cos(lambda_2_rad)/(b*const.m_to_ft))**0.5)*n_ult**0.55*(b*const.m_to_ft*S*const.m_to_ft**2/(M_MZF*const.kg_to_lbs*cos(lambda_2_rad)*(Cr_t*const.m_to_ft)))**0.3
    
def get_fuselage_mass(V_dive, l_h, w_fus, h_fus, S_fus):
    return 0.021*const.K_f*(V_dive*const.ms_to_kts*l_h*const.m_to_ft/(w_fus*const.m_to_ft+h_fus*const.m_to_ft))**0.5*(S_fus*const.m_to_ft**2)**1.2

def get_horizontaltail_mass(K_h,S_h,V_dive,lambda_h_2_rad):
    return K_h*S_h*(3.18*(((S_h*const.m_to_ft**2)**0.2* V_dive*const.ms_to_kts)/(1000*(cos(lambda_h_2_rad))**0.5))-0.287)

def get_verticaltail_mass(K_v,S_v,V_dive,lambda_v_2_rad):
    return K_v*S_v*(3.18*(((S_v*const.m_to_ft**2)**0.2* V_dive*const.ms_to_kts)/(1000*(cos(lambda_v_2_rad))**0.5))-0.287)

def get_nacelle_mass(T_req_TO):
    return 0.065*T_req_TO # this is total weight of ALL nacelles

def get_landinggear_mass(K_gr,Ag,Bg,Cg,Dg,M_TO):
    return K_gr*(Ag+Bg*(M_TO*const.kg_to_lbs)**0.75+Cg*M_TO*const.kg_to_lbs+Dg*(M_TO*const.kg_to_lbs)**1.5)

def get_structural_mass(M_wing,M_fuselage,M_horizontaltail,M_verticaltail,M_nacelle,M_landinggear, M_canard):
    return M_wing+M_fuselage+M_nacelle+M_horizontaltail+M_verticaltail+M_landinggear+M_canard
