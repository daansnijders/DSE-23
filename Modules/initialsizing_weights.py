# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:25:55 2019

@author: Stijn
"""

def get_MTOW(OEW):
    MTOW = [(OEW[0]+1000)/0.582, (OEW[1]+1000)/0.582, (OEW[2]+1000)/0.582]    # [kg]
    return MTOW
    
def get_M_fuel(MTOW,M_ff):
    M_fuel = [MTOW[0]*(1-M_ff[0]), MTOW[1]*(1-M_ff[1]), MTOW[2]*(1-M_ff[2])]    # [kg]    
    return M_fuel

def get_T_req(T_W, MTOW):
    return [T_W[i] * MTOW[i] for i in [0,1,2]]

def get_M_payload(MTOW, OEW, M_fuel):
    return [MTOW[i] - OEW[i] - M_fuel[i] for i in [0,1,2]]