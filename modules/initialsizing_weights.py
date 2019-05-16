# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:25:55 2019

@author: Stijn
"""

def get_MTOW(OEW):
    MTOW = [(OEW[0]+1000)/0.582, (OEW[1]+1000)/0.582, (OEW[2]+1000)/0.582]      # [kg] #delete this one and put every MTOW as an input ?
    return MTOW

def get_TOW(OEW,M_payload, M_ff):
    return [(M_payload[i]+OEW[i])/Mff[i] for i in range(3)]
    
def get_M_fuel(MTOW,M_ff):
    M_fuel = [MTOW[0]*(1-M_ff[0]), MTOW[1]*(1-M_ff[1]), MTOW[2]*(1-M_ff[2])]    # [kg]    
    return M_fuel

def get_T_req(T_W, MTOW):
    return [T_W[i] * MTOW[i] * 9.80665 for i in [0,1,2]]                                  # [N]

def get_M_payload_available(MTOW, OEW, M_fuel):
    return [MTOW[i] - OEW[i] - M_fuel[i] for i in [0,1,2]]                  # define payload between capacity or only passengers
