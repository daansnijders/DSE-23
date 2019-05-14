# -*- coding: utf-8 -*-
"""
Created on Thu May  9 10:58:36 2019

@author: Sybren
"""

""" Fuselage Design """
from math import *
from modules.initialsizing_weights import *
from inputs.constants import *



def get_d_f_inner(N_sa, seat_width, N_aisle, armrest, aisle_width, s_clearance):
    return [N_sa*seat_width + (N_sa+N_aisle+1)*armrest + N_aisle*aisle_width + 2*s_clearance for i in range(3)]

def get_d_f_outer(d_f_inner):
    return [1.045*d_f_inner[i] + 0.084 for i in range(3)]

def get_l_cabin(N_pax,N_sa):
    return [N_pax[i]/N_sa*1.08 for i in range(3)]

def get_l_nose(d_f_outer):
    return [1.85 * d_f_outer[i] for i in range(3)]

def get_l_tail(d_f_outer):
    return [1.6 * d_f_outer[i] for i in range(3)]

def get_l_tailcone(d_f_outer):
    return [3.5 * d_f_outer[i] for i in range(3)]

def get_l_fuselage(l_cockpit, l_cabin, l_tail):
    return [l_cockpit + l_cabin[i] + l_tail[i] for i in range(3)]

def get_overhead_volume(l_cabin):
    return [2*0.2*l_cabin[i]*0.74 for i in range(3)]



def get_cargo_storage(R_f,cabinlength,Npax,V_os):
    p = R_f- h_max - h_floor        	#[m] Distance between the lower point of the inner fuselage and the floor
    phi = 2 * acos(1-p/R_f)             #[rad] Angle between the two connection points of the fuselage and floor
    A_cc = 0.5 * R_f**2*(phi - sin(phi))           #[Cargo hold area]m^2
    V_cc = [cabinlength[i]*0.45*A_cc for i in range(3)]    #[m3] Total available cargo hold volume
    Wtot_carry_on =[ N_pax[i] * W_carry_on  for i in [0,1,2]]  #[kg] Total carry-on weight
    Wtot_check_in =[ N_pax[i] * W_check_in  for i in [0,1,2] ]  #[kg] Total check-in weight 
    V_carry_on[i] = [Wtot_carry_on[i] / rho_lugg for i in range(3) ]             #[m3] Total carry-on volume needed
    V_check_in[i]= [Wtot_check_in[i] / rho_lugg  for i in range(3)]             #[m3] Total check-in volume needed
    V_cargo[i] = [V_cc[i] - (V_carry_on[i] + V_check_in[i] - V_os[i]) for i in [0,1,2]] #[m3] Total available cargo volume
    M_cargo[i] = [V_cargo[i] * rho_cargo  for i in [0,1,2]]
    M_payload[i] = [N_pax[i]* (W_carry_on[i] + W_check_in[i] + W_pax[i]) + M_cargo[i]  for i in [0,1,2]]#[kg] Total payload weight
    return V_cc , M_cargo,M_payload
    
def get_fuselage_length(l_cockpit, l_tail, l_cabin) :       
    return l_cockpit + l_tail + l_cabin




