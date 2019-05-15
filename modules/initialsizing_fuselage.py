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



def get_cargo_volume(R_f,cabinlength):
    p = [R_f[i]- h_max - h_floor for i in range(3)]        	#[m] Distance between the lower point of the inner fuselage and the floor
    phi = [2 * acos(1-p[i]/R_f[i]) for i in range(3)]          #[rad] Angle between the two connection points of the fuselage and floor
    A_cc = [0.5 * R_f[i]**2*(phi [i]- sin(phi[i]))     for i in range(3)]        #[Cargo hold area]m^2
    V_cc = [cabinlength[i]*0.45*A_cc[i] for i in range(3)]  #[m3] Total available cargo hold volume
    return V_cc

def get_masses_volumes(N_pax, V_cc, V_os):
    Wtot_carry_on = [ N_pax[i] * W_carry_on  for i in range(3)]  #[kg] Total carry-on weight
    Wtot_check_in = [ N_pax[i] * W_check_in  for i in range(3) ]  #[kg] Total check-in weight 
    V_carry_on = [Wtot_carry_on[i] / rho_lugg for i in range(3) ]             #[m3] Total carry-on volume needed
    V_check_in= [Wtot_check_in[i] / rho_lugg  for i in range(3)]             #[m3] Total check-in volume needed
    return Wtot_carry_on, Wtot_check_in, V_carry_on, V_check_in


def get_available_cargo_volume(V_cc,V_os,V_carry_on, V_check_in):
    return[V_cc[i] - (V_carry_on[i] + V_check_in[i] - V_os[i]) for i in range(3)] #[m3] Total available cargo volume

def get_cargo_mass(N_pax,M_payload):
    M_cargo =[ M_payload[i]-N_pax[i]*(W_pax+W_carry_on+W_check_in)  for i in range(3)]               #[kg] Total cargo weight
    return M_cargo


    





