# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 22:27:28 2019

@author: Lisa
"""
import numpy as np
import matplotlib.pyplot as plt


fuselage_length=32.5
cabin_lenght=21.19
seatpitch=32*0.0254 # in meters

n_pax=109
n_rows=22

CG_fuselage=19.42   # from nose 
W_fuselage=33.03


CG_wing_F100=1.5155      #from LEMAC
CG_wing_F120=1.5122
Wwing_F100=17.48
Wwing_F120=18.957


MAC=3.682   




cg_fuel=0.4
Xfirst=6.4
Xlast=Xfirst+seatpitch*(n_rows-1)
    
volume_cargo_front=9.5
volume_cargo_aft=7.2
volume_cargo=volume_cargo_front+volume_cargo_aft


OEW_F100=24747
OEW_F120=24747+0.1*6766.137

MTOW=45810
W_fuel=9070

W_payload_total_F100=11993 # bagagge and cargo
W_payload_total_F120=MTOW-OEW_F120-W_fuel

W_passenger=88
W_pass_total=W_passenger*n_pax

W_payload=W_payload_total-W_pass_total
W_payload_front=W_payload*volume_cargo_front/volume_cargo
W_payload_aft=W_payload*volume_cargo_aft/volume_cargo

Xlemac_F100=18.112
Xlemac_F120=18.1163751



def loading_diagrams(Xlemac,W_wing,OEW,W_payload_total,CG_wing):
    W_payload=W_payload_total-W_pass_total
    W_payload_front=W_payload*volume_cargo_front/volume_cargo
    W_payload_aft=W_payload*volume_cargo_aft/volume_cargo
    
    
    xcg_OEW=((CG_fuselage*W_fuselage+(CG_wing+Xlemac)*W_wing)/(W_wing+W_fuselage)-Xlemac)/MAC
    
    xcg1=[xcg_OEW]
    xcg2=[xcg_OEW]
    Weight=[OEW]
    cg_cargo_front=(10.4-Xlemac)/MAC
    cg_cargo_aft=(20.8-Xlemac)/MAC
    cgfirst=(Xfirst-Xlemac)/MAC
    cglast=(Xlast-Xlemac)/MAC
    
    for i in range(2):
        cg_prev1=xcg1[i]
        cg_prev2=xcg2[i]
        W_old=Weight[i]
        if i==0:
            W_old=OEW
            W_new1=W_old+W_payload_front
            W_new2=W_old+W_payload_aft
            cg1=(W_payload_front*cg_cargo_front+W_old*cg_prev1)/(W_new1)
            cg2=(W_payload_aft*cg_cargo_aft+W_old*cg_prev2)/(W_new2)
            
        if i==1:
            W_old1=OEW+W_payload_front
            W_old2=OEW+W_payload_aft
            W_new2=W_old2+W_payload_front
            W_new1=W_old1+W_payload_aft
            cg1=(W_payload_aft*cg_cargo_aft+W_old1*cg_prev1)/(W_new1)
            cg2=(W_payload_front*cg_cargo_front+W_old2*cg_prev2)/(W_new2)
        print(W_new1,W_new2)
        xcg1.append(cg1)
        xcg2.append(cg2)
        Weight.append(W_new1)



    passenger_cg1= np.arange(cgfirst, cglast+seatpitch/MAC,seatpitch/MAC)
    passenger_cg2= np.arange(cglast,cgfirst-seatpitch/MAC,-seatpitch/MAC)
    
    print(passenger_cg1,passenger_cg2)
    
    for i in range(len(passenger_cg1)):
        cg_prev1=xcg1[i+2]
        cg_prev2=xcg2[i+2]
        W_old=Weight[i+2]
        W_new=W_old+W_passenger*2
        cg2=(W_passenger*2*passenger_cg2[i]+W_old*cg_prev2)/(W_new)
        cg1=(W_passenger*2*passenger_cg1[i]+W_old*cg_prev1)/(W_new)
        xcg1.append(cg1)
        xcg2.append(cg2)
        Weight.append(W_new)
        
        
     
    for n in range(len(passenger_cg1), (len(passenger_cg1))*2):
        cg_prev1=xcg1[n+2]
        cg_prev2=xcg2[n+2]
        W_old=Weight[n+2]
        W_new=W_old+W_passenger*2
        cg2=(W_passenger*2*passenger_cg2[n-22]+W_old*cg_prev2)/(W_new)
        cg1=(W_passenger*2*passenger_cg1[n-22]+W_old*cg_prev1)/(W_new)
        xcg1.append(cg1)
        xcg2.append(cg2)
        Weight.append(W_new)   
        
    for p in range(len(passenger_cg1)*2, (len(passenger_cg1))*3-1):
        cg_prev1=xcg1[p+2]
        cg_prev2=xcg2[p+2]
        W_old=Weight[p+2]
        W_new=W_old+W_passenger
        cg2=(W_passenger*passenger_cg2[p+1-44]+W_old*cg_prev2)/(W_new)
        cg1=(W_passenger*passenger_cg1[p-44]+W_old*cg_prev1)/(W_new)
        xcg1.append(cg1)
        xcg2.append(cg2)
        Weight.append(W_new)   
    
    
    Wtot=W_fuel+Weight[-1]   
    
    cgfuel=(Weight[-1]*xcg1[-1]+W_fuel*cg_fuel)/Wtot 
    xcg1.append(cgfuel)
    xcg2.append(cgfuel)
    Weight.append(Wtot)
    
    xcg1=np.array(xcg1)*100
    xcg2=np.array(xcg2)*100
    
    plt.figure()   
    plt.plot(xcg1, Weight, color='blue', marker='o')
    plt.plot(xcg2, Weight, color='green', marker='o')
    #plt.hlines(Weight[23],min(xcg), max(xcg),'r')
    #plt.hlines(Weight[45],min(xcg), max(xcg), 'r')
    plt.hlines(Weight[-1],min(xcg1), max(xcg2), 'r')
    plt.hlines(Weight[-2],min(xcg1), max(xcg2), 'r')
    plt.vlines(min(xcg1)-2,Weight[0], Weight[-1], 'k')
    plt.vlines(max(xcg2)+2,Weight[0], Weight[-1], 'k')
    plt.title('Load diagram for Xlemac at %2.3f m ' %(Xlemac), fontsize=14)
    plt.xlabel('Center of gravity location[% of MAC]', fontsize=12)
    plt.ylabel('Mass[kg]', fontsize=12)
    plt.show()
    
    return xcg1, xcg2,Weight


xcg1_F100,xcg2_F100,WeightF100= loading_diagrams(Xlemac_F100,Wwing_F100,OEW_F100,W_payload_total_F100,CG_wing_F100)
xcg1_F120,xcg2_F120,WeightF120= loading_diagrams(Xlemac_F120,Wwing_F120,OEW_F120,W_payload_total_F120,CG_wing_F120)

#Xlemac_plus=(Xlemac+0.1*fuselage_length)
#Xlemac_min=Xlemac-0.1*fuselage_length
#
#xcg1, xcg2,Weight1= loading_diagrams(Xlemac)
#xcg1_plus, xcg2_plus,Weight2= loading_diagrams(Xlemac_plus)
#xcg1_min, xcg2_min,Weigth3= loading_diagrams(Xlemac_min)


#min_cg=[min(xcg1_min)-2,min(xcg1)-2,min(xcg1_plus)-2]    
#max_cg=[max(xcg2_min)+2,max(xcg2)+2,max(xcg2_plus)+2]        
#long_pos=[ Xlemac_min/fuselage_length,Xlemac/fuselage_length, Xlemac_plus/fuselage_length]
#wing longitudinal position vs cg position 


#plt.figure()
#plt.plot(min_cg,long_pos,label='Most forward CG location')
#plt.plot(max_cg,long_pos,label='Most aft CG location')   
#plt.xlabel('center of gravity range [% of MAC]', fontsize=12) 
#plt.ylabel('Longitudinal wing position (Xlemac/Lfuselage) [-]', fontsize=12)
#plt.title('c.g. range vs. wing longitudinal position diagram ',fontsize=14)
#plt.legend()
##plt.show()











