# -*- coding: utf-8 -*-
"""
Created on Mon May 13 11:44:09 2019

@author: thong

Analysis of Concept 1 fuselage joint.
120pax, 4000km

Assumptions:
1) Fuselage is cylindrical
2) Fuselage structural mass and payload is distributed equally
3) Point load, unclamped, fuel and wing weight conincide with wing lift, tail lift
    conincide with tail weight, canard lift conincide with canard weight
4) Stringers carry longitudinal stress
5) z_CG at the center of fuselage, symmetrical in x,y,zplane
6) Aluminium 7075-T6, yield strength 503MPa
7) Torsion neglected
8) Effect of buckling neglected
9) Calculation all done to meet limit load according the tensile yield strength
10) Thickness of webframes on each side of joint is a quarter of the length of
    reinforcement and the height of the web frames is 3 times the height of stringer
    which is square root of the area.
11) Length of reinforcement takes into account of the fittings and multiple frames
    and increased size of stringers near it for reinforcement
    
Description:
1) The program begins with plugging in the fixed parameters
2) The location and mass of components are inputted into the program
3) Difference in OEW and sum of components is included into the mass of fuselage
4) Force and moment equilibrium is performed to obtain lift for canard and tail
5) Internal Shear Diagram produced
6) Internal Bending moment diagram produced
7) 8) FUNCTION 1 : Moment of inertia created
9) Longtitudinal stress due to pressure difference (37000ft and 2438ft CS25)
10) Bending stress calculated
11) Number of stringers needed at each point of the fuselage is calculated which 
    depends on surface area of the stringers (Idealized boom structure)
12) Skin Thickness calculated due to radial stress by pressure
13) Shear Stress calculated.
    
Input:
1) MTOW @ 90,4000km                 [kg]
2) MTOW @ 120,4000km                [kg]
3) OEW @ 120,4000km                 [kg]
4) Fuselage length @ 120,4000km     [m]
5) Fuselage length @ 90,4000km      [m]
6) Fuselage Diameter @ 120,4000km   [m]
7) Fuselage Mass @ 120,4000km       [kg]
8) Payload Mass @ 120,4000km        [kg]
9) Wing group mass @ 120,4000km     [kg]
10) Wing location @ 120,4000km      [m]
11) Fuel Mass @ 120,4000km          [kg]
12) Tail mass @ 120,4000km          [kg]
13) Tail location @ 120,4000km      [m]
14) Canard location @ 120,4000km    [m]
15) Wing mass @ 120,4000km          [kg]
16) Wing mass @ 90,4000km           [kg]

Output:
List[ Mass increase per joint, number of bolts at front joint, number of bolts at aft joints]

To be done:
1)Verification
2)Influence of fatigue
3)Different loadings
4)Total mass of reinforcement needed (Size of web frames and bolts and fittings)
5)Influence of other failure mode (bearing)
-------------------------------------------------------------------------------
"""

import numpy as np
import matplotlib.pyplot as plt
from isa import isa
#from mpl_toolkits.mplot3d import Axes3D

#"""
#Parameters
#Previous update time: - 
#Current update time: 13/5/2019, 12:11PM
#-------------------------------------------------------------------------------
#"""

def fuselage1_mass_analysis(MTOW1,MTOW2,OEW2,l_fuselage2,l_fuselage1,d_fuselage2,\
                           m_fuselage2,m_payload2,m_winggroup2,x_wing2,m_fuel2,\
                           m_tail2,x_tail2,x_canard2,m_wing2,m_wing1):
    print("********** FUSELAGE MASS ANALYSIS CONCEPT 1 **********")
    # MTOW [kg] @ 90,4000km
    MTOW1 = MTOW1 #49391.28866
    # MTOW [kg] @ 120,4000km
    MTOW2 = MTOW2 #68264.27835
    # OEW [kg] @ 90,4000km
    OEW = OEW2 #38729.81
    # Safety factor
    safety_factor = 2
    # Load Factor
    load_factor = 2.1
    # Gravitation acceleration [m/s^2]
    g = 9.81
    
    # Fuselage length [m]
    l_fuselage = l_fuselage2 #29.335
    l_fuselage_short = l_fuselage1 #19.44        
    # Fuselage diameter [m]
    d_fuselage = d_fuselage2 #3.685
    # Fuselage structure mass [kg]
    m_fuselage = m_fuselage2 #15973.84113
    
    # Payload mass [kg]
    m_payload = m_payload2 #12925.76943
    
    # Wing group mass [kg]
    m_winggroup = m_winggroup2 #11946.24871
    # Wing group location [m]
    x_wing = x_wing2#16.1349902               #16.1349902
    # Fuel mass [kg]
    m_fuel = m_fuel2 #16608.69892
    
    
    # Tail mass [kg]
    m_tail = m_tail2 #1638.34268
    # Tail location [m]
    x_tail = x_tail2 #26.4015            #26.4015
    
    # Wing lift [N]
    L_wing = MTOW1*g
    
    # Canard mass [kg]
    m_canard = m_wing2-m_wing1#(157.4050988-113.8873926)/157.4050988*6280.313608
    # Canard location [m]
    x_canard = x_canard2
    
    """Bolts Details"""
    l_reinforcement = 100           # [mm]
    rho_material = 2810             # [kg/m^3]
    # Tensile Yield Strength [MPa] includes fatigue
    yield_strength = 159
    
    """Pressure at 2438ft [Pa]"""
    P1 = isa(2438/3.281)[1]
    """Pressure at 37000ft [Pa]"""
    P2 = isa(37000/3.281)[1]
    
    """ Fuselage cross-section design"""
    # Fuselage Stringer Area [mm^2]
    stringer_area = 400
    
    """
    Analysis (Step Size)
    -------------------------------------------------------------------------------
    """
    
    n = 1000
    step_size = l_fuselage/n
    
    """
    -------------------------------------------------------------------------------
    Adjusting m_fuselage to include main and nose landing gear
    """
    delta_oew = OEW - (m_fuselage+m_winggroup+m_canard+m_tail)
    print("Delta OEW = ", delta_oew, " kg")
    m_fuselage+=delta_oew
    
    
    """ Force and Moment Equilibrium """
    delta_F = (m_fuselage+m_payload+m_winggroup+m_canard+m_tail+m_fuel)*g-L_wing
    delta_M = -m_fuselage*l_fuselage/2*g-m_payload*l_fuselage/2*g-\
                m_winggroup*x_wing*g-m_canard*g*x_canard-m_tail*g*x_tail-m_fuel*g*x_wing+L_wing*x_wing
    
    # Canard lift [N]
    L_canard = (-delta_M-x_tail*delta_F)/(x_canard-x_tail)
    # Tail lift [N]
    L_tail = delta_F-L_canard
    
    EQ_F = L_wing+L_canard+L_tail-(m_fuselage+m_payload+m_winggroup+m_canard+m_tail+m_fuel)*g
    EQ_M = L_wing*x_wing+L_canard*x_canard+L_tail*x_tail\
            -(m_fuselage*l_fuselage/2+m_payload*l_fuselage/2+\
              m_winggroup*x_wing+m_canard*x_canard+m_tail*x_tail+m_fuel*x_wing)*g
    print("Summation Force = ", EQ_F)
    print("Summation Moment = ", EQ_M)
    
    """
    Point Loads: Canard weight, Canard Lift,Fuel Weight, Wing Group weight, Wing Lift, Tail weight, Tail Lift
    Distributed Loads: Fuselage, Payload
    """
    P_loads = np.array([-m_canard*g, L_canard, -m_fuel*g, -m_winggroup*g, L_wing, -m_tail*g, L_tail])*load_factor
    P_x = np.array([x_canard,x_canard,x_wing,x_wing,x_wing,x_tail,x_tail])
    D_loads = np.array([-m_fuselage*g/l_fuselage, -m_payload*g/l_fuselage])*load_factor
    D_x = np.linspace(0,l_fuselage,n)
    
    """
    Internal Shear Diagram
    """
    
    x = 0
    x_lst = []
    shear_lst = []
    shear = 0
    for i in range(n):
        x += step_size
        for j in range(len(P_loads)):
            if x-step_size/2<P_x[j]<=x+step_size/2:
                shear += P_loads[j]
        shear += sum(D_loads)*step_size
        x_lst.append(x)
        shear_lst.append(shear)
    
    
    plt.figure("Internal Shear")
    plt.plot(x_lst,shear_lst)
    plt.title("Internal Shear Diagram")
    plt.xlabel("x [m]")
    plt.ylabel("Internal Shear [Nm]")
    plt.show()
    
    """
    Internal Moment Diagram
    Macaulay's Step Function
    """
    moment_lst = []
    moment = 0
    for i in range(n):
        x_loc = i*step_size
        moment = sum(D_loads)*0.5*x_loc**2
        number = 0
        for j in range(len(P_x)):
            if x_loc>P_x[j]:
                moment +=  (x_loc-P_x[j])*P_loads[j]
                number += 1
    
        moment_lst.append(moment)
    print("End internal moment = ", moment_lst[-1])
    
    
    plt.figure("Internal Moment Diagram")
    plt.plot(x_lst,(moment_lst))
    plt.title("Internal Moment Diagram")
    plt.xlabel("x [m]")
    plt.ylabel("Internal Moment [Nm]")
    plt.show()
    
    
    
    """Calculate MMOI of fuselage"""
    def mmoi(stringer_number,radius,stringer_area):   #[-,mm,mm^2]
        angle = 2*np.pi/stringer_number
        stringer_lst = []
        for i in range(stringer_number):
            y = np.sin(angle*i)*radius
            z = np.cos(angle*i)*radius
            stringer_lst.append([y,z,stringer_area])
        stringer_lst = np.array(stringer_lst)
        mean_y = sum(stringer_lst[:,2]*stringer_lst[:,1])/sum(stringer_lst[:,2])
        MMOI = sum((stringer_lst[:,1])**2*stringer_lst[:,2])
        return [MMOI,mean_y]
    
    """Calculate number of stringers needed"""
    fuselage_radius = d_fuselage/2     #[m]
    stringer_number = 2
    stringer_lst = []
    stress_lst = []
    
    i = 0
    while i < len(moment_lst):
        k = moment_lst[i]
        long_stress = (P1-P2)*np.pi*(fuselage_radius)**2/(stringer_number*stringer_area)
        bending_stress = (abs(k)*1000*fuselage_radius*1000/(mmoi(stringer_number,fuselage_radius*1000,stringer_area)[0]))
        normal_stress = ((bending_stress)+long_stress)*safety_factor
        
        if yield_strength<normal_stress:
            stringer_number += 1
        else:
            stringer_lst.append(stringer_number)
            stress_lst.append(normal_stress)
            i += 1
            stringer_number = 1
            
    
    plt.figure("Number of Stringers")
    plt.plot(x_lst,stringer_lst)
    plt.title("Number of Stringers")
    plt.xlabel("x location(m)")
    plt.ylabel("Number of Stringers")
    plt.show()
    
    """Number of bolts at the joint"""
    
    x_joint1 = -(l_fuselage-l_fuselage_short)/2+x_canard         #[m]
    x_joint2 = (l_fuselage-l_fuselage_short)/2+x_canard          #[m]
    joint1 = min(range(len(x_lst)), key=lambda i: abs(x_lst[i]-x_joint1))
    joint2 = min(range(len(x_lst)), key=lambda i: abs(x_lst[i]-x_joint2))
    n_joint1 = stringer_lst[joint1]
    n_joint2 = stringer_lst[joint2]
    
    
    """ Thickness of skin """
    t_skin = (P1-P2)*safety_factor*fuselage_radius/(yield_strength*10**6)
    print("Skin thickness= ",t_skin*10**3, " mm")
    
    """ Maximum Normal Stress """
    print("Maximum Normal Stress = ", max(stress_lst), " MPa")
        
    """ Shear Stress """
    def shear(shear_lst, stringer_lst,fuselage_radius, stringer_area):
        radius = fuselage_radius*1000
        shear_stress = []
        for i in range(len(shear_lst)):
            stringer = []
            for j in range(stringer_lst[i]):
                angle = 2*np.pi/stringer_lst[i]
                y = np.cos(angle*j)*radius
                z = np.sin(angle*j)*radius
                Q = -radius*radius*np.cos(angle*j) #<------ thickness removed
                q = - shear_lst[i] * safety_factor*Q/mmoi(stringer_lst[i],radius,stringer_area)[0]
                stringer = [shear_lst[i],x_lst[i],y,z,Q,q]
                shear_stress.append(stringer)
        shear_stress = np.array(shear_stress)
        return shear_stress
    
    shear_stress_all = shear(shear_lst, stringer_lst,fuselage_radius, stringer_area)
    #fig = plt.figure()
    #ax = fig.add_subplot(111,projection='3d')
    #ax.scatter(shear_flow_all[:,1],shear_flow_all[:,2],shear_flow_all[:,3],c=abs(shear_flow_all[:,-1]))
    #ax.set_xlabel('x')
    #ax.set_ylabel('y')
    #ax.set_zlabel('z')
        
    joint1_normal_stress = stress_lst[joint1]
    joint2_normal_stress = stress_lst[joint2]
    
    joint1_shear_idx = np.where(shear_stress_all[:,1] == shear_stress_all[np.abs(shear_stress_all[:,1] - x_joint1).argmin(),1])
    joint1_shear = (shear_stress_all[joint1_shear_idx,:])
    joint1_maxshear = np.max(joint1_shear[:,:,-1])
    
    joint2_shear_idx = np.where(shear_stress_all[:,1] == shear_stress_all[np.abs(shear_stress_all[:,1] - x_joint2).argmin(),1])
    joint2_shear = (shear_stress_all[joint2_shear_idx,:])
    joint2_maxshear = np.max(joint2_shear[:,:,-1])
    print("")
    print("Joint 1 Forward")
    print("---------------")
    print("Number of bolts = ", n_joint1)
    print("Normal Stress = ", joint1_normal_stress, " MPa")
    print("Shear Stress = ", joint1_maxshear, " MPa")
    
    print("Joint 2 Aft")
    print("---------------")
    print("Number of bolts = ", n_joint2)
    print("Normal Stress = ", joint2_normal_stress, " MPa")
    print("Shear Stress = ", joint2_maxshear, " MPa")
    
    
    bolt_mass = (2*n_joint2)*stringer_area/10**6*l_reinforcement/1000*rho_material
    
    print("Bolt mass =", bolt_mass, " kg")
    
    web_mass = (d_fuselage*np.pi)*(stringer_area**0.5/1000)*3*l_reinforcement/1000*rho_material
    print("Web mass = ", web_mass, " kg")
    
    skin_mass = t_skin*d_fuselage*1000*np.pi*l_reinforcement*10**-6*rho_material
    print("Skin mass = ", skin_mass, " kg")
    
    total_mass = web_mass+bolt_mass+skin_mass
    print("Total mass per joint = ", web_mass+bolt_mass+skin_mass, " kg")
    
    return [total_mass,n_joint1,n_joint2]
    
    
