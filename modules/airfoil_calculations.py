# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:01:38 2019

@author: Anique
"""
import math
import inputs.constants as const
import inputs.performance_inputs as perf

#### airfoil calculates different Reynold numbers for different locations 
#### and calculates Cl design for the airfoils and CL max for the wing
def airfoil( Ct, Cr, MTOW, FF1, FF2, FF3, FF4, FF5, S, sweep_le, sweep_c2, b, Taper, A, Cl_max):
    avgC = (Cr + Ct)/2    # Average chord [m]
    #Re during cruise
    Re1 = (const.rho*const.V_cruise*Cr)/const.mu_37   # Reynolds number at root chord [-]
    Re2 = (const.rho*const.V_cruise*avgC)/const.mu_37    # Reynolds number at avg chord [-]
    Re3 = (const.rho*const.V_cruise*Ct)/const.mu_37      # Reynolds number at tip chord [-]
    
    #Re during take-off
    Reto1 = (const.rho_0*perf.Vto1*const.kts_to_ms*Cr)/const.mu_sl     # Reynolds number at root chord [-]
    Reto2 = (const.rho_0*perf.Vto1*const.kts_to_ms*avgC)/const.mu_sl    # Reynolds number at root chord [-]
    Reto3 = (const.rho_0*perf.Vto1*const.kts_to_ms*Ct)/const.mu_sl      # Reynolds number at root chord [-]
    q = 0.5*const.rho*const.V_cruise**2        # Dynamic pressure [Pa]
    
    WSbegin = [(MTOW[i]*g*FF1*FF2*FF3*FF4)/S   for i in range(3)]        # Wing Loading begin cruise  
    WSend = [(MTOW[i]*g*FF1*FF2*FF3*FF4*FF5)/S  for i in range(3)]        # Wing loading end cruise
    CLdes =[ 1.1*(1/q)*(0.5*(WSbegin[i] + WSend[i]))   for i in range(3)] # Wing CL design [-]
                                             
    Veff = const.V_cruise*math.cos(sweep_le)             # Effective velocity [m/s]
    qeff = 0.5*const.rho*Veff**2                  # Effective Dynamic Pressure [Pa]
    Cl_des = [(q*CLdes[i])/qeff for i in range(3)]                   # Airfoil Cl design [-]
    
    beta = math.sqrt(1-const.M_cruise**2)      # Prandtl-Glauert compressibility correction factor
    CL_alpha = (2*math.pi*A)/(2 + math.sqrt(4+(A*beta/0.95)**2*(1+(math.tan(sweep_c2)**2)/beta**2)))
    
    #Wing CLmax for different Re numbers for cruise
    CLmax = [0.8*Cl_max[i]-0.24 for i in range(3)]
    
    CLmaxto = [0.8*Cl_max[i] for i in range(3)]
    return(Reto1, Re1, Reto3, CLdes, Cl_des, CL_alpha, CLmax, CLmaxto)
       

#### drag1 calculates CD0 for each concept and configuration as well as CDcruise FOR CONCEPT 1   
def drag1(A, S, S_h, S_v, l_nose, l_tailcone, l_fuselage, D, Dnacel, Lnacel, sweep_le, CLdes):
    Wing = []
    
    for i in range(3):
        if i == 0:
            Wing.append(1.07*2*S*0.003)
        if i == 1:
            Wing.append(1.07*2*S*0.003)
        if i == 2:
            Wing.append(1.07*2*S*0.003 + (((45*const.inch_to_m)/2)**2*math.pi*2 + (45*const.inch_to_m*math.pi*308.4*const.inch_to_m))*0.006)
            
    l2 = [l_fuselage[i] - l_nose - l_tailcone for i in range(3)]
    Fuselage = [0.0024*(math.pi*D/4)*(1/(3*l_nose**2)*((4*l_nose**2+D**2/4)**(1.5)-D**3/8)-D+4*l2[i]+2*math.sqrt(l_tailcone**2+D**2/4)) for i in range(3)]
    Nacelle = (2*(Dnacel/2)**2*math.pi + Dnacel*math.pi*Lnacel)*0.0060
    Tailplane = [1.05*2*(S_h[i] + S_v[i])*0.0025  for i in range(3)]
    CD0 = [1/S*(Wing[i] + Fuselage[i] + Nacelle + Tailplane[i])*1.1  for i in range(3)]
    
    e = 4.61*(1-0.045*A**0.68)*(math.cos(sweep_le))**0.15 -3.1
    
    CDcruise = [CD0[i] + 1/(math.pi*A*e)*CLdes[i]**2 for i in range(3)]
    
    LoverD = [CLdes[i]/CDcruise[i] for i in range(3)]
    return(CD0, CDcruise, LoverD, Wing, Fuselage, Nacelle, Tailplane)
    

#### drag calculates CD0 for each concept and configuration as well as CDcruise FOR CONCEPT 2  
def drag2(A, S, S_h, S_v, l_nose, l_tailcone, l_fuselage, D, Dnacel, Lnacel, sweep_le, CLdes):
    Wing = []
    for i in range(3):
        if i == 0:
            Wing.append(1.07*2*147.5107162*0.003)
        if i == 1:
            Wing.append(1.07*2*S[i]*0.003)
        if i == 2:
            Wing.append(1.07*2*S[i]*0.003 + (((45*const.inch_to_m)/2)**2*math.pi*2 + (45*const.inch_to_m*math.pi*308.4*const.inch_to_m))*0.006)
        
    l2 = [l_fuselage[i] - l_nose[i] - l_tailcone[i] for i in range(3)]
    Fuselage = [0.0024*(math.pi*D[i]/4)*(1/(3*l_nose[i]**2)*((4*l_nose[i]**2+D[i]**2/4)**(1.5)-D[i]**3/8)-D[i]+4*l2[i]+2*math.sqrt(l_tailcone[i]**2+D[i]**2/4)) for i in range(3)]
    Nacelle = (2*(Dnacel/2)**2*math.pi + Dnacel*math.pi*Lnacel)*0.0060
    Tailplane = [1.05*2*(S_h[i] + S_v[i])*0.0025  for i in range(3)]
    CD0 = [1/147.5107162*(Wing[i] + Fuselage[i] + Nacelle + Tailplane[i])*1.1  for i in range(3)]
    
    e = [4.61*(1-0.045*A[i]**0.68)*(math.cos(sweep_le[i]))**0.15 -3.1 for i in range(3)]
    
    CDcruise = [CD0[i] + 1/(math.pi*A[i]*e[i])*CLdes[i]**2 for i in range(3)]
    
    LoverD = [CLdes[i]/CDcruise[i] for i in range(3)]
    return(CD0, CDcruise, LoverD)
    
#### drag calculates CD0 for each concept and configuration as well as CDcruise FOR CONCEPT 3  
def drag3(A, S, S_h, S_v, l_nose, l_tailcone, l_fuselage, D, Dnacel, Lnacel, sweep_le, CLdes):
    Wing = [1.07*2*S[i]*0.003 for i in range(3)]  
    l2 = [l_fuselage[i] - l_nose[i] - l_tailcone[i] for i in range(3)]
    Fuselage = [0.0024*(math.pi*D[i]/4)*(1/(3*l_nose[i]**2)*((4*l_nose[i]**2+D[i]**2/4)**(1.5)-D[i]**3/8)-D[i]+4*l2[i]+2*math.sqrt(l_tailcone[i]**2+D[i]**2/4)) for i in range(3)]
    Nacelle = []
    for i in range(3):
        if i == 2:
            Nacelle.append((2*(Dnacel/2)**2*math.pi + Dnacel*math.pi*Lnacel)*0.006 + (((45*const.inch_to_m)/2)**2*math.pi*2 + (45*const.inch_to_m*math.pi*308.4*const.inch_to_m))*0.006)
        else:
            Nacelle.append((2*(Dnacel/2)**2*math.pi + Dnacel*math.pi*Lnacel)*0.0060)
            
    Tailplane = [1.05*2*(S_h[i] + S_v[i])*0.0025  for i in range(3)]
    CD0 = [1/S[0]*(Wing[i] + Fuselage[i] + Nacelle[i] + Tailplane[i])*1.1  for i in range(3)]
    
    e = [4.61*(1-0.045*A[i]**0.68)*(math.cos(sweep_le[i]))**0.15 -3.1 for i in range(3)]
    
    CDcruise = [CD0[i] + 1/(math.pi*A[i]*e[i])*CLdes[i]**2 for i in range(3)]
    
    LoverD = [CLdes[i]/CDcruise[i] for i in range(3)]
    return(CD0, CDcruise, LoverD)