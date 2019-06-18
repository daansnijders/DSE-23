# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:29:12 2019

@author: Lisa
"""

'plot the loading diagram FOR THE ITERATIONS'
import modules.initialsizing_loading as loadingdiagram


import inputs.concept_1 as conc1
import inputs.performance_inputs as perf
import inputs.constants as const
#print('done')
from Output.connection_departments import CL_TO, CL_cruise1, CL_cruise2, CL_cruise3, CL_land, V_climb1, V_climb2,V_climb3, f1,f2,f3, CD0_1, CD0_2, CD0_3

config1_plotdiagram=loadingdiagram.plot_loadingdiagram(perf.Sland*const.m_to_ft,Cl_TO,Cl_cruise1,Cl_land,V_climb1,perf.c,f1,perf.sigma, perf.TOP, CD0_1,conc1.A,conc1.e,1000,7000,100)
#config2_plotdiagram=loadingdiagram.plot_loadingdiagram(perf.Sland*const.m_to_ft,detailedsizing.Cl_TO,detailedsizing.Cl_cruise2,detailedsizing.Cl_land,detailedsizing.V_climb2,perf.c,detailedsizing.f2,perf.sigma, perf.TOP, detailedsizing.CD0_2,conc1.A,conc1.e,1000,7000,100)
#config3_plotdiagram=loadingdiagram.plot_loadingdiagram(perf.Sland*const.m_to_ft,detailedsizing.Cl_TO,detailedsizing.Cl_cruise3,detailedsizing.Cl_land,detailedsizing.V_climb3,perf.c,detailedsizing.f3,perf.sigma, perf.TOP, detailedsizing.CD0_3,conc1.A,conc1.e,1000,7000,100)