# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:29:12 2019

@author: Lisa
"""

'plot the loading diagram FOR THE ITERATIONS'
import modules.initialsizing_loading as loadingdiagram

import inputs.concept_1 as conc1
import inputs.performance_inputs as perf
import inputs.constant as const

import outputs.connection_departments as detailedsizing

plotdiagram=loadingdiagram.plot_loadingdiagram((perf.Sland,detailedsizing.Cl_TO,deltailedsizing.Cl_cruise,detailedsizing.Cl_land,detailedsizing.V_climb,detailedsizing.c,detailedsizing.f,perf.sigma, perf.TOP, CD0,conc1.A,conc1.e,1000,7000,100)