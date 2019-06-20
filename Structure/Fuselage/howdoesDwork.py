# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 13:40:26 2019

@author: ahmadmahmoud
"""

import numpy as np
import matplotlib.pyplot as plt

D = np.linspace(0.01,0.2, 1000)
sigma = 350616614.36523831
delta = 0.0001
t = 0.0049
A = 0.00028812
b = 0.5 *t +0.25 *t + delta
Mom = sigma*A*(b/2)
MS = [0]*1000
for i in range(1000):
    
    bolt_d = D[i] - delta*2
    bolt_rupture = 1020000000
    I_bolt = np.pi/4 * (bolt_d/2)**4
    MS[i] = ((Mom*bolt_d)/(2*I_bolt))/bolt_rupture -1


plt.figure()
plt.plot(D, MS)
plt.show()
