# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 10:20:45 2019

@author: ahmadmahmoud
"""
import numpy as np
import matplotlib.pyplot as plt
n = 1000000
xr = np.linspace(0,1,n)
#Engr = 5.914 * x**1.2 *(1-x)**2
#ME = 2.316 * (x/0.85)**1.5 * (1-(x/0.85))**2
#Supp = 0.1197 * x**0.5 * (1-x)**0.5
#TD = 10.11 * ((x-0.22)/0.45)**2.5 * (1-((x-0.22)/0.45))**2
#TF = 20.88 * ((x-0.27)/0.50)**2 * (1-((x-0.27)/0.50))**2
yr =[]
for x in xr:
    
    if x>=0 and x<0.22:
        y = 5.914 * x**1.2 *(1-x)**2 + 2.316 * (x/0.85)**1.5 * (1-(x/0.85))**2 + 0.1197 * x**0.5 * (1-x)**0.5
    elif x>=0.22 and x<0.27:
        y = 5.914 * x**1.2 *(1-x)**2 + 2.316 * (x/0.85)**1.5 * (1-(x/0.85))**2 + 0.1197 * x**0.5 * (1-x)**0.5 + 10.11 * ((x-0.22)/0.45)**2.5 * (1-((x-0.22)/0.45))**2
    elif x>=0.27 and x<0.67:
        y = 5.914 * x**1.2 *(1-x)**2 + 2.316 * (x/0.85)**1.5 * (1-(x/0.85))**2 + 0.1197 * x**0.5 * (1-x)**0.5 + 10.11 * ((x-0.22)/0.45)**2.5 * (1-((x-0.22)/0.45))**2 + 20.88 * ((x-0.27)/0.50)**2 * (1-((x-0.27)/0.50))**2
    elif x>=0.67 and x<0.77:
        y = 5.914 * x**1.2 *(1-x)**2 + 2.316 * (x/0.85)**1.5 * (1-(x/0.85))**2 + 0.1197 * x**0.5 * (1-x)**0.5 + 20.88 * ((x-0.27)/0.50)**2 * (1-((x-0.27)/0.50))**2
    elif x>=0.77 and x<0.85:
        y = 5.914 * x**1.2 *(1-x)**2 + 2.316 * (x/0.85)**1.5 * (1-(x/0.85))**2 + 0.1197 * x**0.5 * (1-x)**0.5
    elif x>=0.85 and x<=1:
        y = 5.914 * x**1.2 *(1-x)**2 + 0.1197 * x**0.5 * (1-x)**0.5
    else:
        y=0
        print 0
    yr.append(y)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xr,yr)
ax.set(xlabel = 'normalised time',ylabel='normalised cost')
plt.show

#area =np.trapz(yr, xr)
#print (area)
a_tot =0
no_of_years = 8
n_per_year = n/no_of_years
for i in range(no_of_years):
    a = np.trapz(yr[(i*n_per_year): ((i+1)*n_per_year)], xr[(i*n_per_year): ((i+1)*n_per_year)])
    a_tot+= a
for i in range(no_of_years):
    a = np.trapz(yr[(i*n_per_year): ((i+1)*n_per_year)], xr[(i*n_per_year): ((i+1)*n_per_year)])
    print ("cost of year", i+1, "=", a/a_tot, "total cost")
    