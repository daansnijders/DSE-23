# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:44:29 2019

@author: Niels
"""
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


#def integrand(x,a,b):
#    return a*x**2+b

#a=2
#b=1
#I=quad(integrand,1,2,args=(a,b))
#print (I)

beta    = 1.263389
n       = 3867.6847 

def integrand(beta,n,t):
        return (beta / n) * (t / n)**(beta - 1) * np.exp(-(t / n)**beta)*t

t=np.linspace(0,6000,12000)
lol=integrand(beta,n,t)


plt.plot(t,lol)
#plt.show()
def integral(tp):
    f = lambda t:(beta / n) * (t / n)**(beta - 1) * np.exp(-(t / n)**beta)*t
    return quad(f,0,tp)
print(integral(1000))