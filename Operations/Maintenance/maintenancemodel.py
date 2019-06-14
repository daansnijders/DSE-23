# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 09:21:25 2019

@author: Niels
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
#maintenance modelling
def failure_rate(beta, n, t):
    return beta/n *(t/n)**(beta-1)
    
def failure_function(beta, n, t):
    return 1-(np.exp(-(t/n)**beta))

def reliability_function(beta,n,t):
    return (np.exp(-(t/n)**beta))

def pdf(beta, n, t):
    return (beta / n) * (t / n)**(beta - 1) * np.exp(-(t / n)**beta)


def N_failures(beta,n,t):
    return (t/n)**beta

def integral_function(beta, n,tp):
    f = lambda t:(beta / n) * (t / n)**(beta - 1) * np.exp(-(t / n)**beta)*t
    return quad(f,0,tp)[0]

    
beta    = 1.263389  
n       = 3867.6847 
Cp      = 20000
Cc      = Cp*30
Ci      = Cp/5

#start=0
#end=10000
#steps = 1000000

#x = np.linspace(start,end,steps)
#problist = pdf(beta,n,x)
#failure = failure_function(beta,n,x)
#reliable = reliability_function(beta,n,x)
#failurerate = failure_rate(beta,n,x)
#n_fault = N_failures(beta,n,x)

def constant_interval(Cp, Cc, t):
    TEClist=[]
    for i in range(len(t)):
        tp=t[i]
        TEC = (Cp + Cc *N_failures(beta,n,tp))/tp
        TEClist.append(TEC)
    return TEClist

def age_based_replacement(Cp, Cc, t):
    TEClist=[]
    for i in range(len(t)):
        tp=t[i]
        TEC = (Cp *reliability_function(beta,n,tp)+Cc*failure_function(beta,n,tp))/(tp*reliability_function(beta,n,tp)+integral_function(beta,n,tp))
        TEClist.append(TEC)
    return TEClist

def inspection_based_replacement(Cp,Cc,Ci,t):
    
    TEClist=[]
    for i in range(len(t)):
        tp=t[i]
        n_number=0
        exlist=[]
        for i in range(100):
            n_number=n_number+1
            expected = reliability_function(beta,n,n_number*tp)*((Ci) + Cp*failure_function(beta,n,tp))  
            exlist.append(expected)
        exlistsum=sum(exlist)+Cc*failure_function(beta, n, tp)
        ttff = (tp*reliability_function(beta,n,tp)+integral_function(beta,n,tp))
        TEC = exlistsum/ttff
        TEClist.append(TEC)
        
          
        
    return TEClist
    

tp = np.linspace(100,2000,4000)
TEC_ci = constant_interval(Cp,Cc,tp)                    #constant interval
TEC_ab = age_based_replacement(Cp,Cc,tp)                #age based
TEC_ib = inspection_based_replacement(Cp,Cc,Ci,tp)      #inspection based


y_limits = [100,180]
y_base = y_limits[0]

TEC_ci_per=[]
TEC_ab_per=[]
TEC_ib_per=[]
for x in range(len(TEC_ci)):
    TEC_ci_per.append(TEC_ci[x]/y_base*100-100)
    TEC_ab_per.append(TEC_ab[x]/y_base*100-100)
    TEC_ib_per.append(TEC_ib[x]/y_base*100-100)
    


ci_min = min(TEC_ci)
ab_min = min(TEC_ab)
ib_min = min(TEC_ib)

ci_min_index = TEC_ci.index(min(TEC_ci))
ab_min_index = TEC_ab.index(min(TEC_ab))
ib_min_index = TEC_ib.index(min(TEC_ib))


total_list = [ci_min, ab_min, ib_min]
total_list_index = [ci_min_index, ab_min_index, ib_min_index]
name_list  = ['Constant interval', 'Age based replacement', 'Inspection based replacement']

total_min = total_list.index(min(total_list))

print ('The most optimal maintenance strategy is: ', name_list[total_min], 'with time between repair of ', round(tp[(total_list_index)[total_min]],0) )


#plt.subplot(211)
#plt.plot(tp, TEC_ci, tp, TEC_ab, tp, TEC_ib)
#plt.legend(['Constant interval', 'Age based replacement', 'Inspection based replacement'])
#plt.title('Maintenance strategy')
#plt.xlabel('Time between two preventive repairs')
#plt.ylabel('Total expected cost per unit time')
#plt.ylim(y_limits)
#plt.show()

#plt.subplot(212)
plt.plot(tp,TEC_ci_per, tp,TEC_ab_per,tp,TEC_ib_per)
plt.legend(['Constant interval', 'Age based replacement', 'Inspection based replacement'])
plt.title('Maintenance strategy')
plt.xlabel('Time between two preventive repairs (hours)')
plt.ylabel('Total expected cost per unit time (%)')
plt.ylim([0,80])
plt.grid()
plt.show()




    
    