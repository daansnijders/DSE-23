# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 09:49:52 2019

@author: Niels
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as s



def weib(x,n,a):
    return (a / n) * (x / n)**(a - 1) * np.exp(-(x / n)**a)

   
def deboarding():
    
    weibullist=[]
    xlist=[]
    a=2.23619
    b=19.2339
    
    x=np.linspace(0,40,100)
    for i in range(len(x)):
        if x[i]>0:
            weibullist.append(weib(x[i],b,a))
            xlist.append(x[i]+10)
        
        
    #plt.plot(xlist,weibullist)
    #plt.grid()
    #plt.show()
    
    return xlist, weibullist

def boarding():
    weibullist=[]
    xlist=[]
    a=2.29263
    b=19.5224
    x=np.linspace(0,50,100)
    for i in range(len(x)):
        weibullist.append(weib(x[i],b,a))
        xlist.append(x[i])
     
        
    #plt.plot(xlist,weibullist)
    #plt.grid()
    #plt.show()
    
    return xlist, weibullist

def cleaning():  #USE THIS
    weibullist=[]
    xlist=[]
    a=2.16
    b=6.76
    
    x=np.linspace(0,45,100)
    for i in range(len(x)):
        if x[i]>0:
            weibullist.append(weib(x[i],b,a))
            xlist.append(x[i]+5)
    #plt.plot(xlist,weibullist)
    return xlist, weibullist

def fuel():
    gammalist=[]
    xlist=[]
    x=np.linspace(0,48,100)
    for i in range(len(x)):
        if x[i]>0:
            gammalist.append(s.gamma.pdf(x[i],a=1.64, scale=9.12))
            xlist.append(x[i]+2)
    #plt.plot(xlist,gammalist)
    return xlist, gammalist


def catering():
    weibullist=[]
    xlist=[]
    a=2.18
    b=17.37
    x=np.linspace(0,50,100)
    for i in range(len(x)):
        weibullist.append(weib(x[i],b,a))
        xlist.append(x[i])
     
        
    #plt.plot(xlist,weibullist)
    #plt.grid()
    #plt.show()
    
    return xlist, weibullist
    
#print (catering())
#print (cleaning())    

def deboard_time(pax):
    xlist,weibullist = deboarding()
    paxspeed=[] 
    
    for i in range(len(xlist)): 
        if xlist[i]>10:
            paxspeed.append(pax/xlist[i])
        if xlist[i]<10: 
            paxspeed.append(0)
    
    #plt.plot(paxspeed,weibullist)
    #plt.show()
    
    return paxspeed,weibullist

def deboardingplot():
    paxspeed90,weibullist90 = deboard_time(90)
    paxspeed120,weibullist120 = deboard_time(120)
    plt.plot(paxspeed90,weibullist90,paxspeed120,weibullist120)
    plt.grid()
    plt.legend(['90 passenger configuration', '120 passenger configuration'])
    plt.xlabel('Time (min)')
    plt.ylabel('Probability density')
    plt.title('Deboarding time')
    plt.show()



def board_time(pax):
    xlist,weibullist = boarding()
    paxspeed=[]
    weibull=[]
    for i in range(len(xlist)):
        if xlist[i]>1:
            paxspeed.append(pax/xlist[i])
            weibull.append(weibullist[i])
    return paxspeed, weibull

def boardingplot():
    paxspeed90,weibullist90 = board_time(90)
    paxspeed120,weibullist120 = board_time(120)
    plt.plot(paxspeed90,weibullist90,paxspeed120,weibullist120)
    plt.grid()
    plt.legend(['90 passenger configuration', '120 passenger configuration'])
    plt.xlabel('Time (min)')
    plt.ylabel('Probability density')
    #plt.axis([0,15,0,0.06])
    plt.title('Boarding time')
    plt.show()



def total_turnaround_time(pax):
    totalprob=[]
    totaltime=[]
    deboardtime,deboardprob = deboard_time(pax)
    boardtime,boardprob = board_time(pax)
    fueltime,fuelprob = fuel()
    for i in range(len(deboardprob)):
        for k in range(len(boardprob)):
            for m in range(len(fuelprob)):
                totalprob.append(deboardprob[i]*boardprob[k]*fuelprob[m])
                totaltime.append(deboardtime[i]+boardtime[k]+fueltime[m])
                
    return totalprob,totaltime
#print (total_turnaround_time(90))
def histogram(pax):
    problist=[]
    timelist=[]
    b,a=total_turnaround_time(pax)
    plt.subplot(221)
    plt.title('Histrogram turn around time')
    plt.xlabel('Time (min)')
    plt.ylabel('Probability density')
    plt.xlim([6,100])
    n, x, _ = plt.hist(a, bins=200, density=True, range=(0, max(a)))
    plt.clf()
    for i in range(len(n)):
        if n[i]!=0:
            problist.append(n[i]/sum(n))
            timelist.append(x[i]+(x[2]-x[1])/2)
    
    cummulative=[]      
    for y in range(len(problist)):
        cummulative.append(sum(problist[0:y+1]))
   
    
    return timelist, problist, cummulative

#print (histogram(90))

def plotting():
    timelist90, problist90, cummulative90 = histogram(90)
    timelist120, problist120, cummulative120 = histogram(120)
    
    plt.subplot(121)
    plt.plot(timelist90,problist90, timelist120, problist120)
    plt.xlabel('Time (min)')
    plt.ylabel('Normalized probability density')
    plt.title('Probability density distribution')
    plt.legend(['90 passenger configuration', '120 passenger configuration'])
    plt.xlim([0,100])
    plt.grid()
    
    
    plt.subplot(122)
    plt.plot(timelist90,cummulative90, timelist120, cummulative120)
    plt.xlabel('Time (min)')
    plt.ylabel('Accumulated probability')
    plt.title('Accumulated probability distribution')
    plt.legend(['90 passenger configuration', '120 passenger configuration'])
    plt.xlim([0,100])
    plt.grid()
    plt.suptitle('Turn around performance')
    return                             
print (plotting())

