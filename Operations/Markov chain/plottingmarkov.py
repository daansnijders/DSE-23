# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 11:37:08 2019

@author: Niels
"""
from ConversionTime90 import conversiontimeto90
from ConversionTime120 import conversiontimeto120
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

timeallowed = 8*60
#plottinglist = conversiontimeto90(100000)
plottinglist = conversiontimeto120(1000000)
binsamount=1000

plt.suptitle('Conversion time 90 to 120 pax')
plt.subplot(221)
plt.grid()
n, x, _ = plt.hist(plottinglist, bins=binsamount, density=True, range=(min(plottinglist), max(plottinglist)))
plt.axis([min(x),timeallowed+100,0,max(n)+0.003])
plt.xlabel('Time required (min)')
plt.ylabel('Probability density')
plt.title('Conversion time histogram')

problist=[]
timelist=[]



for i in range(len(n)):
    if n[i]!=0:
        problist.append(n[i])
        timelist.append(x[i]+(x[2]-x[1])/2)
        #timelist.append(x[i])
        
plt.subplot(222)     
plt.plot(timelist,problist)
plt.grid()
plt.xlim([min(x),timeallowed+100])
plt.xlabel('Time required (min)')
plt.ylabel('Normalized probability density')
plt.title('Conversion time probability density distribution')

cummulative=[]      
for y in range(len(problist)):
    cummulative.append(sum(problist[0:y+1]))

plt.subplot(223)
plt.grid()
plt.plot(timelist,cummulative) 
plt.xlabel('Time required (min)')
plt.ylabel('Accumulated probability')  
plt.xlim([min(x),timeallowed+100]) 
plt.axvline(timeallowed, linewidth = 3, color = 'r')
plt.legend(['Conversion time', '8-hour requirement'])
plt.title('Conversion time accumulated probability distribution')

def timepercentage():
    faults=0
    for i in range (len(plottinglist)):
        if plottinglist[i]>timeallowed:
            faults=faults+1
    percentage = faults/len(plottinglist)*100
    return percentage



print (timepercentage())
#print (max(plottinglist))
#plottinglist.sort(reverse=True)
#print (plottinglist)

#newprob=[]
#newtime=[]
#timemoment = min(plottinglist)-10
#
#while timemoment<=max(plottinglist):
#    number =0
#    for i in range(len(plottinglist)):
#        if plottinglist[i]<timemoment:
#            number=number+1
#    percent = number/len(plottinglist)
#    newprob.append(percent)
#    newtime.append(timemoment)
#    timemoment = timemoment+1
#    #print (timemoment)
#    


    
#plt.subplot(224)    
#plt.plot(newtime,newprob)
#plt.axvline(timeallowed, linewidth = 3, color = 'r')
#plottinglist.sort(reverse=True)
#print (plottinglist)



        