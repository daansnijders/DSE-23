#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 09:47:37 2019

@author: coenminderhoud
"""

from math import *
import random
def conversiontimeto90(x):
    N=0
    A=0.1 # chance of faliure of step 1.2
    B=0.16 # chance of faliure of step 2.2 t/m 2.5
    B1=0.015 #chance of damage fuselage 2.2 t/m 2.5
    B2=0.94 #chance of repeated damage 2.2 t/m 2.5
    C=0.05 # chance of failure of step 2.6
    C1=0.005 #chance of damage of step 2.6
    C2=0.96 #chance of repeated damage of step 2.6
    D=0.05 # chance of failure of step 2.7
    E=0.05 # change of failure of step 2.8
    E1=0.005 # chance of damage of step 2.8
    E2=0.96 # chance of repeated damage of step 2.8
    F=0.025 # chance of failure of step 2.9
    F1=0.005 # chance of damage of step 2.9
    F2=0.94 # chance of repeated damage of step 2.9
    G=0.15 # chance of faliure of step 3.2
    G1=0.05 #chance of damage fuselage 3.2
    G2=0.94 #chance of repeated damage 3.2
    H=0.05 # chance of failure of step 3.3
    H1=0.005 # chance of damage of step 3.3
    H2=0.96 # chance of repeated damage of step 3.3
    I=0.11 # chance of failure of step 3.4 t/m 3.6
    I1=0.01 # chance of damage of step 3.4 t/m 3.6
    I2=0.94 # chance of repeated damge of step 3.4 t/m 3.6
    Delay1=15
    Delay2=15
    Delay3=15
    Delay4=15
    Delay5=15
    Delay6=15
    Delay7=15
    P=[]
    while N < x:
        Time=0
        # step 1.2    
        a = random.random()
        Time=Time+30
        while a<A:
            Time=Time+30
            a = random.random()
        # step 2.2 t/m 2.5
        b = random.random()
        Time=Time+30
        while b>B and b<B+B1:
            Time=Time+Delay1
            b2 = random.random()
            b  = random.random()
            while b2<B2:
                Time=Time+Delay1
                b2 = random.random()
        while b<B:
            Time=Time+30
            b = random.random() 
        # step 2.6
        c = random.random()
        Time=Time+15
        while c<C:
            Time=Time+15
            c = random.random()
        while c>C and c<C+C1:
            Time=Time+Delay2
            c2 = random.random()
            c  = random.random()
            while c2<C2:
                Time=Time+Delay2
                c2 = random.random()
        # step 2.7
        d = random.random()
        Time=Time+15
        while d<D:
            Time=Time+15
            d = random.random()    
        # step 2.8
        e = random.random()
        Time=Time+15
        while e<E:
            Time=Time+15
            e = random.random()
        while e>E and e<E+E1:
            Time=Time+Delay3
            e2 = random.random()
            e  = random.random()
            while e2<E2:
                Time=Time+Delay3
                e2 = random.random()
        # step 2.9
        f = random.random()
        Time=Time+15
        while f<F:
            Time=Time+15
            f = random.random() 
        while f>F and f<F+F1:
            Time=Time+Delay4
            f2 = random.random()
            f  = random.random()
            while f2<F2:
                Time=Time+Delay4
                f2 = random.random()
        # step 3.2
        g = random.random()
        Time=Time+30
        while g>G and g<G+G1:
            Time=Time+Delay5
            g2 = random.random()
            g  = random.random()
            while g2<G2:
                Time=Time+Delay5
                g2 = random.random()
        while g<G:
            Time=Time+30
            g = random.random() 
        #step 3.3
        h = random.random()
        Time=Time+15
        while h<H:
            Time=Time+15
            h = random.random() 
        while h>H and h<H+H1:
            Time=Time+Delay6
            h2 = random.random()
            h  = random.random()
            while h2<H2:
                Time=Time+Delay6
                h2 = random.random()
        #step 3.4 t/m 3.6
        i = random.random()
        Time=Time+15
        while i<I:
            Time=Time+15
            i = random.random() 
        while i>I and i<I+I1:
            Time=Time+Delay7
            i2 = random.random()
            i  = random.random()
            while i2<I2:
                Time=Time+Delay7
                i2 = random.random()
        P.append(Time)
        N=N+1       
    return P
