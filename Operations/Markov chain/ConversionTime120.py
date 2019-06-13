#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:08:35 2019

@author: coenminderhoud
"""


from math import *
import random
def conversiontimeto120(x):
    N=0
    A=0.1 # chance of faliure of step 1.2
    B=0.11 # chance of faliure of step 2.10 t/m 2.12
    B1=0.01 #chance of damage fuselage 2.10 t/m 2.12
    B2=0.94 #chance of repeated damage 2.10 t/m 2.12
    C=0.05 # chance of failure of step 2.13
    C1=0.005 #chance of damage of step 2.13
    C2=0.96 #chance of repeated damage of step 2.13
    D=0.025 # chance of failure of step 2.14
    D1=0.005 # chance of damage of step 2.14
    D2=0.94 # chance of repeated damage of step 2.14
    E=0.025 # chance of failure of step 3.7
    E1=0.005 # chance of damage of step 3.7
    E2=0.94 # chance of repeated damage of step 3.7
    F=0.15 # chance of failure of step 3.8
    F1=0.05 # chance of damage of step 3.8
    F2=0.94 # chance of repeated damage of step 3.8
    G=0.05 # chance of failure of step 3.9
    G1=0.005 # chance of damage of step 3.9
    G2=0.96 # chance of repeated damage of step 3.9
    H=0.05 # chance of failure of step 3.10
    I=0.15 # chance of failure of step 3.11
    I1=0.05 # chance of damage of step 3.11
    I2=0.94 # chance of repeated damage of step 3.11
    J=0.05 # chance of failure of step 3.12
    J1=0.005 # chance of damage of step 3.12
    J2=0.96 # chance of repeated damage of step 3.12
    K=0.16 # chance of failure of step 3.13 t/m 3.16
    K1=0.015 # chance of damage of step 3.13 t/m 3.16
    K2=0.94 # chance of repeated damage of step 3.13 t/m 3.16    
    Delay1=15
    Delay2=15
    Delay3=15
    Delay4=15
    Delay5=15
    Delay6=15
    Delay7=15
    Delay8=15
    Delay9=15
    
    P=[]
    while N < x:
        Time=0
        # step 1.2    
        a = random.random()
        Time=Time+30
        while a<A:
            Time=Time+30
            a = random.random()
        # step 2.12 t/m 2.14
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
        # step 2.15
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
        # step 2.16
        d = random.random()
        Time=Time+15
        while d<D:
            Time=Time+15
            d = random.random() 
        while d>D and d<D+D1:
            Time=Time+Delay3
            d2 = random.random()
            d  = random.random()
            while d2<D2:
                Time=Time+Delay3
                d2 = random.random()
        # step 3.9
        e = random.random()
        Time=Time+15
        while e>E and e<E+E1:
            Time=Time+Delay4
            e2 = random.random()
            e  = random.random()
            while e2<E2:
                Time=Time+Delay4
                e2 = random.random()
        while e<E:
            Time=Time+15
            e = random.random()
        # step 3.10
        f = random.random()
        Time=Time+30
        while f>F and f<F+F1:
            Time=Time+Delay5
            f2 = random.random()
            f  = random.random()
            while f2<F2:
                Time=Time+Delay5
                f2 = random.random()
        while f<F:
            Time=Time+30
            f = random.random()
        # step 3.11
        g = random.random()
        Time=Time+15
        while g<G:
            Time=Time+15
            g = random.random() 
        while G>g and g<G+G1:
            Time=Time+Delay6
            g2 = random.random()
            g  = random.random()
            while g2<G2:
                Time=Time+Delay6
                g2 = random.random()
        # step 3.12
        h = random.random()
        Time=Time+15
        while h<H:
            Time=Time+15
            h = random.random()
        # step 3.13
        i = random.random()
        Time=Time+30
        while i>I and i<I+I1:
            Time=Time+Delay7
            i2 = random.random()
            i  = random.random()
            while i2<I2:
                Time=Time+Delay7
                i2 = random.random()
        while i<I:
            Time=Time+30
            i = random.random()
        # step 3.14
        j = random.random()
        Time=Time+15
        while g<G:
            Time=Time+15
            j = random.random() 
        while J>j and j<J+J1:
            Time=Time+Delay8
            j2 = random.random()
            j  = random.random()
            while j2<J2:
                Time=Time+Delay8
                j2 = random.random()
        # step 3.15
        k = random.random()
        Time=Time+15
        while k<K:
            Time=Time+15
            k = random.random() 
        while K>k and k<K+K1:
            Time=Time+Delay9
            k2 = random.random()
            k  = random.random()
            while k2<K2:
                Time=Time+Delay9
                k2 = random.random()                
        P.append(Time)
        N=N+1
    return P
