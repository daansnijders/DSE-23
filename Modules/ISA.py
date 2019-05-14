# author @R.J. Bouwmeester

from math import *

a = -0.0065
T0 = 288.15
g0 = 9.80665
p0 = 101325.
R = 287.
rho0 = 1.225

def isa_meters(h):
    rho = 0
    T = 0
    p = 0
    #check if within limits
    if h > 20000:
        print("This tool can calculate for up to 20000 meters.")
    elif h < 0:
        print("Invalid altitude.")
    else:
    #
    #calculations for 0 < h < 11000
    #
        if h < 11000:
            # print("Calculating...")

            #temperature
            T = T0 + a * h
            # print("Temperature at given altitude is",T,"degrees Kelvin.")
            #pressure
            exp = -g0 / a / R
            p1 = p0 * (T / T0)**exp
            # print("Pressure at given altitude is",p1,"Pa.")
            #density
            rho1 = rho0 * (T / T0)**(exp-1)
            rho = rho1
            p = p1
            # print("Density at given altitude is",rho1,"kg/m3.")
    #
    #calculations for 20000 > h > 11000
    #
        else:
            # print("Calculating")
            #temperature
            T = T0 + a * 11000
            # print("Temperature at given altitude is",T,"degrees Kelvin.")
            #pressure
            exp = -g0 / a / R
            p1 = p0 * (T / T0)**exp
            exp1 = -g0 / R / T * (h - 11000)
            p2 = p1 * e**exp1
            # print("Pressure at given altitude is",p2,"Pa.")
            #density
            rho1 = rho0 * (T / T0)**(exp-1)
            rho2 = rho1 * e **exp1
            rho = rho2
            p = p2
            # print("Density at given altitude is",rho2,"kg/m3.")
    return T, p, rho

#p->h
# if x==2:
#     #checking for invalid entries
#     from math import *
#     if p1 < 0.0000000001:
#         print "Invalid request."
#         continue
#     elif p1 > p0:
#         print "Too large a pressure."
#         continue
#     elif p1 < 22625.7914896:
#         exp = 1 / (-g0/a/R)
#         T = T0 * (p1 / p0)**exp
#         hm = (T - T0)/a
#         hf = hm / 0.3048
#         hfl = hf / 4500
#         print "Altitude for given pressure is",hm,"meters,",hf,"feet or",hfl,"FL."
#     else:
#         exp = 1 / (-g0/a/R)
#         T = T0 * (22625.7914896 / p0)**exp
#         hpre = (T - T0)/a
#         hm = hpre + log(p1/p0)*R*T/(-g0)
#         hf = hm / 0.0348
#         hfl = hf / 4500
#         print "Altitude for given pressure is",hm,"meters,",hf,"feet or",hfl,"FL."
# dummy = raw_input("Press enter to continue")
#
