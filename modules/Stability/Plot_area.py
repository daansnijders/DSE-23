# -*- coding: utf-8 -*-
"""
Created on Sat Jun 22 14:17:06 2019

@author: Daan
"""
plot = True
S = [132.21108287003085, 126.89626562809143,  120.83275487601703, 119.54654308369308, 117.61438984614605, 120.46981723317646]
S_v = [27.88510912244353, 27.230549857723517, 25.442543141882005, 25.33736252053676, 25.02741692395045, 26.47115748721427]
S_h = [17.194899048328345, 12.736971135874501, 11.804027075592344, 28.0 , 28.0, 29.0]
S_c = [35.69699237490833, 26.815402706195854, 26.982035800985738, 30.281128048680518 , 35.18387347778334, 37.13336136411162]
num = [1, 2, 3, 4, 5, 6]

if plot:
    fig = plt.figure()
    plt.title('Area of lifting surfaces development over iterations')
    ax1 = fig.add_subplot(111)
#            ax1.plot(x_cg_min1, x_le_MAC_range_perc)
#            ax1.plot(x_cg_max1, x_le_MAC_range_perc)s
    ax1.scatter(num, S)
    ax1.scatter(num, S_v)
    ax1.scatter(num, S_h)
    ax1.scatter(num, S_c)
    ax1.plot(num, S, label = "Wing area")
    ax1.plot(num, S_v, label = "Vertical tail area")
    ax1.plot(num, S_h, label = "Horizontal tail area")
    ax1.plot(num, S_c, label = "Canard area")
    ax1.set(ylim=[0., 140.],ylabel =  'Area $[m^2]$', xlabel = 'Iteration $[-]$')
    
#            ax1.scatter([f_min(y),f_max(y)],[y,y], color = 'b')
#            ax1.plot([f_min(y),f_max(y)],[y,y], color = 'b')

    ax1.legend()
    plt.show()