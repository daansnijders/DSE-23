# -*- coding: utf-8 -*-
"""
Created on Sat May  4 18:52:16 2019

@author: Stijn
"""
# importing modules required
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

# Loading data
#path = sys.argv[1]
path = r"D:\Bibliotheken\Documenten\OneDrive\Stijn\Documenten\TU\3th Year\DSE\DSE-23\Graphing\Coen.xlsx"

def read_excel(path):
    df = pd.read_excel(path)
    return df

data = read_excel(path)

def length(data):
    i = 5
    while data.isnull().iloc[i,1] == False:
        i += 1
    return i
l = length(data)

# Data gathering
title = data.axes[1][1]
x_label = data.iloc[2,1]
y_label = data.iloc[2,2]
x = np.array(data.iloc[5:l,1])
y = np.array(data.iloc[5:l,2:])[:,0]
label = data.iloc[l+1,2:]
graph_type = data.iloc[l+2,1]
legend = data.iloc[l+3,1]
save_path = data.iloc[28,1]

# Building the figure
fig = plt.figure(figsize = (12,4))
ax = fig.add_subplot(111)

plt.gcf().subplots_adjust(bottom=0.15)

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

ax.plot(x,y*100)    
ax.set(xlabel = x_label, ylabel = y_label, title = title)  
#ax.set(xlim = (0,np.max(x)*1.1), ylim = (0,np.max(y)*1.1))  
#ax.grid()
#if legend != 'no':
#    ax.legend(loc = legend)
#
fig.savefig('Seasonalfluctuations.pdf')
