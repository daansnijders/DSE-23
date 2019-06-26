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
path = r"D:\Bibliotheken\Documenten\OneDrive\Stijn\Documenten\TU\3th Year\DSE\DSE-23\Graphing\Moosha2.xlsx"

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

x = np.array(x,dtype = 'f8')
y = np.array(y,dtype = 'f8')

x_current = 2550
step = 100
index = 0

x_new = []
y_new = []

while x_current<5500:
    x_new.append(0)
    for i in range(len(x)):
        if x_current - step < x[i] < x_current + step:

            x_new[index] += y[i]
    index+=1
    x_current+=step
    y_new.append(x_current-0.5*step)
            
fig = plt.figure(figsize = (12,6))
ax = fig.add_subplot(111)
ax.bar(y_new, x_new,width = step-5, align = 'center')
#ax.set_xticks(range(50,160,10))


for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

ax.set(xlabel = x_label, ylabel = y_label, title = title)  


fig.savefig('barplot2.pdf')
