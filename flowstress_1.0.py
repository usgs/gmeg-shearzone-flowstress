%matplotlib nbagg

import matplotlib
import numpy as n
import matplotlib.pyplot as plt
import os
import math
from matplotlib import cm

 #Calculate differential stress
d = 9.2# grain size in microns
B = 2451 #2451 (holyoke 2010) #3631 (stipp 2003)
m = -1.26 # 
s = []# (MPa) Differential Stress

part = (math.log(d)-math.log(B))/m
s = math.exp(part)
s


#Define Constants
R = 8.3144598 #Gas constant
w = 30000 # (mm) shear zone width

#Define flow laws
FL = []#Citation, A, N, Q (in that order)
FL = [
      #['GT95wom', 0.00011, 4, 137000],
      ['KT84', 0.0000022, 2.7, 120000],
      ['GT95wm', 0.000000018, 4, 137000],
      ['J84', 0.00288, 1.8, 151000],
      ['K89', 0.0000011, 2.7, 134000],
      ['HC82', 0.0199, 1.8, 167000],
      ['LP92g', 0.000000066, 3.1, 135000],
      ['LP92a', 3.98E-10, 4, 135000],
      ['H01', 6.3E-12, 4, 135000],
      ['RB04',0.000012, 2.97, 242000]
      ]

#Set temperature range
Tc = []
for i in range(300, 602, 25):#min, max, step interval
    T = i
    Tc.append(T)
Tk = []
for i in range(len(Tc)): #Convert to Kelvin
    T = Tc[i]+273.15
    Tk.append(T)

#Calculate strain rate for all flow laws over temperature range, store in list
e = []
for j in range(0,len(FL)):
    f = []
    for i in range(0, len(Tk)):
        element = FL[j][1]*s**(FL[j][2])*math.exp(-FL[j][3]/(R*Tk[i]))#Flow law equation
        f.append(element)
        if len(f) == len(Tk):
            e.append(f)
            del f

#Include effect of water fugacity
F = 300 #MPa (Water Fugacity)
e = []
for j in range(0,len(FL)):
    f = []
    for i in range(0, len(Tk)):
        element = FL[j][1]*s**(FL[j][2])*F*math.exp(-FL[j][3]/(R*Tk[i]))#Flow law equation
        f.append(element)
        if len(f) == len(Tk):
            e.append(f)
            del f

#Calculate slip rate
V = []
for i in range(0, len(FL)):
    f = []
    for j in range(0, len(Tk)):
        element = w*31540000*e[i][j]#slip rate calculation
        f.append(element)
        if len(f) == len(Tk):
            V.append(f)
            del f

###PLOT Strain Rate
colors = iter(cm.jet(np.linspace(0, 1, len(e))))#set color ramp to jet with a line spacing of length e

plt.gca().set_yscale('log') #set y axis to log 
 
labels = []
for i in range(0,len(FL)):
    l = FL[i][0]
    labels.append(l)
    plt.plot(Tc, e[i], color=next(colors))
    

    
plt.xlabel('Temperature ( C)')
plt.ylabel('Strain Rate (1/s)')
plt.title('Flow Laws')
#for i in range(0, len(Tc)):
 #   plt.text(450,e[i][6],labels[i], fontsize=10, rotation = 8)
plt.legend(labels, fontsize=10, loc=4)
plt.show()


###PLOT Slip Rate
colors = iter(cm.jet(np.linspace(0, 1, len(V))))#set color ramp to jet with a line spacing of length e

plt.gca().set_yscale('log') #set y axis to log 
 
labels = []
for i in range(0,len(FL)):
    l = FL[i][0]
    labels.append(l)
    plt.plot(Tc, V[i], color=next(colors))
    
plt.xlabel('Temperature ( C)')
plt.ylabel('Velocity (mm/yr)')
plt.title('Slip Rates')
#for i in range(0, len(Tc)):
 #   plt.text(450,e[i][6],labels[i], fontsize=10, rotation = 8)
plt.legend(labels, fontsize=10, loc=4)
plt.show()


#Save figure to pdf
from matplotlib.backends.backend_pdf import PdfPages

with PdfPages('SlipRatesNEW.pdf') as pdf:#Enter file name
    pdf.savefig()  # saves the current figure into a pdf page

