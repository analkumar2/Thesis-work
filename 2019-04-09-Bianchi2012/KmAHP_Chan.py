#exec(open('/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work/Bianchi2002_WT_CA1/KmAHP_Chan.py').read())

import numpy as np
import csv

tbA = "/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work/Bianchi2002_WT_CA1/KmAHP_Chan_tbA.csv"
tbB = "/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work/Bianchi2002_WT_CA1/KmAHP_Chan_tbB.csv"

SOMA_A = 3.32e-9
F = 96485.3329
R = 8.314
Temp = 307.15
dt = 0.05e-3
ENa = 0.050
EK = -0.080 #EK for KDR is set to -0.077
Eh = -0.010
ECa = 0.140
Em = -0.070
Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
dV = (Vmax-Vmin)/Vdivs
Camin = 0
Camax = 5e-3
Cadivs = 5000
dCa = (Camax-Camin)/Cadivs
celsius = 20
gkbar = 0.01e4
d1 =1
d2 = 1.5
k1 = 0.18
k2 = 0.011
bbar = 0.28e3
abar = 0.48e3
F_KC = F/1000.0

def alp(v,c):
    return c*abar/(c + k1*np.exp(-2*d1*F*v/R/(273.15 + celsius)))
def bet(v,c):
    return bbar/(1 + c/k2*np.exp(-2*d2*F*v/R/(273.15 + celsius)))
    
tableA = np.zeros([Vdivs+1,Cadivs+2])
tableB = np.zeros([Vdivs+1,Cadivs+2])
x = Vmin
for i in np.arange(Vdivs+1):
    tableA[i] = [alp(x,y) for y in np.arange(Camin,Camax+dCa,dCa)]
    tableB[i] = tableA[i] + [bet(x, y) for y in np.arange(Camin,Camax+dCa,dCa)]
    x = x + dV
    print(str(i), end='\r')

with open(tbA, mode='w') as file:
    writer = csv.writer(file)
    writer.writerows(tableA)
    
with open(tbB, mode='w') as file:
    writer = csv.writer(file)
    writer.writerows(tableB)
    
tableA = []
tableB = [] 
i=0   
with open(tbA) as file:
    reader = csv.reader(file)
    for row in reader:
        tableA.append(row)
        print(i, end='\r')
        i+=1

i=0        
with open(tbB) as file:
    reader = csv.reader(file)
    for row in reader:
        np.append(tableB,row,axis =0)
    
    
    
    
    
    
    
    
    
    
    