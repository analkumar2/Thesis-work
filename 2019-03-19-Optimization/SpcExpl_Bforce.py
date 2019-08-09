#Aim is to load the input file, then run the simulation, and calculate sq error cost. Then a file is written with all the parameters and the cost.
#exec(open('/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/SpcExpl_Bforce.py').read())


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import csv
import random
import sys
import time

inputfile = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/inputV_trace.csv'
inputV_trace = []
with open(inputfile) as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        inputV_trace.append(row)
inputV_trace = np.transpose(inputV_trace)
inputV_trace = inputV_trace[20000:50000].flatten().astype(np.float)



Em_list = np.linspace(-0.080, -0.065, 150) #-0.075
sm_area_list = np.linspace(0.01e-12,100000e-12,1000) #927e-12
RM_list = np.linspace(0.01,20,1000) #2
CM_list = np.linspace(0.001,0.2,1000) #0.01
Na_SGbar_list = np.linspace(1, 1000, 1000) #350
KDR_SGbar_list = np.linspace(1,500, 500) #150
KA_SGbar_list = np.linspace(0.1, 20, 200) #5
KMGbar_list = np.linspace(0.1, 30, 200) #10
hGbar_list =  np.linspace(0.01, 5, 500) #0.18
CaTGbar_list = np.linspace(0.01, 5, 500) #0.5
CaRGbar_list = np.linspace(0.1, 10, 500) #1
CaLGbar_list = np.linspace(0.1, 20, 200) #5
KsAHPGbar_list = np.linspace(0.01, 500, 500) #150
KmAHPGbar_list = np.linspace(0.1, 5000, 5000) #2475

parmdone = []
exec(open('/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/CA1_WT_withoutchange.py').read())
# moose.reinit()
# moose.start( 6 )

# outputV_trace = moose.element('/model/graphs/plot0').vector
# outputV_trace = outputV_trace[20000:50000]

# start = False
# sp_start = []
# sp_end = []
# for j in range(0,len(outputV_trace)):
    # if outputV_trace[j]>0 and start == False:
        # start = True
        # sp_start.append(j)
    # if outputV_trace[j]<0 and start == True:
        # start = False
        # sp_end.append(j)
# sp_peak = np.add(sp_start,sp_end)/2
# err = np.sum((outputV_trace - inputV_trace)**2)
# parmdone.append([Em, sm_area, RM, CM, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KsAHPGbar, KmAHPGbar, len(sp_peak), err])


someprob = []
sec = time.time()
for i in range(5000):
    try:
        Em = random.choice(Em_list)
        sm_area = random.choice(sm_area_list)
        RM = random.choice(RM_list)
        CM = random.choice(CM_list)
        Na_SGbar = random.choice(Na_SGbar_list)
        KDR_SGbar = random.choice(KDR_SGbar_list)
        KA_SGbar = random.choice(KA_SGbar_list)
        KMGbar = random.choice(KMGbar_list)
        hGbar =  random.choice(hGbar_list)
        CaTGbar = random.choice(CaTGbar_list)
        CaRGbar = random.choice(CaRGbar_list)
        CaLGbar = random.choice(CaLGbar_list)
        KsAHPGbar = random.choice(KsAHPGbar_list)
        KmAHPGbar = random.choice(KmAHPGbar_list)

        moose.element('/model/elec/soma').Em = Em
        moose.element('/model/elec/soma').Rm = RM/sm_area
        moose.element('/model/elec/soma').Cm = CM*sm_area
        moose.element('/model/elec/soma').length = sm_area/(sm_diam*np.pi)
        moose.element('/model/elec/soma/Na_Schan').Gbar = Na_SGbar*sm_area
        moose.element('/model/elec/soma/KDR_Schan').Gbar = KDR_SGbar*sm_area
        moose.element('/model/elec/soma/KA_Schan').Gbar = KA_SGbar*sm_area
        moose.element('/model/elec/soma/KM_chan').Gbar = KMGbar*sm_area
        moose.element('/model/elec/soma/h_chan').Gbar = hGbar*sm_area
        moose.element('/model/elec/soma/CaT_chan').Gbar = CaTGbar*sm_area
        moose.element('/model/elec/soma/CaR_Schan').Gbar = CaRGbar*sm_area
        moose.element('/model/elec/soma/CaL_Schan').Gbar = CaLGbar*sm_area
        moose.element('/model/elec/soma/KsAHP_chan').Gbar = KsAHPGbar*sm_area
        moose.element('/model/elec/soma/KmAHP_chan').Gbar = KmAHPGbar*sm_area
        moose.element('/model/graphs/plot0').vector = np.array([])

        moose.reinit()
        moose.start( 6 )

        outputV_trace = moose.element('/model/graphs/plot0').vector
        outputV_trace = outputV_trace[20000:50000]

        start = False
        sp_start = []
        sp_end = []
        for j in range(0,len(outputV_trace)):
            if outputV_trace[j]>0 and start == False:
                start = True
                sp_start.append(j)
            if outputV_trace[j]<0 and start == True:
                start = False
                sp_end.append(j)
        if len(sp_start) != len(sp_end):
            sp_end.append(50000)
        sp_peak = np.add(sp_start,sp_end)/2
        err = np.sum((outputV_trace - inputV_trace)**2)
        if len(sp_peak)<36 and len(sp_peak)>20:
            parmdone.append([Em, sm_area, RM, CM, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KsAHPGbar, KmAHPGbar, len(sp_peak), err])
        # parmdone.append([Em, sm_area, RM, CM, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KsAHPGbar, KmAHPGbar, len(sp_peak), err])
        print(i, end='\r')
        # rdes.display()
    except:
        someprob.append([Em, sm_area, RM, CM, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KsAHPGbar, KmAHPGbar])
        continue

outputfile = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/parmdone.csv'
with open(outputfile,'a') as file:
    writer = csv.writer(file)
    # writer.writerow(['Em', 'sm_area', 'RM', 'CM', 'Na_SGbar', 'KDR_SGbar', 'KA_SGbar', 'KMGbar', 'hGbar', 'CaTGbar', 'CaRGbar', 'CaLGbar', 'KsAHPGbar', 'KmAHPGbar', 'len(sp_peak)', 'err'])
    writer.writerows(parmdone)

file.close()
print(time.time()-sec)

def plotwparms(parms):
    exec(open('/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/CA1_WT_withoutchange.py').read())
    Em = parms[0]
    sm_area = parms[1]
    RM = parms[2]
    CM = parms[3]
    Na_SGbar = parms[4]
    KDR_SGbar = parms[5]
    KA_SGbar = parms[6]
    KMGbar = parms[7]
    hGbar =  parms[8]
    CaTGbar = parms[9]
    CaRGbar = parms[10]
    CaLGbar = parms[11]
    KsAHPGbar = parms[12]
    KmAHPGbar = parms[13]

    moose.element('/model/elec/soma').Em = Em
    moose.element('/model/elec/soma').Rm = RM/sm_area
    moose.element('/model/elec/soma').Cm = CM*sm_area
    moose.element('/model/elec/soma').length = sm_area/(sm_diam*np.pi)
    moose.element('/model/elec/soma/Na_Schan').Gbar = Na_SGbar*sm_area
    moose.element('/model/elec/soma/KDR_Schan').Gbar = KDR_SGbar*sm_area
    moose.element('/model/elec/soma/KA_Schan').Gbar = KA_SGbar*sm_area
    moose.element('/model/elec/soma/KM_chan').Gbar = KMGbar*sm_area
    moose.element('/model/elec/soma/h_chan').Gbar = hGbar*sm_area
    moose.element('/model/elec/soma/CaT_chan').Gbar = CaTGbar*sm_area
    moose.element('/model/elec/soma/CaR_Schan').Gbar = CaRGbar*sm_area
    moose.element('/model/elec/soma/CaL_Schan').Gbar = CaLGbar*sm_area
    moose.element('/model/elec/soma/KsAHP_chan').Gbar = KsAHPGbar*sm_area
    moose.element('/model/elec/soma/KmAHP_chan').Gbar = KmAHPGbar*sm_area
    moose.element('/model/graphs/plot0').vector = np.array([])

    moose.reinit()
    moose.start( 6 )

    outputV_trace = moose.element('/model/graphs/plot0').vector
    plt.plot(np.linspace(0,6,len(outputV_trace)),outputV_trace)
    plt.show()
