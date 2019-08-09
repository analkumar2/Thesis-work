#exec(open('/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/temp.py').read())
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
import numpy.random


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
    return outputV_trace[20000:50000]
    # outputV_trace = outputV_trace[20000:50000]
    # plt.plot(np.linspace(0,6,len(outputV_trace)),outputV_trace, label = 'Best fit')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Membrane potential (V)')
    # plt.show()

inputfile = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/best_citizen.csv'
best_citizen = []
with open(inputfile) as file:
    reader = csv.reader(file, delimiter = ',')
    a = next(reader)
    for row in reader:
        best_citizen.append(row)
best_citizen = np.array(best_citizen, dtype = np.float)

inputfile = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/parmdone.csv'
with open(inputfile) as file:
    reader = csv.reader(file, delimiter = ',')
    a = next(reader)
    for row in reader:
        objective = row
        break
objective = np.array(objective, dtype = np.float)
