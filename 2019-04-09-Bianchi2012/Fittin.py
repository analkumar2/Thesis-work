#exec(open('/home/analkumar2/Thesis/Thesis work/Bianchi2002_WT_CA1/Fittin.py').read())

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import time

Em_list = -0.075
sm_area_list = 927e-12
# sm_area = 53750e-12
RM_list = 2
Cm_list = 120e-12
Na_SGbar_list = 220
KDR_SGbar_list = 50
KA_SGbar_list = 5
KMGbar_list = 9
hGbar_list =  0.18
CaTGbar_list = 0.5
CaRGbar_list = 0.5
CaLGbar_list = 1
KsAHPGbar_list = 150
KmAHPGbar_list = 10

j = 1
for RM in RM_list:
    exec(open('/home/analkumar2/Thesis/Thesis work/Bianchi2002_WT_CA1/CA1_WT.py').read())
    V_trace = moose.element('/model/graphs/plot0').vector
    start = False
    sp_start = []
    sp_end = []
    for i in range(0,len(V_trace)):
        if V_trace[i]>0 and start == False:
            start = True
            sp_start.append(V_trace[i])
        if V_trace[i]<0 and start == True:
            start = False
            sp_end.append(V_trace[i])
    sp_peak = np.add(sp_start,sp_end)/2
    #if len(sp_peak)>2:
       # ax =plt.subplot(1,1,1)
       # ax.plot(V_trace)
       # plt.draw()
        #plt.pause(0.01)
    ax = plt.subplot(1,1,1)
    ln, = ax.plot(V_trace)
    plt.draw()
    plt.pause(0.00001)
    ln.remove()

    print(j, end='\r')
    j = j+1
