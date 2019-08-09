#exec(open('Optimization/Custom/Real/SpcExpl_Galgo.py').read())

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

def child_sperr(parms):
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
    err = np.sum(np.abs(outputV_trace - inputV_trace))
    return [len(sp_peak),err]


input25pA = []
with open('Optimization/Custom/Real/D225118_25pA_WT_preAp.csv') as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        input25pA.append(row)
input25pA = np.transpose(input25pA)
input25pA = input25pA.flatten().astype(np.float)[:15000]

input150pA = []
with open('Optimization/Custom/Real/D225118_150pA_WT_preAp.csv') as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        input150pA.append(row)
input150pA = np.transpose(input150pA)
input150pA = input150pA.flatten().astype(np.float)[:15000]

input300pA = []
with open('Optimization/Custom/Real/D225118_300pA_WT_preAp.csv') as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        input300pA.append(row)
input300pA = np.transpose(input300pA)
input300pA = input300pA.flatten().astype(np.float)[:15000]

sm_area_list = np.linspace(0.01e-12,100000e-12,1000) #1477.4e-12
sm_vol_list = np.linspace(1e-18, 100000e-18, 1000) #467.5e-18
Na_SGbar_list = np.linspace(1, 1000, 1000) #350
KDR_SGbar_list = np.linspace(1,1000, 500) #150
KA_SGbar_list = np.linspace(0.1, 100, 200) #35
KMGbar_list = np.linspace(0.1, 50, 200) #10
hGbar_list =  np.linspace(0.001, 5, 500) #0.18
CaTGbar_list = np.linspace(0.01, 10, 500) #0.5
CaRGbar_list = np.linspace(0.1, 10, 500) #1
CaLGbar_list = np.linspace(0.01, 20, 500) #5
KSKGbar_list = np.linspace(0.01, 20, 500) #150
KBKGbar_list = np.linspace(0.1, 5000, 5000) #2475

Paramlist = [sm_area_list, sm_vol_list, Na_SGbar_list, KDR_SGbar_list, KA_SGbar_list, KMGbar_list, hGbar_list, CaTGbar_list, CaRGbar_list, CaLGbar_list, KSKGbar_list, KBKGbar_list]

def mutate(childd,mutprob):
    child = childd
    tomut_list = random.choices([True, False], cum_weights=[mutprob,1], k = len(Paramlist)-2)
    for i in range(len(tomut_list)):
        if tomut_list[i] == True:
            child[i] = random.choice(Paramlist[i])
    return child


exec(open('Optimization/Custom/Real/CA1_custom.py').read())

numpop = 1000
numparent = numpop/10
numgenerations = 50
mutprob = 0.25
best_citizen = []
population = []
sec = time.time()
while len(population)<numpop:
    sm_area = random.choice(sm_area_list)
    sm_vol = random.choice(sm_vol_list)
    Na_SGbar = random.choice(Na_SGbar_list)
    KDR_SGbar = random.choice(KDR_SGbar_list)
    KA_SGbar = random.choice(KA_SGbar_list)
    KMGbar = random.choice(KMGbar_list)
    hGbar =  random.choice(hGbar_list)
    CaTGbar = random.choice(CaTGbar_list)
    CaRGbar = random.choice(CaRGbar_list)
    CaLGbar = random.choice(CaLGbar_list)
    KSKGbar = random.choice(KsAHPGbar_list)
    KBKGbar = random.choice(KmAHPGbar_list)

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
    err = np.sum(np.abs(outputV_trace - inputV_trace))
    if len(population)==0:
        population = np.concatenate((population, [Em, sm_area, RM, CM, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KsAHPGbar, KmAHPGbar, len(sp_peak), err]))
    else:
        population = np.vstack((population, [Em, sm_area, RM, CM, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KsAHPGbar, KmAHPGbar, len(sp_peak), err]))
    print('The shape of the present population is '+str(population.shape), end='\r')
    # rdes.display()

# curr_pop_file = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/curr_pop.csv'
# with open(curr_pop_file) as file:
#     reader = csv.reader(file)
#     a = next(reader)
#     for row in reader:
#         population.append(row)
# population  = population[-1*numpop:]
# print(len(population))

print('Initial population done')
population = np.array(population, dtype=np.float)
min_err = np.min(population[:,-1])
print(min_err)
# best_citizen = population[np.where(population[:,15] == min_err)][0]
inputfile = 'Optimizatiom/Custom/Real/best_citizen.csv'
best_citizen = []
with open(inputfile) as file:
    reader = csv.reader(file, delimiter = ',')
    a = next(reader)
    for row in reader:
        best_citizen.append(row)
best_citizen = np.array(best_citizen, dtype = np.float)

for i in range(numgenerations):
    # prconst= numpop/np.sum(1/population[:,15])
    # nxtpopulation = np.random.choice(numpop,size = int(numparent), replace = False, p = prconst/numpop/population[:,15])
    # population = population[nxtpopulation]
    ind = np.argsort(population[:,15])
    population = population[ind][:int(numparent)]
    print('Generation '+str(i))
    while len(population)< numpop:
        nxtparent12 = np.random.choice(int(numparent), size = 2)
        parent12 = population[nxtparent12]
        child = [random.choice(pm) for pm in np.stack((parent12[0],parent12[1]), axis = 1)]
        child = mutate(child,mutprob)
        child[14],child[15] = child_sperr(child)
        population = np.vstack((population,child))
        print('Member #'+str(len(population)), end='\r')
        # print('Member #'+len(population))
    min_err = np.min(population[:,15])
    print('Least error in Generation '+str(i)+ ' is '+str(min_err))
    best_citizen = np.vstack((best_citizen,population[np.where(population[:,15] == min_err)][0]))
    outputfile = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/curr_pop.csv'
    with open(outputfile,'w') as file:
        writer = csv.writer(file)
        writer.writerow(['Em', 'sm_area', 'RM', 'CM', 'Na_SGbar', 'KDR_SGbar', 'KA_SGbar', 'KMGbar', 'hGbar', 'CaTGbar', 'CaRGbar', 'CaLGbar', 'KsAHPGbar', 'KmAHPGbar', 'len(sp_peak)', 'err'])
        writer.writerows(population)
    file.close()

    outputfile = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/best_citizen.csv'
    with open(outputfile,'w') as file:
        writer = csv.writer(file)
        writer.writerow(['Em', 'sm_area', 'RM', 'CM', 'Na_SGbar', 'KDR_SGbar', 'KA_SGbar', 'KMGbar', 'hGbar', 'CaTGbar', 'CaRGbar', 'CaLGbar', 'KsAHPGbar', 'KmAHPGbar', 'len(sp_peak)', 'err'])
        writer.writerows(best_citizen)
    file.close()

print(time.time()-sec)
