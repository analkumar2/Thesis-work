#exec(open('Optimization/Custom/Real/SpcExpl_Galgo_2.py').read())

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

#Fixed model parameters
celsius = 34
ENa = 0.050
EK = -0.080
EKDR = -0.077
Eh = -0.010
ECa = 0.140
F = 96485.3329
R_input = 152e6
E_rest = -0.063
Cm = 134e-12
CaConc_tau = 20e-3
a1 = 4.25e-6
a2 = 1.8e-4
a3 = 1.07e-1
a4 = 6.52e-4
a5 = 5.22e-3
a6 = 1.48e-4
a7 = 9.75e-6
a8 = 1.42e-2
a9 = 1.60e-5
a10 = 6.49e-6

def child_err(parms):
    sm_area = parms[0]
    sm_vol = parms[1]
    Na_SGbar = parms[2]
    KDR_SGbar = parms[3]
    KA_SGbar = parms[4]
    KMGbar = parms[5]
    hGbar =  parms[6]
    CaTGbar = parms[7]
    CaRGbar = parms[8]
    CaLGbar = parms[9]
    KSKGbar = parms[10]
    KBKGbar = parms[11]
    CaConc_B = 100000/2/F/0.05/18/sm_area

    #Derived model parameters
    sm_diam = 4*sm_vol/sm_area
    sm_len = sm_area**2/4/sm_vol/np.pi
    RM = np.maximum(1e-12,1/(1/(R_input*sm_area) - Na_SGbar*a1 - KDR_SGbar*a2 - hGbar*a3 - KA_SGbar*a4 - KMGbar*a5 - CaLGbar*a6 - CaTGbar*a7 - CaRGbar*a8 - KSKGbar*a9 - KBKGbar*a10))
    Em = (E_rest/(R_input*sm_area) - Na_SGbar*ENa*a1 - KDR_SGbar*EKDR*a2 - hGbar*Eh*a3 - KA_SGbar*EK*a4 - KMGbar*EK*a5 - CaLGbar*ECa*a6 - CaTGbar*ECa*a7 - CaRGbar*ECa*a8 - KSKGbar*EK*a9 - KBKGbar*EK*a10)/(1/(R_input*sm_area) - Na_SGbar*a1 - KDR_SGbar*a2 - hGbar*a3 - KA_SGbar*a4 - KMGbar*a5 - CaLGbar*a6 - CaTGbar*a7 - CaRGbar*a8 - KSKGbar*a9 - KBKGbar*a10)

    #ChannelProtos file
    ChP = 'Optimization/Custom/Real/ChannelProtos_Combe2018'

    err = 0
    for q in [0,1,2]:
        #Deleting any previous run of the model
        try:
            # [moose.delete(x) for x in ['/model', '/library']]
            [moose.delete(x) for x in ['/model']]
        except:
            pass

        if q == 0:
            inp_curr = 25e-12
        elif q == 1:
            inp_curr = 150e-12
        elif q == 2:
            inp_curr = 300e-12
        rdes = rd.rdesigneur(
            elecDt = 10e-6,
            elecPlotDt = 50e-6,
            cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
            chanProto = [
                [ChP+'.Na_SChan()', 'Na_Schan'],
                [ChP+'.KDR_SChan()', 'KDR_Schan'],
                [ChP+'.KA_SChan()', 'KA_Schan'],
                    [ChP+'.KM_Chan()', 'KM_chan'],
                [ChP+'.h_Chan()', 'h_chan'],
                [ChP+'.CaT_Chan()', 'CaT_chan'],
            [ChP+'.CaR_SChan()', 'CaR_Schan'],
                    [ChP+'.CaL_SChan()', 'CaL_Schan'],
            [ChP+'.KSK_Chan()', 'KSK_chan'],
                [ChP+'.KBK_Chan()', 'KBK_chan'],
                [ChP+'.Ca_Conc()', 'Ca_conc'],
            ],
            passiveDistrib = [
                ['soma', 'RM', str(RM), 'RA', '1.5', 'Cm', str(Cm), 'initVm', str(E_rest), 'Em', str(Em)],
            ],
            chanDistrib = [
                ['Na_Schan', 'soma', 'Gbar', str(Na_SGbar)],
                    ['KDR_Schan', 'soma', 'Gbar', str(KDR_SGbar)],
            ['KA_Schan', 'soma', 'Gbar', str(KA_SGbar)],
                ['KM_chan', 'soma', 'Gbar', str(KMGbar)],
                ['h_chan', 'soma', 'Gbar', str(hGbar)],
                ['CaT_chan', 'soma', 'Gbar', str(CaTGbar)],
                ['CaR_Schan', 'soma', 'Gbar', str(CaRGbar)],
                ['CaL_Schan', 'soma', 'Gbar', str(CaLGbar)],
                ['KBK_chan', 'soma', 'Gbar', str(KBKGbar)],
                ['KSK_chan', 'soma', 'Gbar', str(KSKGbar)],
                ['Ca_conc', 'soma', 'thick', '177.9e-6'],
            ],
            stimList = [
                ['soma', '1', '.', 'inject', f'(t>=3.0813 && t<=3.5813) ? {inp_curr} : 0'],
            ],
            plotList = [
                ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
            ],
        )
        rdes.buildModel()
        moose.element('/model/elec/soma/Ca_conc').B = CaConc_B
        moose.element('/model/elec/soma/Ca_conc').tau = CaConc_tau
        moose.reinit()
        moose.start( 4 )
        outputV_trace = moose.element('/model/graphs/plot0').vector
        outputV_trace = outputV_trace[60000:75000]
        if q == 0:
            err = err + np.sum(np.abs(outputV_trace - input25pA))
        elif q == 1:
            err = err + np.sum(np.abs(outputV_trace - input150pA))
        elif q == 2:
            err = err + np.sum(np.abs(outputV_trace - input300pA))
    return err


input25pA = []
with open('Optimization/Custom/Real/D225118_25pA_WT_preAp.csv') as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        input25pA.append(row)
input25pA = np.transpose(input25pA)
input25pA = input25pA.flatten().astype(np.float)[:15000]/1000

input150pA = []
with open('Optimization/Custom/Real/D225118_150pA_WT_preAp.csv') as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        input150pA.append(row)
input150pA = np.transpose(input150pA)
input150pA = input150pA.flatten().astype(np.float)[:15000]/1000

input300pA = []
with open('Optimization/Custom/Real/D225118_300pA_WT_preAp.csv') as file:
    reader = csv.reader(file, delimiter = ',')
    for row in reader:
        input300pA.append(row)
input300pA = np.transpose(input300pA)
input300pA = input300pA.flatten().astype(np.float)[:15000]/1000

sm_area_list = np.linspace(0.01e-12,100000e-12,1000) #1477.4e-12
sm_vol_list = np.linspace(1e-18, 100000e-18, 1000) #467.5e-18
Na_SGbar_list = np.linspace(1, 2000, 1000) #350
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
    tomut_list = random.choices([True, False], cum_weights=[mutprob,1], k = len(Paramlist)-1)
    for i in range(len(tomut_list)):
        if tomut_list[i] == True:
            child[i] = random.choice(Paramlist[i])
    return child


# exec(open('Optimization/Custom/Real/CA1_custom.py').read())

numpop = 1000
numparent = numpop/10
numgenerations = 50
mutprob = 0.2
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
    KSKGbar = random.choice(KSKGbar_list)
    KBKGbar = random.choice(KBKGbar_list)
    CaConc_B = 100000/2/F/0.05/18/sm_area

    #Derived model parameters
    sm_diam = 4*sm_vol/sm_area
    sm_len = sm_area**2/4/sm_vol/np.pi
    RM = np.maximum(1e-12,1/(1/(R_input*sm_area) - Na_SGbar*a1 - KDR_SGbar*a2 - hGbar*a3 - KA_SGbar*a4 - KMGbar*a5 - CaLGbar*a6 - CaTGbar*a7 - CaRGbar*a8 - KSKGbar*a9 - KBKGbar*a10))
    Em = (E_rest/(R_input*sm_area) - Na_SGbar*ENa*a1 - KDR_SGbar*EKDR*a2 - hGbar*Eh*a3 - KA_SGbar*EK*a4 - KMGbar*EK*a5 - CaLGbar*ECa*a6 - CaTGbar*ECa*a7 - CaRGbar*ECa*a8 - KSKGbar*EK*a9 - KBKGbar*EK*a10)/(1/(R_input*sm_area) - Na_SGbar*a1 - KDR_SGbar*a2 - hGbar*a3 - KA_SGbar*a4 - KMGbar*a5 - CaLGbar*a6 - CaTGbar*a7 - CaRGbar*a8 - KSKGbar*a9 - KBKGbar*a10)

    #ChannelProtos file
    ChP = 'Optimization/Custom/Real/ChannelProtos_Combe2018'

    err = 0
    for q in [0,1,2]:
        #Deleting any previous run of the model
        try:
            # [moose.delete(x) for x in ['/model', '/library']]
            [moose.delete(x) for x in ['/model']]
        except:
            pass

        if q == 0:
            inp_curr = 25e-12
        elif q == 1:
            inp_curr = 150e-12
        elif q == 2:
            inp_curr = 300e-12
        rdes = rd.rdesigneur(
            elecDt = 10e-6,
            elecPlotDt = 50e-6,
            cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
            chanProto = [
                [ChP+'.Na_SChan()', 'Na_Schan'],
                [ChP+'.KDR_SChan()', 'KDR_Schan'],
                [ChP+'.KA_SChan()', 'KA_Schan'],
                [ChP+'.KM_Chan()', 'KM_chan'],
                [ChP+'.h_Chan()', 'h_chan'],
                [ChP+'.CaT_Chan()', 'CaT_chan'],
                [ChP+'.CaR_SChan()', 'CaR_Schan'],
                [ChP+'.CaL_SChan()', 'CaL_Schan'],
                [ChP+'.KSK_Chan()', 'KSK_chan'],
                [ChP+'.KBK_Chan()', 'KBK_chan'],
                [ChP+'.Ca_Conc()', 'Ca_conc'],
            ],
            passiveDistrib = [
                ['soma', 'RM', str(RM), 'RA', '1.5', 'Cm', str(Cm), 'initVm', str(E_rest), 'Em', str(Em)],
            ],
            chanDistrib = [
                ['Na_Schan', 'soma', 'Gbar', str(Na_SGbar)],
                ['KDR_Schan', 'soma', 'Gbar', str(KDR_SGbar)],
                ['KA_Schan', 'soma', 'Gbar', str(KA_SGbar)],
                ['KM_chan', 'soma', 'Gbar', str(KMGbar)],
                ['h_chan', 'soma', 'Gbar', str(hGbar)],
                ['CaT_chan', 'soma', 'Gbar', str(CaTGbar)],
                ['CaR_Schan', 'soma', 'Gbar', str(CaRGbar)],
                ['CaL_Schan', 'soma', 'Gbar', str(CaLGbar)],
                ['KBK_chan', 'soma', 'Gbar', str(KBKGbar)],
                ['KSK_chan', 'soma', 'Gbar', str(KSKGbar)],
                ['Ca_conc', 'soma', 'thick', '177.9e-6'],
            ],
            stimList = [
                ['soma', '1', '.', 'inject', f'(t>=3.0813 && t<=3.5813) ? {inp_curr} : 0'],
            ],
            plotList = [
                ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
            ],
        )
        rdes.buildModel()
        moose.element('/model/elec/soma/Ca_conc').B = CaConc_B
        moose.element('/model/elec/soma/Ca_conc').tau = CaConc_tau
        moose.reinit()
        moose.start( 4 )
        outputV_trace = moose.element('/model/graphs/plot0').vector
        outputV_trace = outputV_trace[60000:75000]
        if q == 0:
            err = err + np.sum(np.abs(outputV_trace - input25pA))
        elif q == 1:
            err = err + np.sum(np.abs(outputV_trace - input150pA))
        elif q == 2:
            err = err + np.sum(np.abs(outputV_trace - input300pA))

    if len(population)==0:
        population = np.array([sm_area, sm_vol, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KSKGbar, KBKGbar, err])
    else:
        population = np.vstack((population, [sm_area, sm_vol, Na_SGbar, KDR_SGbar, KA_SGbar, KMGbar, hGbar, CaTGbar, CaRGbar, CaLGbar, KSKGbar, KBKGbar, err]))
    print('The shape of the present population is '+str(population.shape), end='\r')

# curr_pop_file = 'Optimization/Custom/Real/curr_pop.csv'
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
best_citizen = population[np.where(population[:,-1] == min_err)][0]
# inputfile = 'Optimization/Custom/Real/best_citizen.csv'
# best_citizen = []
# with open(inputfile) as file:
#     reader = csv.reader(file, delimiter = ',')
#     a = next(reader)
#     for row in reader:
#         best_citizen.append(row)
# best_citizen = np.array(best_citizen, dtype = np.float)

for i in range(numgenerations):
    ind = np.argsort(population[:,-1])
    population = population[ind][:int(numparent)]
    print('Generation '+str(i))
    while len(population)< numpop:
        nxtparent12 = np.random.choice(int(numparent), size = 2)
        parent12 = population[nxtparent12]
        child = [random.choice(pm) for pm in np.stack((parent12[0],parent12[1]), axis = 1)]
        child = mutate(child,mutprob)
        child[-1] = child_err(child)
        population = np.vstack((population,child))
        print('Member #'+str(len(population)), end='\r')
        # print('Member #'+len(population))
    min_err = np.min(population[:,-1])
    print('Least error in Generation '+str(i)+ ' is '+str(min_err))
    best_citizen = np.vstack((best_citizen,population[np.where(population[:,-1] == min_err)][0]))
    outputfile = 'Optimization/Custom/Real/curr_pop.csv'
    with open(outputfile,'w') as file:
        writer = csv.writer(file)
        writer.writerow(['sm_area', 'sm_vol', 'Na_SGbar', 'KDR_SGbar', 'KA_SGbar', 'KMGbar', 'hGbar', 'CaTGbar', 'CaRGbar', 'CaLGbar', 'KSKGbar', 'KBKGbar', 'err'])
        writer.writerows(population)
    file.close()

    outputfile = 'Optimization/Custom/Real/best_citizen.csv'
    with open(outputfile,'w') as file:
        writer = csv.writer(file)
        writer.writerow(['sm_area', 'sm_vol', 'Na_SGbar', 'KDR_SGbar', 'KA_SGbar', 'KMGbar', 'hGbar', 'CaTGbar', 'CaRGbar', 'CaLGbar', 'KSKGbar', 'KBKGbar', 'err'])
        writer.writerows(best_citizen)
    file.close()

print(time.time()-sec)
