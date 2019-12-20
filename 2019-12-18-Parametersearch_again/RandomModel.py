#exec(open('RandomModel.py').read())
#Do check that the stimuli starting time is the same in the experimental cell and the model
#Before running this script, set the nr seed


import os
import sys
from neo.io import AxonIO
import matplotlib.pyplot as plt
import numpy as np
import moose
import rdesigneur as rd
import time
import numpy.random as nr
import xmltodict

import plotexp
import MOOSEModel

############################
num_of_iterations = 100
errorthreshold = 6
seeed = nr.randint(2**32 - 1)
nr.seed(seeed)
print(seeed)
##############################


# Define constants not to be changed
F = 96485.3329
elecPlotDt = 0.00005
elecDt = 0.00005

# List parameters and errors
Modellist = []

# Getting the experimental data
expData = {}
expData['25pA'] = plotexp.expdata('Experimental recordings/cell 4 of 61016.abf', 25e-12)
expData['25pA'][0], expData['25pA'][1] = expData['25pA'][0][expData['25pA'][0]<1] + (1-0.1391), expData['25pA'][1][expData['25pA'][0]<1]

expData['50pA'] = plotexp.expdata('Experimental recordings/cell 4 of 61016.abf', 50e-12)
expData['50pA'][0], expData['50pA'][1] = expData['50pA'][0][expData['50pA'][0]<1] + (1-0.1391), expData['50pA'][1][expData['50pA'][0]<1]

expData['150pA'] = plotexp.expdata('Experimental recordings/cell 4 of 61016.abf', 150e-12)
expData['150pA'][0], expData['150pA'][1] = expData['150pA'][0][expData['150pA'][0]<1] + (1-0.1391), expData['150pA'][1][expData['150pA'][0]<1]

# Generating base model
rdes = MOOSEModel.generateModel('Modelparameters/dummyModel.xml', 25e-12) #Dummy model to load the channel kinetics
Model = {}

## Start iterating
for i in range(num_of_iterations):
    print(i,end='\r' )
    sys.stdout = open(os.devnull, 'w')
    sm_area = 14.2e-9 #Does not matter as a random value of Gbars is choosen anyway
    moose.element('/model/elec/soma').Em = nr.uniform(-0.1, -0.02)
    moose.element('/model/elec/soma').Rm = nr.uniform(150e6, 1000e6)
    moose.element('/model/elec/soma').Cm = nr.uniform(80e-12, 200e-12)
    moose.element("/model/elec/soma/Na_Chan").Gbar = nr.uniform(7, 1000)*sm_area
    moose.element("/model/elec/soma/Na_P_Chan").Gbar = nr.uniform(0.0001, 1)*sm_area
    moose.element("/model/elec/soma/K_DR_Chan").Gbar = nr.uniform(0.00001, 1000)*sm_area
    moose.element("/model/elec/soma/K_D_Chan").Gbar = nr.uniform(0.00001, 0.5)*sm_area
    moose.element("/model/elec/soma/K_A_Chan").Gbar = nr.uniform(0.00001, 1000)*sm_area
    moose.element("/model/elec/soma/K_M_Chan").Gbar = nr.uniform(0.00001, 11)*sm_area
    moose.element("/model/elec/soma/K_SK_Chan").Gbar = nr.uniform(0.00001, 11)*sm_area
    moose.element("/model/elec/soma/K_BK_Chan").Gbar = nr.uniform(0.00001, 8)*sm_area
    moose.element("/model/elec/soma/h_Chan").Gbar = nr.uniform(0.00001, 2.5)*sm_area
    moose.element("/model/elec/soma/Ca_T_Chan").Gbar = nr.uniform(0.00001, 20)*sm_area
    moose.element("/model/elec/soma/Ca_L_Chan").Gbar = nr.uniform(0.00001, 20)*sm_area
    moose.element("/model/elec/soma/Ca_N_Chan").Gbar = nr.uniform(0.00001, 10)*sm_area
    moose.element("/model/elec/soma/Ca_conc").B = nr.uniform(1.8e7, 1.8e11)
    moose.element("/model/elec/soma/Ca_conc").tau = nr.uniform(0.001, 0.200)
    moose.element("/model/elec/soma/Ca_conc").Ca_base = nr.uniform(0.01e-3, 1e-3)
    moose.element("/model/elec/soma/Ca_conc").CaBasal = moose.element("/model/elec/soma/Ca_conc").Ca_base

    moose.element('/model/stims/stim0').expr = '(t>=1 && t<=1.5) ? 25e-12 : 0'
    moose.reinit()
    moose.start( 3 )
    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector
    tvec25pA = tvec[np.logical_and(tvec>=expData['25pA'][0][0],tvec<=expData['25pA'][0][-1])]
    Vmvec25pA = Vmvec[np.logical_and(tvec>=expData['25pA'][0][0],tvec<=expData['25pA'][0][-1])]
    error25pA = np.sum((Vmvec25pA-expData['25pA'][1])**2)

    moose.element('/model/stims/stim0').expr = '(t>=1 && t<=1.5) ? 50e-12 : 0'
    moose.reinit()
    moose.start( 3 )
    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector
    tvec50pA = tvec[np.logical_and(tvec>=expData['50pA'][0][0],tvec<=expData['50pA'][0][-1])]
    Vmvec50pA = Vmvec[np.logical_and(tvec>=expData['50pA'][0][0],tvec<=expData['50pA'][0][-1])]
    error50pA = np.sum((Vmvec50pA-expData['50pA'][1])**2)

    moose.element('/model/stims/stim0').expr = '(t>=1 && t<=1.5) ? 150e-12 : 0'
    moose.reinit()
    moose.start( 3 )
    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector
    tvec150pA = tvec[np.logical_and(tvec>=expData['150pA'][0][0],tvec<=expData['150pA'][0][-1])]
    Vmvec150pA = Vmvec[np.logical_and(tvec>=expData['150pA'][0][0],tvec<=expData['150pA'][0][-1])]
    error150pA = np.sum((Vmvec150pA-expData['150pA'][1])**2)

    # Calculate total error
    Model['Error'] = error25pA + error50pA + error150pA
    Model['parameters']={}
    Model['parameters']['notes']=''
    Model['parameters']['Morphology']={}
    Model['parameters']['Morphology']['sm_len'] = moose.element('/model/elec/soma').length
    Model['parameters']['Morphology']['sm_diam'] = moose.element('/model/elec/soma').diameter
    Model['parameters']['Passive']={}
    Model['parameters']['Passive']['Cm'] = moose.element('/model/elec/soma').Cm
    Model['parameters']['Passive']['Rm'] = moose.element('/model/elec/soma').Rm
    Model['parameters']['Passive']['Em'] = moose.element('/model/elec/soma').Em
    Model['parameters']['Channels']={}
    for chan in moose.wildcardFind('/model/elec/soma/#[CLASS==HHChannel2D]') + moose.wildcardFind('/model/elec/soma/#[CLASS==ZombieHHChannel]'):
        Model['parameters']['Channels'][chan.name] = {}
        Model['parameters']['Channels'][chan.name]['Gbar'] = moose.element(f"/model/elec/soma/{chan.name}").Gbar
        Model['parameters']['Channels'][chan.name]['Kinetics'] = [i for i in np.ravel(rdes.chanProtoList) if chan.name+'()' in i][0].split('.')[0]
    Model['parameters']['Ca_Conc'] = {}
    Model['parameters']['Ca_Conc']['Ca_B'] = moose.element("/model/elec/soma/Ca_conc").B
    Model['parameters']['Ca_Conc']['Ca_tau'] = moose.element("/model/elec/soma/Ca_conc").tau
    Model['parameters']['Ca_Conc']['Ca_base'] = moose.element("/model/elec/soma/Ca_conc").Ca_base
    Model['parameters']['Ca_Conc']['Kinetics'] = [i for i in np.ravel(rdes.chanProtoList) if 'Ca_Conc()' in i][0].split('.')[0]

    sys.stdout = sys.__stdout__
    if Model['Error']<errorthreshold:
        print(Model)
        print('\n')
        plt.plot(tvec25pA, Vmvec25pA, label='Model25pA')
        plt.plot(expData['25pA'][0], expData['25pA'][1], label='Exp25pA')
        plt.plot(tvec50pA, Vmvec50pA, label='Model50pA')
        plt.plot(expData['50pA'][0], expData['50pA'][1], label='Exp50pA')
        plt.plot(tvec150pA, Vmvec150pA, label='Model150pA')
        plt.plot(expData['150pA'][0], expData['150pA'][1], label='Exp150pA')
        plt.legend()
        plt.show()
