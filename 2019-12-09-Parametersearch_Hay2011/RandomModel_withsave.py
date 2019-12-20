#exec(open('RandomModel_withsave.py').read())
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

import xmltodict

import plotexp
import MOOSEModel_2
exec(open('Modelparameters/dummyModels.py').read())




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
rdes = MOOSEModel_2.generateModel(Models['Model1'], 25e-12) #Dummy model to load the channel kinetics
Model = {}


## Start iterating
for i in range(num_of_iterations):
    print(i,end='\r' )
    Modeltemp = {}
    sys.stdout = open(os.devnull, 'w')
    lenn = moose.element('/model/elec/soma').length
    diamm = moose.element('/model/elec/soma').diameter

    sm_area = np.pi*lenn*diamm #Does not matter as a random value of Gbars is choosen anyway

    ##########Initializing each parameter randomly##############
    Modeltemp['Error'] = 0
    Modeltemp['parameters']={}
    Modeltemp['parameters']['notes']=''
    Modeltemp['parameters']['Morphology']={}
    Modeltemp['parameters']['Morphology']['sm_len'] = moose.element('/model/elec/soma').length
    Modeltemp['parameters']['Morphology']['sm_diam'] = moose.element('/model/elec/soma').diameter
    Modeltemp['parameters']['Passive']={}
    Modeltemp['parameters']['Passive']['Cm'] = nr.uniform(80e-12, 200e-12)
    Modeltemp['parameters']['Passive']['Rm'] = nr.uniform(150e6, 1000e6)
    Modeltemp['parameters']['Passive']['Em'] = nr.uniform(-0.1, -0.02)
    Modeltemp['parameters']['Channels']={}
    for chan in moose.wildcardFind('/model/elec/soma/#[CLASS==HHChannel2D]') + moose.wildcardFind('/model/elec/soma/#[CLASS==ZombieHHChannel]'):
        Modeltemp['parameters']['Channels'][chan.name] = {}
        Modeltemp['parameters']['Channels'][chan.name]['Gbar'] = moose.element(f"/model/elec/soma/{chan.name}").Gbar
        Modeltemp['parameters']['Channels'][chan.name]['Erev'] = moose.element(f"/model/elec/soma/{chan.name}").Ek
        Modeltemp['parameters']['Channels'][chan.name]['Kinetics'] = [i for i in np.ravel(rdes.chanProtoList) if chan.name+'()' in i][0][::-1].partition('.')[2][::-1]
    Modeltemp['parameters']['Ca_Conc'] = {}
    Modeltemp['parameters']['Ca_Conc']['Ca_B'] = nr.uniform(1.8e7, 1.8e11)
    Modeltemp['parameters']['Ca_Conc']['Ca_tau'] = nr.uniform(0.001, 0.200)
    Modeltemp['parameters']['Ca_Conc']['Ca_base'] = nr.uniform(0.01e-3, 1e-3)
    Modeltemp['parameters']['Ca_Conc']['Kinetics'] = [i for i in np.ravel(rdes.chanProtoList) if 'Ca_Conc()' in i][0][::-1].partition('.')[2][::-1]
    # Modeltemp['parameters']['Channels']['K_BK_Chan']['Gbar'] = nr.uniform(0.00001, 8)*sm_area
    Modeltemp['parameters']['Channels']['Ca_HVA_Chan']['Gbar'] = nr.uniform(0.00001, 10)*sm_area
    Modeltemp['parameters']['Channels']['Ca_LVAst_Chan']['Gbar'] = nr.uniform(0.00001, 100)*sm_area
    Modeltemp['parameters']['Channels']['Na_T_Chan']['Gbar'] = nr.uniform(7, 50000)*sm_area
    Modeltemp['parameters']['Channels']['Na_P_Chan']['Gbar'] = nr.uniform(0.0001, 100)*sm_area
    Modeltemp['parameters']['Channels']['K_DR_Chan']['Gbar'] = nr.uniform(0.00001, 20000)*sm_area
    Modeltemp['parameters']['Channels']['K_Pst_Chan']['Gbar'] = nr.uniform(0.00001, 10000)*sm_area
    Modeltemp['parameters']['Channels']['K_Tst_Chan']['Gbar'] = nr.uniform(0.00001, 1000)*sm_area
    Modeltemp['parameters']['Channels']['K_M_Chan']['Gbar'] = nr.uniform(0.00001, 11)*sm_area
    Modeltemp['parameters']['Channels']['K_SK_Chan']['Gbar'] = nr.uniform(0.00001, 1000)*sm_area
    Modeltemp['parameters']['Channels']['h_Chan']['Gbar'] = nr.uniform(0.00001, 2.5)*sm_area
    EK = nr.uniform(-0.100, -0.080)
    Eh = nr.uniform(-0.050, -0.030)
    ENa = nr.uniform(0.050, 0.100)
    ECa = nr.uniform(0.120, 0.140)
    # Modeltemp['parameters']['Channels']['K_BK_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['Ca_HVA_Chan']['Erev'] = ECa
    Modeltemp['parameters']['Channels']['Ca_LVAst_Chan']['Erev'] = ECa
    Modeltemp['parameters']['Channels']['Na_T_Chan']['Erev'] = ENa
    Modeltemp['parameters']['Channels']['Na_P_Chan']['Erev'] = ENa
    Modeltemp['parameters']['Channels']['K_DR_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_Pst_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_Tst_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_M_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_SK_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['h_Chan']['Erev'] = Eh
    #####################################################################

    rdes = MOOSEModel_2.generateModel(Modeltemp, 25e-12)
    moose.reinit()
    moose.start( 2.5 )
    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector
    tvec25pA = tvec[np.logical_and(tvec>=expData['25pA'][0][0],tvec<=expData['25pA'][0][-1])]
    Vmvec25pA = Vmvec[np.logical_and(tvec>=expData['25pA'][0][0],tvec<=expData['25pA'][0][-1])]
    error25pA = np.sum((Vmvec25pA-expData['25pA'][1])**2)

    rdes = MOOSEModel_2.generateModel(Modeltemp, 50e-12)
    moose.reinit()
    moose.start( 2.5 )
    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector
    tvec50pA = tvec[np.logical_and(tvec>=expData['50pA'][0][0],tvec<=expData['50pA'][0][-1])]
    Vmvec50pA = Vmvec[np.logical_and(tvec>=expData['50pA'][0][0],tvec<=expData['50pA'][0][-1])]
    error50pA = np.sum((Vmvec50pA-expData['50pA'][1])**2)

    rdes = MOOSEModel_2.generateModel(Modeltemp, 150e-12)
    moose.reinit()
    moose.start( 2.5 )
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
        Model['parameters']['Channels'][chan.name]['Kinetics'] = [i for i in np.ravel(rdes.chanProtoList) if chan.name+'()' in i][0][::-1].partition('.')[2][::-1]
        Model['parameters']['Channels'][chan.name]['Erev'] = moose.element(f"/model/elec/soma/{chan.name}").Ek
    Model['parameters']['Ca_Conc'] = {}
    Model['parameters']['Ca_Conc']['Ca_B'] = moose.element("/model/elec/soma/Ca_conc").B
    Model['parameters']['Ca_Conc']['Ca_tau'] = moose.element("/model/elec/soma/Ca_conc").tau
    Model['parameters']['Ca_Conc']['Ca_base'] = moose.element("/model/elec/soma/Ca_conc").Ca_base
    Model['parameters']['Ca_Conc']['Kinetics'] = [i for i in np.ravel(rdes.chanProtoList) if 'Ca_Conc()' in i][0][::-1].partition('.')[2][::-1]

    sys.stdout = sys.__stdout__
    if Model['Error']<errorthreshold:
        print(Model)
        print('\n')
        f.write("Models['Model" + str(mnum) +"'] = " + str(Model) + "\n\n")
        mnum = mnum+1
        # plt.plot(tvec25pA, Vmvec25pA, label='Model25pA')
        # plt.plot(expData['25pA'][0], expData['25pA'][1], label='Exp25pA')
        # plt.plot(tvec50pA, Vmvec50pA, label='Model50pA')
        # plt.plot(expData['50pA'][0], expData['50pA'][1], label='Exp50pA')
        # plt.plot(tvec150pA, Vmvec150pA, label='Model150pA')
        # plt.plot(expData['150pA'][0], expData['150pA'][1], label='Exp150pA')
        # plt.legend()
        # plt.show()

f.close()
