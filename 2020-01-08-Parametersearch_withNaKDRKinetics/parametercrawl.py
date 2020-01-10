## exec(open('parametercrawl.py').read())
## randomly selects a free parameter value and then stores the model with the scores

import os
import sys
from neo.io import AxonIO
import matplotlib.pyplot as plt
import numpy as np
import moose
import rdesigneur as rd
import time
import pandas as pd
import plotexp
import MOOSEModel_4
import numpy.random as nr
import featuresv5 as fts

# Setting seed
seeed = nr.randint(2**32 - 1)
# seeed = 110043517
nr.seed(seeed)
print(seeed)
f = open("OutputModels/outputModels_dict_"+str(seeed)+".py","a+")
f.write("# exec(open('OutputModels/outputModels_dict_"+str(seeed)+".py').read())\n\n")
f.write("Models = {} \n\n")


# Importing the means and std
features_ess_WT = pd.read_csv('features_ess_WT.csv',sep='\t', index_col=0)
meanfeatures = features_ess_WT.mean()
stdfeatures = features_ess_WT.std()

#dummyModel
dummyModel = {'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.437332664019589e-10, 'Rm': 714010979.6046907, 'Em': -0.026818968916777236}, 'Channels': {'K_BK_Chan': {'Gbar': 1.1203559577421247e-09, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_BK_Chan_(Migliore2018)'}, 'Ca_T_Chan': {'Gbar': 1.747076151712841e-07, 'Erev':0.140, 'Kinetics': '../../Compilations/Kinetics/Ca_T_Chan_(Migliore2018)'}, 'Ca_L_Chan': {'Gbar': 2.8133433185169764e-07, 'Erev':0.140, 'Kinetics': '../../Compilations/Kinetics/Ca_L_Chan_(Migliore2018)'}, 'Ca_N_Chan': {'Gbar': 9.154289883517049e-08, 'Erev':0.140, 'Kinetics': '../../Compilations/Kinetics/Ca_N_Chan_(Migliore2018)'}, 'Na_Chan': {'Gbar': 1.1181775516241626e-05, 'Erev':0.092, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'gateX':[-0.038,0.0072,-0.0365,0.020,0.0161,0.0547,0.0311,0.00064], 'gateY':[-0.050,-0.004, -0.0456,0.00433,0.01198,0.0262,0.00854,0.039]}, 'Na_P_Chan': {'Gbar': 8.695807678039012e-09, 'Erev':0.092, 'Kinetics': '../../Compilations/Kinetics/Na_P_Chan_(Migliore2018)'}, 'K_DR_Chan': {'Gbar': 3.1760459929263415e-07, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'gateX':[0.013,0.0088, 0.0125,0.0173,0,0,0.0341,0.1022]}, 'K_D_Chan': {'Gbar': 2.595004779072101e-09, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_D_Chan_(Migliore2018)'}, 'K_A_Chan': {'Gbar': 3.2555273173827533e-06, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_A_Chan_(Migliore2018)_ghk'}, 'K_M_Chan': {'Gbar': 1.332926217655261e-07, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_(Migliore2018)'}, 'K_SK_Chan': {'Gbar': 5.901861633366117e-10, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_(Migliore2018)'}, 'h_Chan': {'Gbar': 3.0214985435717124e-08, 'Erev':-0.030, 'Kinetics': '../../Compilations/Kinetics/h_Chan_(Migliore2018)'}}, 'Ca_Conc': {'Ca_B': 1278970793.3584003, 'Ca_tau': 0.19983256449900522, 'Ca_base': 5e-05, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Common)'}}}

# Define constants not to be changed
F = 96485.3329
elecPlotDt = 0.00005
elecDt = 0.00005

# Generating base model
sys.stdout = open(os.devnull, 'w')
rdes = MOOSEModel_4.generateModel(dummyModel, 150e-12) #Dummy model to load the channel kinetics
rdes.buildModel()
sys.stdout = sys.__stdout__

# Actual for loop
num_of_iterations = 5000
mnum = 1
for i in range(num_of_iterations):
    print(f'iternum = {i}')
    sys.stdout = open(os.devnull, 'w')
    Modeltemp = {}
    sm_area = 14.2e-9 #Does not matter as a random value of Gbars is choosen anyway

    ##########Initializing each parameter randomly##############
    Modeltemp['Scores'] = {}
    Modeltemp['parameters']={}
    Modeltemp['parameters']['notes']=''
    Modeltemp['parameters']['Morphology']={}
    Modeltemp['parameters']['Morphology']['sm_len'] = moose.element('/model/elec/soma').length
    Modeltemp['parameters']['Morphology']['sm_diam'] = moose.element('/model/elec/soma').diameter
    Modeltemp['parameters']['Passive']={}
    Modeltemp['parameters']['Passive']['Cm'] = nr.uniform(80e-12, 200e-12)
    Modeltemp['parameters']['Passive']['Rm'] = nr.uniform(150e6, 1000e6)
    Modeltemp['parameters']['Passive']['Em'] = nr.uniform(-0.1, -0.05)
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
    Modeltemp['parameters']['Channels']['K_BK_Chan']['Gbar'] = nr.uniform(0.00001, 8)*sm_area
    Modeltemp['parameters']['Channels']['Ca_T_Chan']['Gbar'] = nr.uniform(0.00001, 1)*sm_area
    Modeltemp['parameters']['Channels']['Ca_L_Chan']['Gbar'] = nr.uniform(0.00001, 1)*sm_area
    Modeltemp['parameters']['Channels']['Ca_N_Chan']['Gbar'] = nr.uniform(0.00001, 1)*sm_area
    Modeltemp['parameters']['Channels']['Na_Chan']['Gbar'] = nr.uniform(7, 1000)*sm_area
    Modeltemp['parameters']['Channels']['Na_P_Chan']['Gbar'] = nr.uniform(0.0001, 1)*sm_area
    Modeltemp['parameters']['Channels']['K_DR_Chan']['Gbar'] = nr.uniform(0.00001, 1000)*sm_area
    Modeltemp['parameters']['Channels']['K_D_Chan']['Gbar'] = nr.uniform(0.00001, 0.5)*sm_area
    Modeltemp['parameters']['Channels']['K_A_Chan']['Gbar'] = nr.uniform(0.00001, 1000)*sm_area
    Modeltemp['parameters']['Channels']['K_M_Chan']['Gbar'] = nr.uniform(0.00001, 11)*sm_area
    Modeltemp['parameters']['Channels']['K_SK_Chan']['Gbar'] = nr.uniform(0.00001, 11)*sm_area
    Modeltemp['parameters']['Channels']['h_Chan']['Gbar'] = nr.uniform(0.00001, 2.5)*sm_area
    EK = nr.uniform(-0.100, -0.080)
    Eh = nr.uniform(-0.050, -0.030)
    ENa = nr.uniform(0.050, 0.100)
    ECa = nr.uniform(0.120, 0.140)
    Modeltemp['parameters']['Channels']['K_BK_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['Ca_T_Chan']['Erev'] = ECa
    Modeltemp['parameters']['Channels']['Ca_L_Chan']['Erev'] = ECa
    Modeltemp['parameters']['Channels']['Ca_N_Chan']['Erev'] = ECa
    Modeltemp['parameters']['Channels']['Na_Chan']['Erev'] = ENa
    Modeltemp['parameters']['Channels']['Na_P_Chan']['Erev'] = ENa
    Modeltemp['parameters']['Channels']['K_DR_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_D_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_A_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_M_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['K_SK_Chan']['Erev'] = EK
    Modeltemp['parameters']['Channels']['h_Chan']['Erev'] = Eh
    Modeltemp['parameters']['Channels']['Na_Chan']['gateX'] = [nr.uniform(-0.038,-0.024),nr.uniform(0.002,0.015), nr.uniform(-0.050,-0.017),nr.uniform(0.002,0.040),nr.uniform(0,0.1),nr.uniform(0,0.1),nr.uniform(0.010,0.050),nr.uniform(0.0004,0.002)]
    Modeltemp['parameters']['Channels']['Na_Chan']['gateY'] = [nr.uniform(-0.050,-0.030),nr.uniform(-0.010,-0.002), nr.uniform(-0.060,-0.035),nr.uniform(0.002,0.020),nr.uniform(0,0.1),nr.uniform(0,0.1),nr.uniform(0.002,0.030),nr.uniform(0.002,0.07)]
    Modeltemp['parameters']['Channels']['K_DR_Chan']['gateX'] = [nr.uniform(-0.012,0.030),nr.uniform(0.002,0.012), nr.uniform(-0.008,-0.030),nr.uniform(0.003,0.045),nr.uniform(0,0.1),nr.uniform(0,0.1),nr.uniform(0.005,0.055),nr.uniform(0.002,0.11)]
    #####################################################################

    Scores = fts.modelscore(Modeltemp, meanfeatures, stdfeatures, modelfeature=None)
    # print(Scores)
    # print('\n')
    if Scores==False:
        sys.stdout = sys.__stdout__
        continue

    # print(Modeltemp)

    #####################################################################
    Model = {}
    Model['Score']=Scores
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
    Model['parameters']['Channels']['Na_Chan']['gateX'] = Modeltemp['parameters']['Channels']['Na_Chan']['gateX']
    Model['parameters']['Channels']['Na_Chan']['gateY'] = Modeltemp['parameters']['Channels']['Na_Chan']['gateY']
    Model['parameters']['Channels']['K_DR_Chan']['gateX'] = Modeltemp['parameters']['Channels']['K_DR_Chan']['gateX']

    f.write("Models['Model" + str(mnum) +"'] = " + str(Model) + "\n\n")
    mnum = mnum+1
    sys.stdout = sys.__stdout__

f.close()
