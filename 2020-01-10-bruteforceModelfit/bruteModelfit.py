## exec(open('bruteModelfit.py').read())

## First make a function which takes in the free parameters and then do brute_scifit on it
## Migliore2018 free klinetics

import time
import numpy.random as nr
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import os
import sys
import moose
import brute_curvefit as bcf
import MOOSEModel_4 as mm
import plotexp
import featuresv5 as fts
from pprint import *
import warnings
import pickle
warnings.simplefilter(action='ignore', category=FutureWarning)

# Define constants not to be changed
F = 96485.3329
elecPlotDt = 0.00005
elecDt = 0.00005

exec(open('Modelparameters/dummyModels.py').read())

features_ess_WT = pd.read_csv('features_ess_WT.csv',sep='\t', index_col=0)
meanfeatures = features_ess_WT.mean()
stdfeatures = features_ess_WT.std()

bn = 0

def tolist_modelforbrute(modelasdict):
    Cm = modelasdict['parameters']['Passive']['Cm']
    Rm = modelasdict['parameters']['Passive']['Rm']
    Em = modelasdict['parameters']['Passive']['Em']
    EK = modelasdict['parameters']['Channels']['K_DR_Chan']['Erev']
    ENa = modelasdict['parameters']['Channels']['Na_Chan']['Erev']
    ECa = modelasdict['parameters']['Channels']['Ca_T_Chan']['Erev']
    Eh = modelasdict['parameters']['Channels']['h_Chan']['Erev']

    K_BK_Chan_Gbar = modelasdict['parameters']['Channels']['K_BK_Chan']['Gbar']
    Ca_T_Chan_Gbar = modelasdict['parameters']['Channels']['Ca_T_Chan']['Gbar']
    Ca_L_Chan_Gbar = modelasdict['parameters']['Channels']['Ca_L_Chan']['Gbar']
    Ca_N_Chan_Gbar = modelasdict['parameters']['Channels']['Ca_N_Chan']['Gbar']
    Na_Chan_Gbar = modelasdict['parameters']['Channels']['Na_Chan']['Gbar']

    Na_Chan_X_vhalf = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][0]
    Na_Chan_X_slope = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][1]
    Na_Chan_X_A = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][2]
    Na_Chan_X_B = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][3]
    Na_Chan_X_C = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][4]
    Na_Chan_X_D = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][5]
    Na_Chan_X_E = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][6]
    Na_Chan_X_F = modelasdict['parameters']['Channels']['Na_Chan']['gateX'][7]

    Na_Chan_Y_vhalf = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][0]
    Na_Chan_Y_slope = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][1]
    Na_Chan_Y_A = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][2]
    Na_Chan_Y_B = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][3]
    Na_Chan_Y_C = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][4]
    Na_Chan_Y_D = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][5]
    Na_Chan_Y_E = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][6]
    Na_Chan_Y_F = modelasdict['parameters']['Channels']['Na_Chan']['gateY'][7]

    Na_P_Chan_Gbar = modelasdict['parameters']['Channels']['Na_P_Chan']['Gbar']
    K_DR_Chan_Gbar = modelasdict['parameters']['Channels']['K_DR_Chan']['Gbar']

    K_DR_Chan_X_vhalf = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][0]
    K_DR_Chan_X_slope = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][1]
    K_DR_Chan_X_A = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][2]
    K_DR_Chan_X_B = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][3]
    K_DR_Chan_X_C = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][4]
    K_DR_Chan_X_D = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][5]
    K_DR_Chan_X_E = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][6]
    K_DR_Chan_X_F = modelasdict['parameters']['Channels']['K_DR_Chan']['gateX'][7]

    K_D_Chan_Gbar = modelasdict['parameters']['Channels']['K_DR_Chan']['Gbar']
    K_A_Chan_Gbar = modelasdict['parameters']['Channels']['K_DR_Chan']['Gbar']
    K_M_Chan_Gbar = modelasdict['parameters']['Channels']['K_DR_Chan']['Gbar']
    K_SK_Chan_Gbar = modelasdict['parameters']['Channels']['K_DR_Chan']['Gbar']
    h_Chan_Gbar = modelasdict['parameters']['Channels']['K_DR_Chan']['Gbar']

    Ca_Conc_B = modelasdict['parameters']['Ca_Conc']['Ca_B']
    Ca_Conc_tau = modelasdict['parameters']['Ca_Conc']['Ca_tau']
    Ca_Conc_base = modelasdict['parameters']['Ca_Conc']['Ca_base']

    return [Cm, Rm, Em, EK, ENa, ECa, Eh,
K_BK_Chan_Gbar, Ca_T_Chan_Gbar, Ca_L_Chan_Gbar, Ca_N_Chan_Gbar, Na_Chan_Gbar,
Na_Chan_X_vhalf, Na_Chan_X_slope, Na_Chan_X_A, Na_Chan_X_B, Na_Chan_X_C, Na_Chan_X_D, Na_Chan_X_E, Na_Chan_X_F,
Na_Chan_Y_vhalf, Na_Chan_Y_slope, Na_Chan_Y_A, Na_Chan_Y_B, Na_Chan_Y_C, Na_Chan_Y_D, Na_Chan_Y_E, Na_Chan_Y_F, Na_P_Chan_Gbar, K_DR_Chan_Gbar,
K_DR_Chan_X_vhalf, K_DR_Chan_X_slope, K_DR_Chan_X_A, K_DR_Chan_X_B, K_DR_Chan_X_C, K_DR_Chan_X_D, K_DR_Chan_X_E, K_DR_Chan_X_F,
K_D_Chan_Gbar, K_A_Chan_Gbar, K_M_Chan_Gbar, K_SK_Chan_Gbar, h_Chan_Gbar,
Ca_Conc_B, Ca_Conc_tau, Ca_Conc_base]

def todict_modelforbrute(t, Cm, Rm, Em, EK, ENa, ECa, Eh,
K_BK_Chan_Gbar, Ca_T_Chan_Gbar, Ca_L_Chan_Gbar, Ca_N_Chan_Gbar, Na_Chan_Gbar,
Na_Chan_X_vhalf, Na_Chan_X_slope, Na_Chan_X_A, Na_Chan_X_B, Na_Chan_X_C, Na_Chan_X_D, Na_Chan_X_E, Na_Chan_X_F,
Na_Chan_Y_vhalf, Na_Chan_Y_slope, Na_Chan_Y_A, Na_Chan_Y_B, Na_Chan_Y_C, Na_Chan_Y_D, Na_Chan_Y_E, Na_Chan_Y_F, Na_P_Chan_Gbar, K_DR_Chan_Gbar,
K_DR_Chan_X_vhalf, K_DR_Chan_X_slope, K_DR_Chan_X_A, K_DR_Chan_X_B, K_DR_Chan_X_C, K_DR_Chan_X_D, K_DR_Chan_X_E, K_DR_Chan_X_F,
K_D_Chan_Gbar, K_A_Chan_Gbar, K_M_Chan_Gbar, K_SK_Chan_Gbar, h_Chan_Gbar,
Ca_Conc_B, Ca_Conc_tau, Ca_Conc_base):
    themodel = {'Error': 5.907811441753849, 'parameters': {'notes': 'Migliore2018 free kinetics', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': Cm, 'Rm': Rm, 'Em': Em},
    'Channels': {
    'K_BK_Chan': {'Gbar': K_BK_Chan_Gbar, 'Erev':EK, 'Kinetics': '../../Compilations/Kinetics/K_BK_Chan_(Migliore2018)'},
    'Ca_T_Chan': {'Gbar': Ca_T_Chan_Gbar, 'Erev':ECa, 'Kinetics': '../../Compilations/Kinetics/Ca_T_Chan_(Migliore2018)'},
    'Ca_L_Chan': {'Gbar': Ca_L_Chan_Gbar, 'Erev':ECa, 'Kinetics': '../../Compilations/Kinetics/Ca_L_Chan_(Migliore2018)'},
    'Ca_N_Chan': {'Gbar': Ca_N_Chan_Gbar, 'Erev':ECa, 'Kinetics': '../../Compilations/Kinetics/Ca_N_Chan_(Migliore2018)'},
    'Na_Chan': {'Gbar': Na_Chan_Gbar, 'Erev':ENa, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'gateX':[Na_Chan_X_vhalf,Na_Chan_X_slope,Na_Chan_X_A,Na_Chan_X_B,Na_Chan_X_C,Na_Chan_X_D,Na_Chan_X_E,Na_Chan_X_F], 'gateY':[Na_Chan_Y_vhalf,Na_Chan_Y_slope,Na_Chan_Y_A,Na_Chan_Y_B,Na_Chan_Y_C,Na_Chan_Y_D,Na_Chan_Y_E,Na_Chan_Y_F]},
    'Na_P_Chan': {'Gbar': Na_P_Chan_Gbar, 'Erev':ENa, 'Kinetics': '../../Compilations/Kinetics/Na_P_Chan_(Migliore2018)'},
    'K_DR_Chan': {'Gbar': K_DR_Chan_Gbar, 'Erev':EK, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'gateX':[K_DR_Chan_X_vhalf,K_DR_Chan_X_slope,K_DR_Chan_X_A,K_DR_Chan_X_B,K_DR_Chan_X_C,K_DR_Chan_X_D,K_DR_Chan_X_E,K_DR_Chan_X_F]},
    'K_D_Chan': {'Gbar': K_D_Chan_Gbar, 'Erev':EK, 'Kinetics': '../../Compilations/Kinetics/K_D_Chan_(Migliore2018)'},
    'K_A_Chan': {'Gbar': K_A_Chan_Gbar, 'Erev':EK, 'Kinetics': '../../Compilations/Kinetics/K_A_Chan_(Migliore2018)_ghk'},
    'K_M_Chan': {'Gbar': K_M_Chan_Gbar, 'Erev':EK, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_(Migliore2018)'},
    'K_SK_Chan': {'Gbar': K_SK_Chan_Gbar, 'Erev':EK, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_(Migliore2018)'},
    'h_Chan': {'Gbar': h_Chan_Gbar, 'Erev':Eh, 'Kinetics': '../../Compilations/Kinetics/h_Chan_(Migliore2018)'}},
    'Ca_Conc': {'Ca_B': Ca_Conc_B, 'Ca_tau': Ca_Conc_tau, 'Ca_base': Ca_Conc_base, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Common)'}}}
    return themodel


def modelforbrute(t, Cm, Rm, Em, EK, ENa, ECa, Eh,
K_BK_Chan_Gbar, Ca_T_Chan_Gbar, Ca_L_Chan_Gbar, Ca_N_Chan_Gbar, Na_Chan_Gbar,
Na_Chan_X_vhalf, Na_Chan_X_slope, Na_Chan_X_A, Na_Chan_X_B, Na_Chan_X_C, Na_Chan_X_D, Na_Chan_X_E, Na_Chan_X_F,
Na_Chan_Y_vhalf, Na_Chan_Y_slope, Na_Chan_Y_A, Na_Chan_Y_B, Na_Chan_Y_C, Na_Chan_Y_D, Na_Chan_Y_E, Na_Chan_Y_F, Na_P_Chan_Gbar, K_DR_Chan_Gbar,
K_DR_Chan_X_vhalf, K_DR_Chan_X_slope, K_DR_Chan_X_A, K_DR_Chan_X_B, K_DR_Chan_X_C, K_DR_Chan_X_D, K_DR_Chan_X_E, K_DR_Chan_X_F,
K_D_Chan_Gbar, K_A_Chan_Gbar, K_M_Chan_Gbar, K_SK_Chan_Gbar, h_Chan_Gbar,
Ca_Conc_B, Ca_Conc_tau, Ca_Conc_base):
    global bn
    bn = bn+1
    print(bn)
    themodel = todict_modelforbrute(t, Cm, Rm, Em, EK, ENa, ECa, Eh,
        K_BK_Chan_Gbar, Ca_T_Chan_Gbar, Ca_L_Chan_Gbar, Ca_N_Chan_Gbar, Na_Chan_Gbar,
        Na_Chan_X_vhalf, Na_Chan_X_slope, Na_Chan_X_A, Na_Chan_X_B, Na_Chan_X_C, Na_Chan_X_D, Na_Chan_X_E, Na_Chan_X_F,
        Na_Chan_Y_vhalf, Na_Chan_Y_slope, Na_Chan_Y_A, Na_Chan_Y_B, Na_Chan_Y_C, Na_Chan_Y_D, Na_Chan_Y_E, Na_Chan_Y_F, Na_P_Chan_Gbar, K_DR_Chan_Gbar,
        K_DR_Chan_X_vhalf, K_DR_Chan_X_slope, K_DR_Chan_X_A, K_DR_Chan_X_B, K_DR_Chan_X_C, K_DR_Chan_X_D, K_DR_Chan_X_E, K_DR_Chan_X_F,
        K_D_Chan_Gbar, K_A_Chan_Gbar, K_M_Chan_Gbar, K_SK_Chan_Gbar, h_Chan_Gbar,
        Ca_Conc_B, Ca_Conc_tau, Ca_Conc_base)
    # pprint(themodel)
    sys.stdout = open(os.devnull, 'w')
    Sccrore = fts.modelscore(themodel, meanfeature, stdfeature, modelfeature=None)
    if Sccrore == False:
        sys.stdout = sys.__stdout__
        return list(1000*np.ones(len(meanfeatures)))
    else:
        sys.stdout = sys.__stdout__
        print(themodel)
        ssscore = list(Sccrore.values())
        print(ssscore)
        return ssscore
        

sm_area = 14.2e-9
restrictparams = [
                [80e-12,150e6,-0.1,-0.100,0.050,0.120,-0.050,
                0.00001*sm_area,0.00001*sm_area,0.00001*sm_area,0.00001*sm_area,7*sm_area,
                -0.050,0.002,-0.050,0.002,0,0,0.010,0.0004,
                -0.050,-0.010,-0.060,0.002,0,0,0.002,0.002,
                0.0001*sm_area,0.00001*sm_area,
                -0.012,0.002,-0.030,0.003,0,0,0.005,0.002,
                0.00001*sm_area,0.00001*sm_area,0.00001*sm_area,0.00001*sm_area,0.00001*sm_area,
                1.8e7,0.001,0.01e-3],

                [200e-12,1000e6,-0.02,-0.080,0.100,0.140,-0.030,
                8*sm_area,1*sm_area,1*sm_area,1*sm_area,1000*sm_area,
                -0.024,0.015,-0.017,0.040,0.1,0.1,0.050,0.002,
                -0.030,-0.002,-0.035,0.020,0.1,0.1,0.030,0.1,
                1*sm_area,1000*sm_area,
                0.030,0.020,-0.008,0.045,0.1,0.1,0.055,0.4,
                0.5*sm_area,1000*sm_area,11*sm_area,11*sm_area,2.5*sm_area,
                1.8e11,0.200,1e-3]
                ]

# paramfitted,error = bcf.brute_scifit(modelforbrute, [0,1,2], np.zeros(len(meanfeatures)),  restrictparams, ntol = 10000, returnnfactor = 0.005, maxfev = 1000)

# mm.plotModel(todict_modelforbrute([0,1,2],*pfsci), 150e-12)


############################
seeed = nr.randint(2**32 - 1)
# seeed = 110043517
nr.seed(seeed)
print(seeed)
##############################
f = open("Modelparameters_temp/outputModels_dict_"+str(seeed)+".py","a+")
f.write("# exec(open('Modelparameters_temp/outputModels_dict_"+str(seeed)+".py').read())\n\n")
f.write("Models = {} \n\n")

expData = {}
expData['150pA'] = plotexp.expdata('Experimental recordings/cell 4 of 61016.abf', 150e-12)
expData['150pA'][0], expData['150pA'][1] = expData['150pA'][0][expData['150pA'][0]<1] + (1-0.1391), expData['150pA'][1][expData['150pA'][0]<1]

mnum = 0

paramsfitted,errors = bcf.bruteforce(modelforbrute, [0,1,2], np.zeros(len(meanfeatures)),  restrictparams, ntol = 5000, returnnfactor = 0.002, printerrors=True)
pmsf_sl = paramsfitted[errors<900]
errors_sl = errors[errors<900]
for ll in range(len(pmsf_sl)):
    mnum = mnum +1
    Modell_dict = todict_modelforbrute([0,1,2] *pmsf_sl[ll])
    Modell_dict['Error'] = errors_sl[ll]
    print(Modell_dict)
    print('\n')
    f.write("Models['Model" + str(mnum) +"'] = " + str(Modell_dict) + "\n\n")
    plt.plot(*mm.runModel(todict_modelforbrute([0,1,2],*Modell_dict), 150e-12), label='Model150pA')
    plt.plot(expData['150pA'][0], expData['150pA'][1], label='Exp150pA')
    plt.xlabel('Time (s)')
    plt.ylabel('Membrane Potential (V)')
    plt.legend()
    plt.savefig(f'Output_temp/{seeed}_Model{mnum}')
    plt.clf()

if len(pmsf_sl)>0:
    pfsci, errsci = bcf.scipy_fit(modelforbrute, [0,1,2], np.zeros(len(meanfeatures)), restrict= restrictparams, p0list=pmsf_sl, maxfev = 1000, printerrors=True)
    for ll in range(len(pfsci)):
        mnum = mnum +1
        Modell_dict = todict_modelforbrute([0,1,2] *pfsci[ll])
        Modell_dict['Error'] = errsci[ll]
        print(Modell_dict)
        print('\n')
        f.write("Models['Model" + str(mnum) +"'] = " + str(Modell_dict) + "\n\n")
        plt.plot(*mm.runModel(todict_modelforbrute([0,1,2],*Modell_dict), 150e-12), label='Model150pA')
        plt.plot(expData['150pA'][0], expData['150pA'][1], label='Exp150pA')
        plt.xlabel('Time (s)')
        plt.ylabel('Membrane Potential (V)')
        plt.legend()
        plt.savefig(f'Output_temp/{seeed}_Model{mnum}')
        plt.clf()

f.close()
