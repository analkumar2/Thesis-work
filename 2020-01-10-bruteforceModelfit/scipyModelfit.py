## exec(open('scipyModelfit.py').read())

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

basemodel = {'Score': {'E_rest_0': 2.609575477307922, 'Input resistance': 0.16018995544705936, 'Cell capacitance': 2.521294362849961, 'AP1_amp_1.5e-10': 4.389197291242428, 'APp_amp_1.5e-10': 4.275196185861353, 'AP1_width_1.5e-10': 0.5065074179381002, 'APp_width_1.5e-10': 0.1654346733227031, 'AP1_thresh_1.5e-10': 0.25560656877057536, 'APp_thresh_1.5e-10': 2.202315040562511, 'AP1_lat_1.5e-10': 2.4724313826673803, 'ISI1_1.5e-10': 0.16132614242429613, 'ISIl_1.5e-10': 2.127859892948867, 'ISIavg_1.5e-10': 1.7602257318590426, 'freq_1.5e-10': 1.812698652742249, 'Adptn_id_1.5e-10': 2.963431620290263, 'fAHP_AP1_amp_1.5e-10': 1.0636321666213355, 'fAHP_APp_amp_1.5e-10': 0.4159559904743107, 'mAHP_stimend_amp_1.5e-10': 0.6338812387481966, 'sAHP_stimend_amp_1.5e-10': 0.04270541826841748, 'AHP_AP1_amp_1.5e-10': 1.1560994266165598, 'AHP_APp_amp_1.5e-10': 0.1258756959177809, 'AHP_AP1_time_1.5e-10': 0.17982554970402462, 'AHP_APp_time_1.5e-10': 0.7911236815574301}, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.907197635658404e-10, 'Rm': 989580484.9444023, 'Em': -0.05499186805752422}, 'Channels': {'K_BK_Chan': {'Gbar': 4.957785758385209e-10, 'Kinetics': '../../Compilations/Kinetics/K_BK_Chan_(Migliore2018)', 'Erev': -0.08228279049680169}, 'Ca_T_Chan': {'Gbar': 1.4076676582853206e-08, 'Kinetics': '../../Compilations/Kinetics/Ca_T_Chan_(Migliore2018)', 'Erev': 0.13807288932591916}, 'Ca_L_Chan': {'Gbar': 5.605210670489204e-09, 'Kinetics': '../../Compilations/Kinetics/Ca_L_Chan_(Migliore2018)', 'Erev': 0.13807288932591916}, 'Ca_N_Chan': {'Gbar': 3.786164728035968e-09, 'Kinetics': '../../Compilations/Kinetics/Ca_N_Chan_(Migliore2018)', 'Erev': 0.13807288932591916}, 'Na_Chan': {'Gbar': 4.5829387167170905e-06, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.0650129639749113, 'gateX': [-0.030469927896272394, 0.011245916779753744, -0.02072235442245578, 0.023200553939386227, 0.042035248163026985, 0.07051816935011511, 0.03582580258960187, 0.0014073888892402904], 'gateY': [-0.03560080487592945, -0.006682508427258004, -0.035046606614666645, 0.0033272129485621845, 0.0014032138535265749, 0.03303149003773242, 0.022516474169238958, 0.05842349005988617]}, 'Na_P_Chan': {'Gbar': 1.362919938087954e-08, 'Kinetics': '../../Compilations/Kinetics/Na_P_Chan_(Migliore2018)', 'Erev': 0.0650129639749113}, 'K_DR_Chan': {'Gbar': 2.945160743565337e-06, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.08228279049680169, 'gateX': [0.0007569025217737148, 0.011714863982989343, -0.018130960644577608, 0.034260312432749664, 0.09253939509020503, 0.0021986566329906367, 0.031007112531068325, 0.004225666071668848]}, 'K_D_Chan': {'Gbar': 4.449110077189487e-09, 'Kinetics': '../../Compilations/Kinetics/K_D_Chan_(Migliore2018)', 'Erev': -0.08228279049680169}, 'K_A_Chan': {'Gbar': 6.373998833345491e-06, 'Kinetics': '../../Compilations/Kinetics/K_A_Chan_(Migliore2018)_ghk', 'Erev': -0.08228279049680169}, 'K_M_Chan': {'Gbar': 1.535404946211686e-08, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_(Migliore2018)', 'Erev': -0.08228279049680169}, 'K_SK_Chan': {'Gbar': 6.264286206901776e-09, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_(Migliore2018)', 'Erev': -0.08228279049680169}, 'h_Chan': {'Gbar': 9.341874528006847e-10, 'Kinetics': '../../Compilations/Kinetics/h_Chan_(Migliore2018)', 'Erev': -0.04183742501449533}}, 'Ca_Conc': {'Ca_B': 63363983214.667854, 'Ca_tau': 0.1104004093193713, 'Ca_base': 0.0009029917502226019, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Common)'}}}


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

x = [0,1,2]
func = modelforbrute
y = np.zeros(len(meanfeatures))
restrict = restrictparams
p0list = [tolist_modelforbrute(basemodel)]
bestModel  = bcf.scipy_fit(func, x, y,  restrict, p0list, maxfev = 1000, printerrors=True)