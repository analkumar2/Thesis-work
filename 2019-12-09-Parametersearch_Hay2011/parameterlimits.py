# exec(open('parameterlimits.py').read())

import moose
import rdesigneur as rd
import matplotlib.pyplot as plt
import numpy as np
import pprint
import sys
import os

MinModel = {'Error': 5.786203961042588, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.6298427057258092e-10, 'Rm': 160962561.99524152, 'Em': -0.027511362706586523}, 'Channels': {'Ca_HVA_Chan': {'Gbar': 7.696749828221082e-08, 'Kinetics': '../../Compilations/Kinetics/Ca_HVA_Chan_(Hay2011)', 'Erev': 0.12130450054058048}, 'Ca_LVAst_Chan': {'Gbar': 2.873155173538973e-07, 'Kinetics': '../../Compilations/Kinetics/Ca_LVAst_Chan_(Hay2011)', 'Erev': 0.12130450054058048}, 'Na_T_Chan': {'Gbar': 0.000677225763028016, 'Kinetics': '../../Compilations/Kinetics/Na_T_Chan_(Hay2011)', 'Erev': 0.056771258417900934}, 'Na_P_Chan': {'Gbar': 6.988157506226489e-07, 'Kinetics': '../../Compilations/Kinetics/Na_P_Chan_(Hay2011)', 'Erev': 0.056771258417900934}, 'K_DR_Chan': {'Gbar': 9.312096069421991e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_Pst_Chan': {'Gbar': 3.475531746522014e-05, 'Kinetics': '../../Compilations/Kinetics/K_Pst_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_Tst_Chan': {'Gbar': 1.4625662314691967e-06, 'Kinetics': '../../Compilations/Kinetics/K_Tst_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_M_Chan': {'Gbar': 1.0917650985519803e-07, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_SK_Chan': {'Gbar': 7.838939751365516e-06, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'h_Chan': {'Gbar': 8.475782518152347e-09, 'Kinetics': '../../Compilations/Kinetics/h_Chan_(Hay2011)', 'Erev': -0.03980264614196384}}, 'Ca_Conc': {'Ca_B': 2380812661.109322, 'Ca_tau': 0.04757716032463525, 'Ca_base': 3.523548538143316e-05, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Hay2011)'}}}

MaxModel = {'Error': 5.786203961042588, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.6298427057258092e-10, 'Rm': 160962561.99524152, 'Em': -0.027511362706586523}, 'Channels': {'Ca_HVA_Chan': {'Gbar': 7.696749828221082e-08, 'Kinetics': '../../Compilations/Kinetics/Ca_HVA_Chan_(Hay2011)', 'Erev': 0.12130450054058048}, 'Ca_LVAst_Chan': {'Gbar': 2.873155173538973e-07, 'Kinetics': '../../Compilations/Kinetics/Ca_LVAst_Chan_(Hay2011)', 'Erev': 0.12130450054058048}, 'Na_T_Chan': {'Gbar': 0.000677225763028016, 'Kinetics': '../../Compilations/Kinetics/Na_T_Chan_(Hay2011)', 'Erev': 0.056771258417900934}, 'Na_P_Chan': {'Gbar': 6.988157506226489e-07, 'Kinetics': '../../Compilations/Kinetics/Na_P_Chan_(Hay2011)', 'Erev': 0.056771258417900934}, 'K_DR_Chan': {'Gbar': 9.312096069421991e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_Pst_Chan': {'Gbar': 3.475531746522014e-05, 'Kinetics': '../../Compilations/Kinetics/K_Pst_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_Tst_Chan': {'Gbar': 1.4625662314691967e-06, 'Kinetics': '../../Compilations/Kinetics/K_Tst_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_M_Chan': {'Gbar': 1.0917650985519803e-07, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'K_SK_Chan': {'Gbar': 7.838939751365516e-06, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_(Hay2011)', 'Erev': -0.09283990346130364}, 'h_Chan': {'Gbar': 8.475782518152347e-09, 'Kinetics': '../../Compilations/Kinetics/h_Chan_(Hay2011)', 'Erev': -0.03980264614196384}}, 'Ca_Conc': {'Ca_B': 2380812661.109322, 'Ca_tau': 0.04757716032463525, 'Ca_base': 3.523548538143316e-05, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Hay2011)'}}}


def minmaxModel(modeldict,minModel,maxModel):
    for keey in modeldict.keys():
        if isinstance(modeldict[keey], dict):
            a = minmaxModel(modeldict[keey],minModel[keey],maxModel[keey])
        elif isinstance(modeldict[keey], float):
            if modeldict[keey]<minModel[keey]:
                minModel[keey] = modeldict[keey]
            elif modeldict[keey]>maxModel[keey]:
                maxModel[keey] = modeldict[keey]
    return [minModel, maxModel]

for filee in os.listdir('Modelparameters_temp'):
    exec(open('Modelparameters_temp/'+filee).read())
    for modell in Models.keys():
        MinModel, MaxModel = minmaxModel(Models[modell],MinModel,MaxModel)
