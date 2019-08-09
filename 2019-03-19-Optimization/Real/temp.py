#exec(open('Optimization/Custom/Real/temp.py').read())
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
            output25pA = outputV_trace
        elif q == 1:
            output150pA = outputV_trace
        elif q == 2:
            output300pA = outputV_trace
    return [output25pA, output150pA, output300pA]
    # outputV_trace = outputV_trace[20000:50000]
    # plt.plot(np.linspace(0,6,len(outputV_trace)),outputV_trace, label = 'Best fit')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Membrane potential (V)')
    # plt.show()
#
# inputfile = 'Optimization/Custom/Real/best_citizen.csv'
# best_citizen = []
# with open(inputfile) as file:
#     reader = csv.reader(file, delimiter = ',')
#     a = next(reader)
#     for row in reader:
#         best_citizen.append(row)
# best_citizen = np.array(best_citizen, dtype = np.float)
#
# inputfile = 'Optimization/Custom/Real/parmdone.csv'
# with open(inputfile) as file:
#     reader = csv.reader(file, delimiter = ',')
#     a = next(reader)
#     for row in reader:
#         objective = row
#         break
# objective = np.array(objective, dtype = np.float)

# try:
#     [moose.delete(x) for x in ['/model', '/library']]
#     # [moose.delete(x) for x in ['/model']]
# except:
#     pass
#
# rdes = rd.rdesigneur(
#     cellProto = [['somaProto', 'soma', 1, 1]],
#     stimList = [['soma', '1', '.', 'vclamp', '-0.065 + (t>0.1 && t<2.2) * 0.065' ]],
#     plotList = [['soma', '1', '.', 'Vm', 'Soma Membrane potential']],
# )
#
# rdes.buildModel()
# moose.reinit()
# moose.start(6)
# rdes.display()
