#exec(open('Kinetics_comparisons/runModel_singleChan2_P1.py').read())
# A model will be run with unit area to know how a channel responds to different kinetic comparison protocol 1.

import io
import os
import sys
from neo.io import AxonIO
import quantities as pq
import csv
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import numpy as np
import warnings
import moose
import pickle
import rdesigneur as rd

ChPr1 =
ChPr2 =
ChN1 =
ChN2 = 

rdes = rd.rdesigneur(
    elecPlotDt = 0.00005,
    # cellProto = [['somaProto', 'soma', 12.76e-6, 0.01e-6]],
    cellProto = [['somaProto', 'soma', 0.5642, 0.5642]],
    chanProto = [
        [ChP+'.Na_Chan()', 'Na_chan'],
        [ChP+'.K_DR_Chan()', 'K_DR_chan'],
        [ChP+'.K_A_Chan()', 'K_A_chan'],
        [ChP+'.K_M_Chan()', 'K_M_chan'],
        [ChP+'.h_Chan()', 'h_chan'],
        [ChP+'.Ca_T_Chan()', 'Ca_T_chan'],
        [ChP+'.Ca_R_Chan()', 'Ca_R_chan'],
        [ChP+'.Ca_L_Chan()', 'Ca_L_chan'],
        [ChPo+'.Ca_N_Chan()', 'Ca_N_chan'],
        [ChP+'.K_BK_Chan()', 'K_BK_chan'],
        [ChP+'.K_SK_Chan()', 'K_SK_chan'],
        [ChP+'.Ca_Conc()', 'Ca_conc'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', '1.5', 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
    ],
    chanDistrib = [
        ['Na_chan', 'soma', 'Gbar', str(Na_Gbar)],
        # ['K_DR_chan', 'soma', 'Gbar', str(K_DR_Gbar)],
        ['K_A_chan', 'soma', 'Gbar', str(K_A_Gbar)],
        # ['K_M_chan', 'soma', 'Gbar', str(K_M_Gbar)],
        # ['h_chan', 'soma', 'Gbar', str(h_Gbar)],
        # ['Ca_T_chan', 'soma', 'Gbar', str(Ca_T_Gbar)],
        # ['Ca_R_chan', 'soma', 'Gbar', str(Ca_R_Gbar)],
        # ['Ca_L_chan', 'soma', 'Gbar', str(Ca_L_Gbar)],
        # ['Ca_N_chan', 'soma', 'Gbar', str(Ca_N_Gbar)],
        # ['K_SK_chan', 'soma', 'Gbar', str(K_SK_Gbar)],
        # ['K_BK_chan', 'soma', 'Gbar', str(K_BK_Gbar)],
        # ['Ca_conc', 'soma', 'thick', '177.9e-6'],
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<6.5) * 0.050' ],
        ['soma', '1', '.', 'inject', f'(t>={stim_start} && t<={stim_end}) ? {curr} : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'KBK_chan', 'Gk', 'Soma KBK conductance'],
        # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
        # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
        # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

        ['soma', '1', 'Na_chan', 'Ik', 'Soma Sodium current'],
        ['soma', '1', 'K_DR_chan', 'Ik', 'Soma Kdr current'],
        ['soma', '1', 'K_A_chan', 'Ik', 'Soma KA current'],
        ['soma', '1', 'K_M_chan', 'Ik', 'Soma KM current'],
        ['soma', '1', 'h_chan', 'Ik', 'Soma h current'],
        ['soma', '1', 'K_SK_chan', 'Ik', 'Soma KSK current'],
        ['soma', '1', 'K_BK_chan', 'Ik', 'Soma KBK current'],
        ['soma', '1', 'Ca_T_chan', 'Ik', 'Soma CaT current'],
        ['soma', '1', 'Ca_R_chan', 'Ik', 'Soma CaR_S current'],
        ['soma', '1', 'Ca_L_chan', 'Ik', 'Soma CaL_S current'],
        ['soma', '1', 'Ca_N_chan', 'Ik', 'Soma CaN_S current'],
    ],
)

exp_tracef(10)

rdes = makeModel()
# text_trap = io.StringIO() # JUst to suppress the inbuilt terminal output
# sys.stdout = text_trap # JUst to suppress the inbuilt terminal output
rdes.buildModel()
# sys.stdout = sys.__stdout__ # restoring terminal output
try:
    moose.element( '/model/elec/soma/vclamp' ).gain *= 0.1
except:
    pass
# moose.element('/model/elec/soma/Ca_conc').CaBasal = Ca_Basal
# moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
# moose.element('/model/elec/soma/Ca_conc').B = Ca_B
moose.reinit()
moose.start( runtime )
rdes.display()
