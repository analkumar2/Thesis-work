#exec(open('Somatic model/myown_SKproto.py').read())

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
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor


#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

#Myown
ChPo = 'Somatic model/ChannelProtos_Sri2015_base'
ChP = 'Somatic model/ChannelProtos_Sri2015_base'
# ChP = 'Somatic model/ChannelProtos_Combe2018'

F = 96485.3329

with open('Somatic model/4_61016_best.pkl', 'rb') as f:
    Em, RM, RA, CM, sm_diam, sm_len, \
    Ca_Basal, Ca_tau, Ca_B, Na_Gbar, K_DR_Gbar, K_A_Gbar, \
    K_M_Gbar, h_Gbar, Ca_T_Gbar, Ca_R_Gbar, \
    Ca_L_Gbar, Ca_N_Gbar, K_SK_Gbar, K_BK_Gbar = pickle.load(f)

Na_Gbar = 0
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
diameter = sm_diam
elecPlotDt = 0.00005
runtime = 3
freeprm_dict = {'Em':Em, 'RM':RM, 'CM':CM, 'Ca_Basal':Ca_Basal, 'Ca_tau':Ca_tau, 'Ca_B':Ca_B, 'Na_Gbar':Na_Gbar, 'K_DR_Gbar':K_DR_Gbar, 'K_A_Gbar':K_A_Gbar, 'K_M_Gbar':K_M_Gbar, 'h_Gbar':h_Gbar, 'Ca_T_Gbar':Ca_T_Gbar, 'Ca_R_Gbar':Ca_R_Gbar, 'Ca_L_Gbar':Ca_L_Gbar, 'Ca_N_Gbar':Ca_N_Gbar, 'K_SK_Gbar':K_SK_Gbar, 'K_BK_Gbar':K_BK_Gbar}
freeprm_list = list(freeprm_dict)

def makeModel(stepV = 0):
    rdes = rd.rdesigneur(
        elecPlotDt = 0.00005,
        # cellProto = [['somaProto', 'soma', 12.76e-6, 0.01e-6]],
        cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
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
            ['K_DR_chan', 'soma', 'Gbar', str(K_DR_Gbar)],
            ['K_A_chan', 'soma', 'Gbar', str(K_A_Gbar)],
            ['K_M_chan', 'soma', 'Gbar', str(K_M_Gbar)],
            ['h_chan', 'soma', 'Gbar', str(h_Gbar)],
            ['Ca_T_chan', 'soma', 'Gbar', str(Ca_T_Gbar)],
            ['Ca_R_chan', 'soma', 'Gbar', str(Ca_R_Gbar)],
            ['Ca_L_chan', 'soma', 'Gbar', str(Ca_L_Gbar)],
            ['Ca_N_chan', 'soma', 'Gbar', str(Ca_N_Gbar)],
            ['K_SK_chan', 'soma', 'Gbar', str(K_SK_Gbar)],
            ['K_BK_chan', 'soma', 'Gbar', str(K_BK_Gbar)],
            ['Ca_conc', 'soma', 'thick', '177.9e-6'],
        ],
        stimList = [
            ['soma', '1', '.', 'vclamp', f'-0.055 + (t>0.96885 && t<=1.06880) * -0.010 + (t>1.06880 && t<1.8688) * {stepV+0.055}' ],
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
            ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
            # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
            # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
            # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
            # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
            # ['soma', '1', 'KBK_chan', 'Gk', 'Soma KBK conductance'],
            # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
            # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
            # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
            # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

            # ['soma', '1', 'Na_chan', 'Ik', 'Soma Sodium current'],
            # ['soma', '1', 'KDR_Schan', 'Ik', 'Soma Kdr current'],
            # ['soma', '1', 'KA_Schan', 'Ik', 'Soma KA current'],
            # ['soma', '1', 'KM_chan', 'Ik', 'Soma KM current'],
            # ['soma', '1', 'h_chan', 'Ik', 'Soma h current'],
            # ['soma', '1', 'KSK_chan', 'Ik', 'Soma KSK current'],
            # ['soma', '1', 'KBK_chan', 'Ik', 'Soma KBK current'],
            # ['soma', '1', 'CaT_chan', 'Ik', 'Soma CaT current'],
            # ['soma', '1', 'CaR_Schan', 'Ik', 'Soma CaR_S current'],
            # ['soma', '1', 'CaL_Schan', 'Ik', 'Soma CaL_S current'],
            # ['soma', '1', 'CaN_Schan', 'Ik', 'Soma CaN_S current'],
        ],
    )
    return rdes

fig, axs = plt.subplots(2, 1, constrained_layout=False)
# fig2, axs2 = plt.subplots(1, 1, constrained_layout=False)

Itrace_col = []
for stepV in np.arange(-0.035, -0.025, 0.050):
    print(stepV)
    #Deleting any previous run of the model
    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        [moose.delete(x) for x in ['/model']]
    except:
        pass
    rdes = makeModel(stepV)
    text_trap = io.StringIO() # JUst to suppress the inbuilt terminal output
    sys.stdout = text_trap # JUst to suppress the inbuilt terminal output
    rdes.buildModel()
    sys.stdout = sys.__stdout__ # restoring terminal output
    try:
        moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).td = 0
    except:
        pass
    moose.element('/model/elec/soma/Ca_conc').CaBasal = Ca_Basal
    moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
    moose.element('/model/elec/soma/Ca_conc').B = Ca_B
    moose.reinit()
    moose.start( runtime )
    Vtrace = moose.element('/model/graphs/plot0' ).vector
    Itrace = moose.element('/model/graphs/plot1' ).vector
    # Nacurr = moose.element('/model/graphs/plot2' ).vector
    Itrace_col.append(Itrace)
    axs[0].plot(np.linspace(0,runtime,len(Itrace)), Itrace, label=f'{stepV:.3f}')
    axs[1].plot(np.linspace(0,runtime,len(Vtrace)), Vtrace, label=f'{stepV:.3f}')
    # axs2.plot(np.linspace(0,runtime,len(Nacurr)), Nacurr, label=f'{stepV:.3f}')
    # axs.plot(np.linspace(0,runtime,len(Itrace)), Itrace, label=f'{stepV:.3f}')

fig.suptitle('Mymodel pre-apamin SK protocol')
axs[0].set_title('Soma holding current')
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Current (A)')
axs[0].set_ylim(-1000e-12,3000e-12)
axs[0].legend()
axs[1].set_title('Soma membrane potential')
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Potential (V)')
axs[1].set_ylim(-0.070,0.040)
axs[1].legend()
# axs2.set_title('Channel currents')
# axs2.set_xlabel('Time (s)')
# axs2.set_ylabel('Current (A)')
# axs2.set_ylim(-1000e-12,3000e-12)
# axs2.legend()
# axs.set_title('Soma holding current')
# axs.set_xlabel('Time (s)')
# axs.set_ylabel('Current (A)')
# axs.set_ylim(-1000e-12,3000e-12)
# axs.legend()
# rdes.display()
plt.show()
