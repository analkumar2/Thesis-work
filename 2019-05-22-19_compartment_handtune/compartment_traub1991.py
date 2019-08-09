# exec(open('manual_traubcompartment/compartment_traub1991.py').read())
#Author: Anal Kumar
import time
start_time = time.time()

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import pickle

#R is the radius of the compartment. B is the scaling term.
#length is the length of the compartment.
#The function finds the thickness of the Ca shell when the scaling
#term is given.
def findThickness(R, length, B):
    thick = R - np.sqrt(R**2 - 1.0/(2*B*F*np.pi*length))
    return(str(thick))

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

with open('6compartment/best_params.pkl', 'rb') as f:
    Em, RM, RA, CM, sm_diam, sm_len, \
    Ca_Basal, Ca_tau, Ca_B, Na_Gbar, K_DR_Gbar, K_A_Gbar, \
    K_M_Gbar, h_Gbar, Ca_T_Gbar, Ca_R_Gbar, \
    Ca_L_Gbar, Ca_N_Gbar, K_SK_Gbar, K_BK_Gbar = pickle.load(f)



ChP = '6compartment/ChannelProtos_Sri2015_base'
F = 96485.3329
diameter = sm_diam
elecPlotDt = 0.00005
preStimTime = 1
injectTime = 0.5
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 50e-12
freeprm_dict = {'Em':Em, 'RM':RM, 'CM':CM, 'Ca_Basal':Ca_Basal, 'Ca_tau':Ca_tau, 'Ca_B':Ca_B, 'Na_Gbar':Na_Gbar, 'K_DR_Gbar':K_DR_Gbar, 'K_A_Gbar':K_A_Gbar, 'K_M_Gbar':K_M_Gbar, 'h_Gbar':h_Gbar, 'Ca_T_Gbar':Ca_T_Gbar, 'Ca_R_Gbar':Ca_R_Gbar, 'Ca_L_Gbar':Ca_L_Gbar, 'Ca_N_Gbar':Ca_N_Gbar, 'K_SK_Gbar':K_SK_Gbar, 'K_BK_Gbar':K_BK_Gbar}
freeprm_list = list(freeprm_dict)

rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    cellProto = [
        ['6compartment/Compartments.swc','elec']
    ],

    chanProto = [
        [ChP+'.Na_Chan()', 'Na_chan'],
        [ChP+'.K_DR_Chan()', 'K_DR_chan'],
        [ChP+'.K_A_Chan()', 'K_A_chan'],
        [ChP+'.K_M_Chan()', 'K_M_chan'],
        [ChP+'.h_Chan()', 'h_chan'],
        [ChP+'.Ca_T_Chan()', 'Ca_T_chan'],
        [ChP+'.Ca_R_Chan()', 'Ca_R_chan'],
        [ChP+'.Ca_L_Chan()', 'Ca_L_chan'],
        [ChP+'.Ca_N_Chan()', 'Ca_N_chan'],
        [ChP+'.K_BK_Chan()', 'K_BK_chan'],
        [ChP+'.K_SK_Chan()', 'K_SK_chan'],
        [ChP+'.Ca_Conc()', 'Ca_conc'],
    ],

    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
        ['apical#', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
        ['dend#', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
    ],

    chanDistrib = [
        # soma
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

        # apical dendrtites
        ['Na_chan', 'apical#', 'Gbar', str(Na_Gbar)],
        ['K_DR_chan', 'apical#', 'Gbar', str(K_DR_Gbar)],
        ['K_A_chan', 'apical#', 'Gbar', str(K_A_Gbar)],
        ['K_M_chan', 'apical#', 'Gbar', str(K_M_Gbar)],
        ['h_chan', 'apical#', 'Gbar', str(h_Gbar)],
        ['Ca_T_chan', 'apical#', 'Gbar', str(Ca_T_Gbar)],
        ['Ca_R_chan', 'apical#', 'Gbar', str(Ca_R_Gbar)],
        ['Ca_L_chan', 'apical#', 'Gbar', str(Ca_L_Gbar)],
        ['Ca_N_chan', 'apical#', 'Gbar', str(Ca_N_Gbar)],
        ['K_SK_chan', 'apical#', 'Gbar', str(K_SK_Gbar)],
        ['K_BK_chan', 'apical#', 'Gbar', str(K_BK_Gbar)],
        ['Ca_conc', 'apical#', 'thick', '177.9e-6'],

        # dendrites
        ['Na_chan', 'dend#', 'Gbar', str(Na_Gbar)],
        ['K_DR_chan', 'dend#', 'Gbar', str(K_DR_Gbar)],
        ['K_A_chan', 'dend#', 'Gbar', str(K_A_Gbar)],
        ['K_M_chan', 'dend#', 'Gbar', str(K_M_Gbar)],
        ['h_chan', 'dend#', 'Gbar', str(h_Gbar)],
        ['Ca_T_chan', 'dend#', 'Gbar', str(Ca_T_Gbar)],
        ['Ca_R_chan', 'dend#', 'Gbar', str(Ca_R_Gbar)],
        ['Ca_L_chan', 'dend#', 'Gbar', str(Ca_L_Gbar)],
        ['Ca_N_chan', 'dend#', 'Gbar', str(Ca_N_Gbar)],
        ['K_SK_chan', 'dend#', 'Gbar', str(K_SK_Gbar)],
        ['K_BK_chan', 'dend#', 'Gbar', str(K_BK_Gbar)],
        ['Ca_conc', 'dend#', 'thick', '177.9e-6'],
    ],

    stimList = [
        # ['soma', '1', '.', 'vclamp', f'-0.065 + (t>{preStimTime} && t<{preStimTime+injectTime}) * 0.050' ],
        ['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {Injectcurr} : 0'],
    ],

    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        # ['soma', '1', '.', 'Im', 'Soma current'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'KBK_chan', 'Gk', 'Soma KBK conductance'],
        # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
        # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
        # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

        # ['soma', '1', 'Na_Schan', 'Ik', 'Soma Sodium current'],
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

        ['apical_1_8', '1', '.', 'Vm', 'apical_1_8 Membrane potential'],
        # ['apical_1_8', '1', '.', 'Im', 'apical_1_8 current'],
        ['dend_2_6', '1', '.', 'Vm', 'dend_2_6 Membrane potential'],
        # ['dend_2_6', '1', '.', 'Im', 'dend_2_6 current'],
    ],
)

rdes.buildModel()

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


rdes.display()
print(time.time() - start_time)
