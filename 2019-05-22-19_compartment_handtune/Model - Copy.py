# exec(open('manual_traubcompartment/Model.py').read())

# Base model
# Author: Anal Kumar
import time
start_time = time.time()

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import pickle

exec(open('manual_traubcompartment/ChDensCA1pyr_Df.py').read()) # Dataframe with parameter values

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

F = 96485.3329
depth = 0.1 # No units. For ca_conc B
elecPlotDt = 0.00005
preStimTime = 1
injectTime = 0.5
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 150e-12

# ChanDistrib parameter value list
ChanDistrib_list = []
for ind in ChDensCA1pyr_Df.index:
    for coul in ['Na_Gbar', 'K_DR_Gbar', 'K_A_Gbar', 'K_M_Gbar', 'h_Gbar', 'Ca_T_Gbar', 'Ca_R_Gbar', 'Ca_L_Gbar', 'Ca_N_Gbar', 'K_SK_Gbar', 'K_BK_Gbar']:
        ChanDistrib_list.append([coul, ind, 'Gbar', str(ChDensCA1pyr_Df[coul][ind])])
    ChanDistrib_list.append(['Ca_conc', ind, 'Ca_Basal', str(ChDensCA1pyr_Df['Ca_Basal'][ind])])

# PassiveDistrib parameter value list
PassiveDistrib_list = []
for ind in ChDensCA1pyr_Df.index:
    PassiveDistrib_list.append([ind, 'RM', str(ChDensCA1pyr_Df['RM'][ind]), 'RA', str(ChDensCA1pyr_Df['RA'][ind]), 'CM', str(ChDensCA1pyr_Df['CM'][ind]), 'initVm', str(ChDensCA1pyr_Df['Em'][ind]), 'Em', str(ChDensCA1pyr_Df['Em'][ind])])


rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    cellProto = [
        ['manual_traubcompartment/Compartments.swc','elec']
    ],

    chanProto = [
        ['Channels/Na_Chan_(Sri2015base).Na_Chan()', 'Na_chan'],
        ['Channels/K_DR_Chan_(Sri2015base).K_DR_Chan()', 'K_DR_chan'],
        ['Channels/K_A_Chan_(Sri2015base).K_A_Chan()', 'K_A_chan'],
        ['Channels/K_M_Chan_(Sri2015base).K_M_Chan()', 'K_M_chan'],
        ['Channels/h_Chan_(Sri2015base).h_Chan()', 'h_chan'],
        ['Channels/Ca_T_Chan_(Sri2015base).Ca_T_Chan()', 'Ca_T_chan'],
        ['Channels/Ca_R_Chan_(Sri2015base).Ca_R_Chan()', 'Ca_R_chan'],
        ['Channels/Ca_L_Chan_(Sri2015base).Ca_L_Chan()', 'Ca_L_chan'],
        ['Channels/Ca_N_Chan_(Sri2015base).Ca_N_Chan()', 'Ca_N_chan'],
        ['Channels/K_BK_Chan_(Sri2015base).K_BK_Chan()', 'K_BK_chan'],
        ['Channels/K_SK_Chan_(Sri2015base).K_SK_Chan()', 'K_SK_chan'],
        ['Channels/Ca_Conc_(Common).Ca_Conc()', 'Ca_conc'],
    ],

    passiveDistrib = [
        PassiveDistrib_list,
        # ['soma', 'RM', str(ChDensCA1pyr_Df.RM[0]), 'RA', str(ChDensCA1pyr_Df.RA[0]), 'CM', str(ChDensCA1pyr_Df.CM[0]), 'initVm', str(ChDensCA1pyr_Df.Em[0]), 'Em', str(ChDensCA1pyr_Df.Em[0])],
        # ['apical_1#', 'RM', str(ChDensCA1pyr_Df.RM[1]), 'RA', str(ChDensCA1pyr_Df.RA[1]), 'CM', str(ChDensCA1pyr_Df.CM[1]), 'initVm', str(ChDensCA1pyr_Df.Em[1]), 'Em', str(ChDensCA1pyr_Df.Em[1])],
        # ['apical_e_1_9', 'RM', str(ChDensCA1pyr_Df.RM[10]), 'RA', str(ChDensCA1pyr_Df.RA[10]), 'CM', str(ChDensCA1pyr_Df.CM[10]), 'initVm', str(ChDensCA1pyr_Df.Em[10]), 'Em', str(ChDensCA1pyr_Df.Em[10])],
        # ['dend_2#', 'RM', str(RM[11]), 'RA', str(RA[11]), 'CM', str(CM[11]), 'initVm', str(Em[11]), 'Em', str(Em[11])],
        # ['dend_e_2_7', 'RM', str(RM[18]), 'RA', str(RA[18]), 'CM', str(CM[18]), 'initVm', str(Em[18]), 'Em', str(Em[18])],
    ],

    chanDistrib = [
        ChanDistrib_list,
        # # soma
        # ['Na_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.Na_Gbar[0])],
        # ['K_DR_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.K_DR_Gbar[0])],
        # ['K_A_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.K_A_Gbar[0])],
        # ['K_M_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.K_M_Gbar[0])],
        # ['h_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.h_Gbar[0])],
        # ['Ca_T_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.Ca_T_Gbar[0])],
        # ['Ca_R_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.Ca_R_Gbar[0])],
        # ['Ca_L_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.Ca_L_Gbar[0])],
        # ['Ca_N_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.Ca_N_Gbar[0])],
        # ['K_SK_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.K_SK_Gbar[0])],
        # ['K_BK_chan', 'soma', 'Gbar', str(ChDensCA1pyr_Df.K_BK_Gbar[0])],
        # ['Ca_conc', 'soma', 'Ca_Basal', str(ChDensCA1pyr_Df.Ca_Basal[0])],
        #
        # # apical dendrtites
        # ['Na_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.Na_Gbar[1])],
        # ['K_DR_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.K_DR_Gbar[1])],
        # ['K_A_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.K_A_Gbar[1])],
        # ['K_M_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.K_M_Gbar[1])],
        # ['h_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.h_Gbar[1])],
        # ['Ca_T_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.Ca_T_Gbar[1])],
        # ['Ca_R_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.Ca_R_Gbar[1])],
        # ['Ca_L_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.Ca_L_Gbar[1])],
        # ['Ca_N_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.Ca_N_Gbar[1])],
        # ['K_SK_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.K_SK_Gbar[1])],
        # ['K_BK_chan', 'apical_1#', 'Gbar', str(ChDensCA1pyr_Df.K_BK_Gbar[1])],
        # ['Ca_conc', 'apical_1#', 'Ca_Basal', str(ChDensCA1pyr_Df.Ca_Basal[0])],
        #
        # ['Na_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.Na_Gbar[10])],
        # ['K_DR_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.K_DR_Gbar[10])],
        # ['K_A_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.K_A_Gbar[10])],
        # ['K_M_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.K_M_Gbar[10])],
        # ['h_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.h_Gbar[10])],
        # ['Ca_T_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_T_Gbar[10])],
        # ['Ca_R_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_R_Gbar[10])],
        # ['Ca_L_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_L_Gbar[10])],
        # ['Ca_N_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_N_Gbar[10])],
        # ['K_SK_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.K_SK_Gbar[10])],
        # ['K_BK_chan', 'apical_e#', 'Gbar', str(ChDensCA1pyr_Df.K_BK_Gbar[10])],
        # ['Ca_conc', 'apical_e#', 'Ca_Basal', str(ChDensCA1pyr_Df.Ca_Basal[10])],
        #
        # # dendrites
        # ['Na_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.Na_Gbar[11])],
        # ['K_DR_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.K_DR_Gbar[11])],
        # ['K_A_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.K_A_Gbar[11])],
        # ['K_M_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.K_M_Gbar[11])],
        # ['h_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.h_Gbar[11])],
        # ['Ca_T_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.Ca_T_Gbar[11])],
        # ['Ca_R_chan', 'dend_2#', 'Gbar',  str(ChDensCA1pyr_Df.Ca_R_Gbar[11])],
        # ['Ca_L_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.Ca_L_Gbar[11])],
        # ['Ca_N_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.Ca_N_Gbar[11])],
        # ['K_SK_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.K_SK_Gbar[11])],
        # ['K_BK_chan', 'dend_2#', 'Gbar', str(ChDensCA1pyr_Df.K_BK_Gbar[11])],
        # ['Ca_conc', 'dend_2#', 'Ca_Basal', str(ChDensCA1pyr_Df.Ca_Basal[11])],
        #
        # ['Na_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.Na_Gbar[18])],
        # ['K_DR_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.K_DR_Gbar[18])],
        # ['K_A_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.K_A_Gbar[18])],
        # ['K_M_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.K_M_Gbar[18])],
        # ['h_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.h_Gbar[18])],
        # ['Ca_T_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_T_Gbar[18])],
        # ['Ca_R_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_R_Gbar[18])],
        # ['Ca_L_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_L_Gbar[18])],
        # ['Ca_N_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.Ca_N_Gbar[18])],
        # ['K_SK_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.K_SK_Gbar[18])],
        # ['K_BK_chan', 'dend_e#', 'Gbar', str(ChDensCA1pyr_Df.K_BK_Gbar[18])],
        # ['Ca_conc', 'dend_e#', 'Ca_Basal', str(ChDensCA1pyr_Df.Ca_Basal[18])],
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

        # ['apical_1_8', '1', '.', 'Vm', 'apical_1_8 Membrane potential'],
        # ['apical_1_8', '1', '.', 'Im', 'apical_1_8 current'],
        # ['dend_2_6', '1', '.', 'Vm', 'dend_2_6 Membrane potential'],
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

try:
    for ind in ChDensCA1pyr_Df.index:
        moose.element('/model/elec/'+ind+'/Ca_conc').B = 1000/(2*F*depth*18*np.pi*ChDensCA1pyr_Df['len'][ind]*ChDensCA1pyr_Df['diam'][ind])
except:
    pass

moose.reinit()


moose.start( runtime )


rdes.display()
print(time.time() - start_time)
