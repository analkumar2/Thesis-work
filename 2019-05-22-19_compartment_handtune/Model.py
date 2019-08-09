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
        ChanDistrib_list.append([coul[:-4]+'chan', ind, 'Gbar', str(ChDensCA1pyr_Df[coul][ind])])
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

    passiveDistrib = PassiveDistrib_list,

    chanDistrib = ChanDistrib_list,

    # chanDistrib = ['Na_Gbar', 'soma', 'Gbar', '300'],

    stimList = [
        ['soma', '1', '.', 'vclamp', f'-0.070 + (t>{preStimTime} && t<{preStimTime+injectTime}) * 0.005' ],
        # ['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {Injectcurr} : 0'],
    ],

    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', '.', 'Im', 'Soma current'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_chan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'K_DR_chan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'K_A_chan', 'Gk', 'Soma KA conductance'],
        # ['soma', '1', 'K_M_chan', 'Gk', 'Soma KM conductance'],
        # ['soma', '1', 'h_chan', 'Gk', 'Soma h conductance'],
        # ['soma', '1', 'K_SK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'K_BK_chan', 'Gk', 'Soma KBK conductance'],
        # ['soma', '1', 'Ca_T_chan', 'Gk', 'Soma CaT conductance'],
        # ['soma', '1', 'Ca_R_chan', 'Gk', 'Soma CaR_S conductance'],
        # ['soma', '1', 'Ca_L_chan', 'Gk', 'Soma CaL_S conductance'],
        # ['soma', '1', 'Ca_N_chan', 'Gk', 'Soma CaN_S conductance'],
        #
        # ['soma', '1', 'Na_chan', 'Ik', 'Soma Sodium current'],
        # ['soma', '1', 'K_DR_chan', 'Ik', 'Soma Kdr current'],
        # ['soma', '1', 'K_A_chan', 'Ik', 'Soma KA current'],
        # ['soma', '1', 'K_M_chan', 'Ik', 'Soma KM current'],
        # ['soma', '1', 'h_chan', 'Ik', 'Soma h current'],
        # ['soma', '1', 'K_SK_chan', 'Ik', 'Soma KSK current'],
        # ['soma', '1', 'K_BK_chan', 'Ik', 'Soma KBK current'],
        # ['soma', '1', 'Ca_T_chan', 'Ik', 'Soma CaT current'],
        # ['soma', '1', 'Ca_R_chan', 'Ik', 'Soma CaR_S current'],
        # ['soma', '1', 'Ca_L_chan', 'Ik', 'Soma CaL_S current'],
        # ['soma', '1', 'Ca_N_chan', 'Ik', 'Soma CaN_S current'],

        # ['apical_1_8', '1', '.', 'Vm', 'apical_1_8 Membrane potential'],
        # ['apical_1_8', '1', '.', 'Im', 'apical_1_8 current'],
        # ['dend_2_6', '1', '.', 'Vm', 'dend_2_6 Membrane potential'],
        # ['dend_2_6', '1', '.', 'Im', 'dend_2_6 current'],
    ],
)

rdes.buildModel()

try:
    moose.element( '/model/elec/soma/vclamp' ).gain = ChDensCA1pyr_Df.CM[0]*np.pi*ChDensCA1pyr_Df.len[0]*ChDensCA1pyr_Df.diam[0]/elecPlotDt
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
