#Author: Anal Kumar
#exec(open('Optimization/Custom/Real_features/CA1_WT_Sri2015.py').read())
#In this file I play with the densities to get the required feature value


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
F = 96485.3329
ChP = 'Optimization/Custom/Real_features/ChannelProtos_Sri2015_base'

sm_diam = 100e-6
sm_len = 100e-6
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
RA = 1.0

Em = -0.065
RM = 3.23
CM = 0.0083
Ca_tau = 0.029
Ca_B = 1/(sm_vol*F*2)
Na_Gbar = 1
K_DR_Gbar = 1
K_A_Gbar = 1
K_M_Gbar = 1
h_Gbar =  1
Ca_T_Gbar = 1
Ca_R_Gbar = 1
Ca_L_Gbar = 1
Ca_N_Gbar = 1
K_SK_Gbar = 1
K_BK_Gbar = 1

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    for ff in moose.le('library'):
        if ff != '/library[0]/KBK_chan':
            moose.delete(ff)
    [moose.delete(x) for x in ['/model']]
except:
    pass

rdes = rd.rdesigneur(
    # cellProto = [['somaProto', 'soma', 12.76e-6, 0.01e-6]],
    cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
    chanProto = [
        [ChP+'.Na_Chan()', 'Na_chan'],
        [ChP+'.KDR_Chan()', 'KDR_chan'],
        [ChP+'.KA_Chan()', 'KA_chan'],
        [ChP+'.KM_Chan()', 'KM_chan'],
        [ChP+'.h_Chan()', 'h_chan'],
        [ChP+'.CaT_Chan()', 'CaT_chan'],
        [ChP+'.CaR_Chan()', 'CaR_chan'],
        [ChP+'.CaL_Chan()', 'CaL_chan'],
        [ChP+'.CaN_Chan()', 'CaN_chan'],
        [ChP+'.KBK_Chan()', 'KBK_chan'],
        [ChP+'.KSK_Chan()', 'KSK_chan'],
        [ChP+'.Ca_Conc()', 'Ca_conc'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', '1.5', 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
    ],
    chanDistrib = [
        ['Na_chan', 'soma', 'Gbar', str(Na_Gbar)],
        ['KDR_chan', 'soma', 'Gbar', str(K_DR_Gbar)],
        ['KA_chan', 'soma', 'Gbar', str(K_A_Gbar)],
        ['KM_chan', 'soma', 'Gbar', str(K_M_Gbar)],
        ['h_chan', 'soma', 'Gbar', str(h_Gbar)],
        ['CaT_chan', 'soma', 'Gbar', str(Ca_T_Gbar)],
        ['CaR_chan', 'soma', 'Gbar', str(Ca_R_Gbar)],
        ['CaL_chan', 'soma', 'Gbar', str(Ca_L_Gbar)],
        ['CaN_chan', 'soma', 'Gbar', str(Ca_N_Gbar)],
        ['KSK_chan', 'soma', 'Gbar', str(K_SK_Gbar)],
        ['KBK_chan', 'soma', 'Gbar', str(K_BK_Gbar)],
        ['Ca_conc', 'soma', 'thick', '177.9e-6'],
    ],
    stimList = [
        ['soma', '1', '.', 'vclamp', '-0.075 + (t>3 && t<6.5) * 0.010' ],
        # ['soma', '1', '.', 'inject', '(t>=2.0 && t<=2.5) ? 150e-12 : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
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
    ],
)

rdes.buildModel()
try:
    moose.element( '/model/elec/soma/vclamp' ).gain *= 0.1
except:
    pass
moose.element('/model/elec/soma/Ca_conc').B = Ca_B
moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
moose.reinit()

# data = moose.Neutral('/data')
# somaNa_SXgate = moose.Table('/data/somaNa_SXgate')
# somaNa_S = moose.element('/model/elec/soma/Na_Schan')
# moose.connect(somaNa_SXgate, 'requestOut', somaNa_S, 'getX')
# somaNa_SYgate = moose.Table('/data/somaNa_SYgate')
# moose.connect(somaNa_SYgate, 'requestOut', somaNa_S, 'getY')
# somaNa_SZgate = moose.Table('/data/somaNa_SZgate')
# moose.connect(somaNa_SZgate, 'requestOut', somaNa_S, 'getZ')

moose.start( 3 )

# plt.figure(100)
# plt.plot(np.linspace(0,6,60000), somaNa_SXgate.vector**3, label='X gate')
# plt.plot(np.linspace(0,6,60000), somaNa_SYgate.vector, label='Y gate')
# plt.plot(np.linspace(0,6,60000), somaNa_SZgate.vector, label='Z gate')
# plt.xlabel('Time (s)')
# plt.ylabel('Gating probabilities')
# plt.title('Na_S gates ModelDB')
# plt.legend()

rdes.display()
