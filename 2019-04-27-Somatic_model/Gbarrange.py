#Author: Anal Kumar
#exec(open('Somatic model/Gbarrange.py').read())
#In this file I play with the densities to get the required feature value


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
F = 96485.3329
ChP = 'Somatic model/ChannelProtos_Sri2015_base'

# sm_diam = 100e-6
# sm_len = 100e-6
# sm_vol = np.pi/4*sm_diam**2*sm_len
# sm_area = np.pi*sm_diam*sm_len
# RA = 1.0

# Em = -0.070
# RM = 3.23
# CM = 0.0083
# Ca_tau = 0.029
# Ca_B = 1/(sm_vol*F*2)
# Na_Gbar = 76
# K_DR_Gbar = 43.4
# K_A_Gbar = 20.7
# K_M_Gbar = 3.88294e-02
# h_Gbar =  0.11
# Ca_T_Gbar = 1.57
# Ca_R_Gbar = 0.07
# Ca_L_Gbar = 0.51
# Ca_N_Gbar = 1.76
# K_SK_Gbar = 4.28006e-02*2
# K_BK_Gbar = 5.08559e-03

with open('Somatic model/best_params.pkl', 'rb') as f:
    Em, RM, RA, CM, sm_diam, sm_len, \
    Ca_Basal, Ca_tau, Ca_B, Na_Gbar, K_DR_Gbar, K_A_Gbar, \
    K_M_Gbar, h_Gbar, Ca_T_Gbar, Ca_R_Gbar, \
    Ca_L_Gbar, Ca_N_Gbar, K_SK_Gbar, K_BK_Gbar = pickle.load(f)

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    for ff in moose.le('library'):
        if ff != '/library[0]/K_BK_chan':
            moose.delete(ff)
    [moose.delete(x) for x in ['/model']]
except:
    pass

def makeModel(stim_start = 2.0, stim_end = 2.5, curr = 150e-12):
    rdes = rd.rdesigneur(
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
            [ChP+'.Ca_N_Chan()', 'Ca_N_chan'],
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
            ['soma', '1', '.', 'vclamp', f'-0.075 + (t>3 && t<3.20) * 0.070' ],
            # ['soma', '1', '.', 'inject', f'(t>={stim_start} && t<={stim_end}) ? {curr} : 0'],
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
            ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
            ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
            # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
            # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
            # ['soma', '1', 'K_SK_chan', 'Gk', 'Soma KSK conductance'],
            # ['soma', '1', 'K_BK_chan', 'Gk', 'Soma KBK conductance'],
            # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
            # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
            # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],
            ['soma', '1', 'Ca_N_chan', 'Gk', 'Soma Ca_N conductance'],

            # ['soma', '1', 'Na_Schan', 'Ik', 'Soma Sodium current'],
            # ['soma', '1', 'KDR_Schan', 'Ik', 'Soma Kdr current'],
            # ['soma', '1', 'KA_Schan', 'Ik', 'Soma KA current'],
            # ['soma', '1', 'KM_chan', 'Ik', 'Soma KM current'],
            # ['soma', '1', 'h_chan', 'Ik', 'Soma h current'],
            # ['soma', '1', 'K_SK_chan', 'Ik', 'Soma K_SK current'],
            # ['soma', '1', 'K_BK_chan', 'Ik', 'Soma KBK current'],
            # ['soma', '1', 'CaT_chan', 'Ik', 'Soma CaT current'],
            # ['soma', '1', 'CaR_Schan', 'Ik', 'Soma CaR_S current'],
            # ['soma', '1', 'CaL_Schan', 'Ik', 'Soma CaL_S current'],
            ['soma', '1', 'Ca_N_chan', 'Ik', 'Soma Ca_N current'],
        ],
    )
    return rdes

rdes = makeModel()
rdes.buildModel()
try:
    moose.element( '/model/elec/soma/vclamp' ).gain *= 0.05
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

moose.start( 5 )

# plt.figure(100)
# plt.plot(np.linspace(0,6,60000), somaNa_SXgate.vector**3, label='X gate')
# plt.plot(np.linspace(0,6,60000), somaNa_SYgate.vector, label='Y gate')
# plt.plot(np.linspace(0,6,60000), somaNa_SZgate.vector, label='Z gate')
# plt.xlabel('Time (s)')
# plt.ylabel('Gating probabilities')
# plt.title('Na_S gates ModelDB')
# plt.legend()

rdes.display()
