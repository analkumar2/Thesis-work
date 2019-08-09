#Author: Anal Kumar
#exec(open('Srikanth2015/CA1_WT_Sri2015.py').read())


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
F = 96485.3329

Em = -0.065
RM = 3.5
CM = 0.01
Na_SGbar = 70
KDR_SGbar = 30
KA_SGbar = 80
KMGbar = 0.01
hGbar =  0.8
CaTGbar = 1
CaRGbar = 1
CaLGbar = 1
CaNGbar = 1
KSKGbar = 0.01
KBKGbar = 0.01

ChP = 'Srikanth2015/ChannelProtos_Sri2015_base'
sm_diam = 100e-6
sm_len = 100e-6

Ca_tau = 0.2/7 #Not sure about this
Ca_B = np.pi*sm_diam*1e10/2/F #Not sure about this

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
        [ChP+'.Na_SChan()', 'Na_Schan'],
        [ChP+'.KDR_SChan()', 'KDR_Schan'],
        [ChP+'.KA_SChan()', 'KA_Schan'],
        [ChP+'.KM_Chan()', 'KM_chan'],
        [ChP+'.h_Chan()', 'h_chan'],
        [ChP+'.CaT_Chan()', 'CaT_chan'],
        [ChP+'.CaR_SChan()', 'CaR_Schan'],
        [ChP+'.CaL_SChan()', 'CaL_Schan'],
        [ChP+'.CaN_SChan()', 'CaN_Schan'],
        [ChP+'.KBK_Chan()', 'KBK_chan'],
        [ChP+'.KSK_Chan()', 'KSK_chan'],
        [ChP+'.Ca_Conc()', 'Ca_conc'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', '1.5', 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
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
        ['CaN_Schan', 'soma', 'Gbar', str(CaNGbar)],
        ['KSK_chan', 'soma', 'Gbar', str(KSKGbar)],
        ['KBK_chan', 'soma', 'Gbar', str(KBKGbar)],
        ['Ca_conc', 'soma', 'thick', '177.9e-6'],
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<6.5) * 0.050' ],
        ['soma', '1', '.', 'inject', '(t>=3.0 && t<=3.5) ? 150e-12 : 0'],
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

moose.start( 6 )

# plt.figure(100)
# plt.plot(np.linspace(0,6,60000), somaNa_SXgate.vector**3, label='X gate')
# plt.plot(np.linspace(0,6,60000), somaNa_SYgate.vector, label='Y gate')
# plt.plot(np.linspace(0,6,60000), somaNa_SZgate.vector, label='Z gate')
# plt.xlabel('Time (s)')
# plt.ylabel('Gating probabilities')
# plt.title('Na_S gates ModelDB')
# plt.legend()

rdes.display()
