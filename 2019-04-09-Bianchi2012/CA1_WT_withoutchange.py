#Author: Anal Kumar
#Modeling R.Iyer et. al., 2017 which shows an increase in firing rate and decrease in variance on blocking SK channel.
#In the paper, Na and KDR were modelled as markov processes
#exec(open('/home/analkumar2/Thesis/Thesis work/Bianchi2002_WT_CA1/CA1_WT_withoutchange.py').read())


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools

Em = -0.075
# sm_area = 1764e-12
sm_area = 927e-12
RM = 2
CM = 0.01
Na_SGbar = 350
KDR_SGbar = 150
KA_SGbar = 5
KMGbar = 10
hGbar =  0.18
CaTGbar = 0.5
CaRGbar = 1
CaLGbar = 5
KsAHPGbar = 150
KmAHPGbar = 2475

ChP = '/home/analkumar2/Thesis/Thesis work/Bianchi2002_WT_CA1/ChannelProtos_Bianchi2012'
sm_diam = 28e-6
sm_len = sm_area/(sm_diam*np.pi)


try:
    # [moose.delete(x) for x in ['/model', '/library']]
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
        [ChP+'.KsAHP_Chan()', 'KsAHP_chan'],
        [ChP+'.KmAHP_Chan()', 'KmAHP_chan'],
        # [ChP+'.KSK_Chan()', 'KSK_chan'],
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
        ['KsAHP_chan', 'soma', 'Gbar', str(KsAHPGbar)],
        ['KmAHP_chan', 'soma', 'Gbar', str(KmAHPGbar)],
        ['Ca_conc', 'soma', 'thick', '177.9e-6'],
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<8) * 0.0' ],
        ['soma', '1', '.', 'inject', '(t>=3.0 && t<=4.00) ? 30e-12 : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'KsAHP_chan', 'Gk', 'Soma KsAHP conductance'],
        # ['soma', '1', 'KmAHP_chan', 'Gk', 'Soma KmAHP conductance'],
        # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
        # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
        # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

        # ['soma', '1', 'Na_Schan', 'Ik', 'Soma Sodium current'],
        # ['soma', '1', 'KDR_Schan', 'Ik', 'Soma Kdr current'],
        # ['soma', '1', 'KsAHP_chan', 'Ik', 'Soma KsAHP current'],
        # ['soma', '1', 'KmAHP_chan', 'Ik', 'Soma KmAHP current'],
        # ['soma', '1', 'KSK_chan', 'Ik', 'Soma KSK current'],
        # ['soma', '1', 'CaT_chan', 'Ik', 'Soma CaT current'],
        # ['soma', '1', 'CaR_Schan', 'Ik', 'Soma CaR_S current'],
        # ['soma', '1', 'CaL_Schan', 'Ik', 'Soma CaL_S current'],
    ],
)

rdes.buildModel()
# moose.element( '/model/elec/soma/vclamp' ).gain *= 0.001
moose.element('/model/elec/soma/Ca_conc').B = 28789637.7
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
