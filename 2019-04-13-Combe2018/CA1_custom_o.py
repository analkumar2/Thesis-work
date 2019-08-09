#Author: Anal Kumar
#pwd to the Thesis work folder, then run this. This way, it can run on both WSL as well as Linux.
#Modeling Combe et. al., 2018 which presumably uses an SK channel and BK channel in a 144 compartment morphology. Here only soma is used.
#exec(open('Combe2018/CA1_custom.py').read())

#importing neccesarry packages
import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools

#Fixed model parameters
celsius = 34
Em = -0.07
ENa = 0.050
EK = -0.080
Eh = -0.010
ECa = 0.140
F = 96485.3329

#Free model parameters
sm_area = 1477.4e-12 #np.pi*(r1+r2)*np.sqrt(h**2+(r1-r2)**2)
sm_vol = 467.5e-18 #np.pi/3*(r1**2+r2**2+r1*r2)
RM = 2
CM = 0.015
Na_SGbar = 1.2*0.035e4
KDR_SGbar = 2.2*0.015e4
hGbar = 1.8e-2
KA_SGbar = 7*0.0005e4
KMGbar = 0*0.001e4
CaLGbar = 0.1*0.00006e4
CaTGbar = 0.0003e4
CaRGbar = 0.00008e4
KSKGbar = 0.5*0.7*4.5*0.0001e4 #kca is mAHP. Corresponds to SK
KBKGbar = 5.5*0.9*1.5*0.03e4 #mykca is fAHP. Corresponds to BK
CaConc_tau = 20e-3
CaConc_B = 10000/2/F/0.05/18/sm_area #Correct the units. Don't know how to correct

#Derived model parameters
sm_diam = 4*sm_vol/sm_area
sm_len = sm_area**2/4/sm_vol

#ChannelProtos file
ChP = 'Combe2018/ChannelProtos_Combe2018'

#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

rdes = rd.rdesigneur(
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
        ['KBK_chan', 'soma', 'Gbar', str(KBKGbar)],
        ['KSK_chan', 'soma', 'Gbar', str(KSKGbar)],
        ['Ca_conc', 'soma', 'thick', '177.9e-6'],
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<8) * 0.0' ],
        ['soma', '1', '.', 'inject', '(t>=3.0 && t<=3.50) ? 300e-12 : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'KsAHP_chan', 'Gk', 'Soma KsAHP conductance'],
        ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
        # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
        # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

        # ['soma', '1', 'Na_Schan', 'Ik', 'Soma Sodium current'],
        # ['soma', '1', 'KDR_Schan', 'Ik', 'Soma Kdr current'],
        # ['soma', '1', 'KsAHP_chan', 'Ik', 'Soma KsAHP current'],
        ['soma', '1', 'KSK_chan', 'Ik', 'Soma KSK current'],
        # ['soma', '1', 'CaT_chan', 'Ik', 'Soma CaT current'],
        # ['soma', '1', 'CaR_Schan', 'Ik', 'Soma CaR_S current'],
        # ['soma', '1', 'CaL_Schan', 'Ik', 'Soma CaL_S current'],
    ],
)

rdes.buildModel()
# moose.element( '/model/elec/soma/vclamp' ).gain *= 0.001
moose.element('/model/elec/soma/Ca_conc').B = CaConc_B
moose.element('/model/elec/soma/Ca_conc').tau = CaConc_tau
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
