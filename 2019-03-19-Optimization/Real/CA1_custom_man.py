#Author: Anal Kumar
#pwd to the Thesis work folder, then run this. This way, it can run on both WSL as well as Linux.
#Modeling Combe et. al., 2018 which presumably uses an SK channel and BK channel in a 144 compartment morphology. Here only soma is used.
#exec(open('Optimization/Custom/Real/CA1_custom_man.py').read())

#importing neccesarry packages
import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import csv

#Fixed model parameters
celsius = 34
ENa = 0.050
EK = -0.080
EKDR = -0.077
Eh = -0.010
ECa = 0.140
F = 96485.3329
R_input = 152e6
E_rest = -0.063
Cm = 134e-12
CaConc_tau = 20e-3
a1 = 4.25e-6
a2 = 1.8e-4
a3 = 1.07e-1
a4 = 6.52e-4
a5 = 5.22e-3
a6 = 1.48e-4
a7 = 9.75e-6
a8 = 1.42e-2
a9 = 1.60e-5
a10 = 6.49e-6

#Free model parameters
sm_area = 1477.4e-12
sm_vol = 467.5e-18
Na_SGbar = 645
KDR_SGbar = 476
hGbar = 0.57
KA_SGbar = 12.8
KMGbar = 19.8
CaLGbar = 16.7
CaTGbar = 0.07
CaRGbar = 1.41
KSKGbar = 0.01
KBKGbar = 3996
CaConc_B = 100000/2/F/0.05/18/sm_area #Correct the units. Don't know how to correct

#Derived model parameters
sm_diam = 4*sm_vol/sm_area
sm_len = sm_area**2/4/sm_vol/np.pi
RM = np.maximum(1e-12,1/(1/(R_input*sm_area) - Na_SGbar*a1 - KDR_SGbar*a2 - hGbar*a3 - KA_SGbar*a4 - KMGbar*a5 - CaLGbar*a6 - CaTGbar*a7 - CaRGbar*a8 - KSKGbar*a9 - KBKGbar*a10))
Em = (E_rest/(R_input*sm_area) - Na_SGbar*ENa*a1 - KDR_SGbar*EKDR*a2 - hGbar*Eh*a3 - KA_SGbar*EK*a4 - KMGbar*EK*a5 - CaLGbar*ECa*a6 - CaTGbar*ECa*a7 - CaRGbar*ECa*a8 - KSKGbar*EK*a9 - KBKGbar*EK*a10)/(1/(R_input*sm_area) - Na_SGbar*a1 - KDR_SGbar*a2 - hGbar*a3 - KA_SGbar*a4 - KMGbar*a5 - CaLGbar*a6 - CaTGbar*a7 - CaRGbar*a8 - KSKGbar*a9 - KBKGbar*a10)

#ChannelProtos file
ChP = 'Optimization/Custom/Real/ChannelProtos_Combe2018'

#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

inp_curr = 300e-12
rdes = rd.rdesigneur(
    elecDt = 10e-6,
    elecPlotDt = 50e-6,
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
        ['soma', 'RM', str(RM), 'RA', '1.5', 'Cm', str(Cm), 'initVm', str(E_rest), 'Em', str(Em)],
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
        # ['soma', '1', '.', 'vclamp', f'-0.063 + (t>3 && t<5) * 0' ],
        ['soma', '1', '.', 'inject', f'(t>=3.0813 && t<=3.5813) ? {inp_curr} : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'KM_chan', 'Gk', 'Soma KM conductance'],
        # ['soma', '1', 'KsAHP_chan', 'Gk', 'Soma KsAHP conductance'],
        # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
        # ['soma', '1', 'KBK_chan', 'Gk', 'Soma KBK conductance'],
        # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
        # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
        # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

        # ['soma', '1', 'Na_Schan', 'Ik', 'Soma Sodium current'],
        # ['soma', '1', 'KDR_Schan', 'Ik', 'Soma Kdr current'],
        # ['soma', '1', 'KsAHP_chan', 'Ik', 'Soma KsAHP current'],
        # ['soma', '1', 'KSK_chan', 'Ik', 'Soma KSK current'],
        # ['soma', '1', 'CaT_chan', 'Ik', 'Soma CaT current'],
        # ['soma', '1', 'CaR_Schan', 'Ik', 'Soma CaR_S current'],
        # ['soma', '1', 'CaL_Schan', 'Ik', 'Soma CaL_S current'],
    ],
)

rdes.buildModel()
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
