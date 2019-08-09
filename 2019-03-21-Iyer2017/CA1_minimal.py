#Replicated Iyer2017 dopamine model
#exec(open('/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/Real/CA1_minimal.py').read())

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
F = 96485.3329

# R_input = 107e12
# Erest = -0.065
#
Em = -0.045
RM = 3.33
#
sm_area = 4*np.pi*0.5e-6**2
sm_vol = 4/3*np.pi*0.5e-6**3
Na_SGbar = 36
KDR_SGbar = 4
KA_SGbar = 40
CaL_SGbar = 50
KSK_Gbar = 50
CaB = 1/(2*F*sm_vol)
Catau = 1.4
curr_inj = 1e-12
#
ChP = '/home/analkumar2/Thesis/Thesis work/Optimizatiom/Custom/Real/ChannelProtos_Iyer2017'
CM = 0.01
sm_diam = 4*sm_vol/sm_area
sm_len = sm_area**2/4/sm_vol

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
        [ChP+'.CaL_SChan()', 'CaL_Schan'],
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
        ['CaL_Schan', 'soma', 'Gbar', str(CaL_SGbar)],
        ['KSK_chan', 'soma', 'Gbar', str(KSK_Gbar)],
        ['Ca_conc', 'soma', 'thick', '177.9e-6'],
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.065' ],
        ['soma', '1', '.', 'inject', f'(t>=3.0814 && t<=3.5814) ? {curr_inj} : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', './Ca_conc', 'Ca', 'Calcium concentration'],
    ],
)

rdes.buildModel()
# moose.element( '/model/elec/soma/vclamp' ).gain *= 0.001
moose.element('/model/elec/soma/Ca_conc').B = CaB
moose.element('/model/elec/soma/Ca_conc').tau = Catau
moose.reinit()

moose.start( 6 )

rdes.display()
