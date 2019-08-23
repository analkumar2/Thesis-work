# exec(open('somacat_MOOSE_ghk.py').read())

import moose
import numpy as np
import matplotlib.pyplot as plt
import rdesigneur as rd

try:
    [moose.delete(x) for x in ['/model', '/library']]
    # [moose.delete(x) for x in ['/model']]
except:
    pass

# ChP = ''

F = 96485.3329
sm_diam = 20e-6
sm_len = 20e-6
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
RA = 35.4
RM = 0.1
CM = 0.01
Em = -0.070
depth = 0.1 # No units. For ca_conc B
elecPlotDt = 0.00005
preStimTime = 1
injectTime = 0.5
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 0.5e-9

rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    cellProto = [
        ['somaProto', 'soma', sm_diam, sm_len],
    ],
    chanProto = [
        ['Ca_T_Chan_(Migliore2018).Ca_T_Chan()', 'Ca_T_chan'],
        ['Ca_Conc_(Common).Ca_Conc()', 'Ca_conc'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
    ],
    chanDistrib = [
        ['Ca_T_chan', 'soma', 'Gbar', str(30)],
        ['Ca_conc', 'soma', 'Ca_Basal', str(0.05e-3)],
    ],
    stimList = [
        ['soma', '1', '.', 'vclamp', f'-0.065 + (t>{preStimTime} && t<{preStimTime+injectTime}) * 0.065' ],
        # ['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {Injectcurr} : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
        # ['soma', '1', ',', 'inject', 'Injected current'],
        ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc'],
        ['soma', '1', 'Ca_T_chan', 'Ik', 'Ca_T channel current'],
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
    moose.element('/model/elec/soma/Ca_conc').B = 1000e3/(2*F*depth*np.pi*sm_diam*sm_len*2)
    # moose.element('/model/elec/soma/Ca_conc').B *= 2
    # moose.element('/model/elec/soma/Ca_conc').B = 0
except:
    pass

moose.reinit()


moose.start( 2 )


rdes.display()