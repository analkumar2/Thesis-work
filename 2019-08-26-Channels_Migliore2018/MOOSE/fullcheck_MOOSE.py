# exec(open('fullcheck_MOOSE.py').read())

import moose
import numpy as np
import matplotlib.pyplot as plt
import rdesigneur as rd

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
    # moose.delete('/library/K_M_chan')
except:
    pass

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
Injectcurr = 7e-9

rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    cellProto = [
        ['somaProto', 'soma', sm_diam, sm_len],
    ],
    chanProto = [
        ['Ca_Conc_(Common).Ca_Conc()', 'Ca_conc'],
        ['Ca_L_Chan_(Migliore2018).Ca_L_Chan()', 'Ca_L_chan'],
        ['Ca_N_Chan_(Migliore2018).Ca_N_Chan()', 'Ca_N_chan'],
        ['Ca_T_Chan_(Migliore2018).Ca_T_Chan()', 'Ca_T_chan'],
        ['h_Chan_(Migliore2018).h_Chan()', 'h_chan'],
        ['K_A_Chan_(Migliore2018).K_A_Chan()', 'K_A_chan'],
        ['K_BK_Chan_(Migliore2018).K_BK_Chan()', 'K_BK_chan'],
        ['K_D_Chan_(Migliore2018).K_D_Chan()', 'K_D_chan'],
        ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
        ['K_M_Chan_(Migliore2018).K_M_Chan()', 'K_M_chan'],
        ['K_SK_Chan_(Migliore2018).K_SK_Chan()', 'K_SK_chan'],
        ['Na_Chan_(Migliore2018).Na_Chan()', 'Na_chan'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
    ],
    chanDistrib = [
        ['Ca_conc', 'soma', 'CaBasal', str(0.05e-3)],
        ['Ca_L_chan', 'soma', 'Gbar', '0.003e4'],
        ['Ca_N_chan', 'soma', 'Gbar', '0.0003e4'],
        ['Ca_T_chan', 'soma', 'Gbar', '0.003e4'],
        ['h_chan', 'soma', 'Gbar', '0.0001e4'],
        ['K_A_chan', 'soma', 'Gbar', '0.008e4'],
        ['K_BK_chan', 'soma', 'Gbar', '0.01e4'],
        ['K_D_chan', 'soma', 'Gbar', '0.01e4'],
        ['K_DR_chan', 'soma', 'Gbar', '0.003e4'],
        ['K_M_chan', 'soma', 'Gbar', '0.0001e4'],
        ['K_SK_chan', 'soma', 'Gbar', '0.01e4'],
        ['Na_chan', 'soma', 'Gbar', '0.010e4'],
    ],
    stimList = [
        ['soma', '1', '.', 'vclamp', f'-0.065 + (t>{preStimTime} && t<{preStimTime+injectTime}) * 0.045' ],
        # ['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {Injectcurr} : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
        # ['soma', '1', ',', 'inject', 'Injected current MOOSE'],
        ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
        # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
        # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
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
