#exec(open('temp.py').read())

import os
import sys
import csv
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import numpy as np
import warnings
import moose
import pickle
import rdesigneur as rd

sm_diam=60e-6; sm_len=60e-6
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    moose.delete('/model')
    # moose.delete('/library/Ca_conc')
except:
    pass

rdes = rd.rdesigneur(
    elecPlotDt = 0.00005,
    cellProto = [
        ['somaProto', 'soma', sm_len, sm_diam],
    ],
    chanProto = [
        ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
        ['Na_Chan_(Migliore2018).Na_Chan()', 'Na_chan'],
    ],
    passiveDistrib = [
        ['soma', 'RM', '4', 'CM', '0.01', 'initVm', '-0.070', 'Em', '-0.070'],
    ],
    chanDistrib = [
        ['K_DR_chan', 'soma', 'Gbar', '30'],
        ['Na_chan', 'soma', 'Gbar', '100'],
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.070 + (t>1 && t<1.5) * 0.070' ],
        ['soma', '1', '.', 'inject', f'(t>=1 && t<=1.5) ? 150e-12 : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
        ['soma', '1', '.', 'inject', 'Injected current MOOSE'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
        # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
        # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
    ],
)

rdes.buildModel()
clk = moose.element('/clock')
plott = moose.Table('/model/graphs/plott')
moose.connect(plott, 'requestOut', clk, 'getCurrentTime')

try:
    moose.element( '/model/elec/soma/vclamp' ).gain = 0.01*sm_area/0.00005
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*0.00005
    moose.element( '/model/elec/soma/vclamp' ).ti = 0.00005
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass

moose.reinit()
moose.start(2)
rdes.display()

# # moose.delete('/model/elec/soma/vclamp')
# vclamp = moose.VClamp('/model/elec/soma/vclamp')
# vclamp.gain = 0.01*sm_area/0.00005
# vclamp.tau = 5*0.00005
# vclamp.ti = 0.00005
# vclamp.td = 0
# vclamp.command = -0.04
# print(vclamp.command)
# # moose.connect(vclamp, 'currentOut', '/model/elec/soma', 'injectMsg')
# # moose.connect('/model/elec/soma', 'VmOut', vclamp, 'sensedIn')
# # command = moose.Function('/model/stims/stim0')
# # moose.element('/model/stims/stim0').expr = '-0.070 + (t>1 && t<1.5) * 0.04'
#
# comp = moose.element('/model/elec/soma')
# # # command = moose.PulseGen('/model/elec/soma/command')
# # # command.delay[0] = 1000.0
# # # command.width[0] = 1500.0
# # # command.level[0] = -10000
# # # command.delay[1] = 1e9
# # # moose.connect(command, 'output', vclamp, 'commandIn')
# moose.connect(vclamp, 'currentOut', comp, 'injectMsg')
# moose.connect(comp, 'VmOut', vclamp, 'sensedIn')
# #
# moose.reinit()
# print(vclamp.command)
# moose.start(2)
# rdes.display()
