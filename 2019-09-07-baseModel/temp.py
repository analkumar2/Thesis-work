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
Pr = {}

try:
    [moose.delete(x) for x in ['/model', '/library']]
    # moose.delete('/model')
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
        ['soma', '1', '.', 'inject', f'(t>=1 && t<=1.5) ? 0e-12 : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
        # ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
        # ['soma', '1', '.', 'inject', 'Injected current MOOSE'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
        ['soma', '1', 'K_DR_chan', 'Ik', 'K_DR_chan current MOOSE'],
        ['soma', '1', 'K_DR_chan', 'Gk', 'K_DR_chan conductance MOOSE'],
    ],
)

# def set_initialgatevalues(Vrest=-0.060, Carest=1e-3):
#     for chan in moose.wildcardFind('/library/#[TYPE=HHChannel]'):
#         if chan.Xpower>0:
#             gate = moose.element( chan.path + '/gateX' )
#             v = np.linspace(gate.min,gate.max,gate.divs)
#             idx = np.argmin(abs(v-Vrest))
#             chan.X = gate.tableA[idx]/gate.tableB[idx]
#         if chan.Ypower>0:
#             gate = moose.element( chan.path + '/gateY' )
#             v = np.linspace(gate.min,gate.max,gate.divs)
#             idx = np.argmin(abs(v-Vrest))
#             chan.Y = gate.tableA[idx]/gate.tableB[idx]
#         if chan.Zpower>0:
#             gate = moose.element( chan.path + '/gateZ' )
#             vorca = np.linspace(gate.min,gate.max,gate.divs)
#             if chan.useConcentration == 1:
#                 idx = np.argmin(abs(vorca-Carest))
#             else:
#                 idx = np.argmin(abs(vorca-Vrest))
#             chan.Z = gate.tableA[idx]/gate.tableB[idx]
# set_initialgatevalues(Vrest=-0.0716, Carest=1e-3)

def calcPr(Vrest=-0.07, Carest=0.13e-3):
    for chan in moose.wildcardFind('/library/#[TYPE=HHChannel]'):
        if chan.Xpower>0:
            gate = moose.element( chan.path + '/gateX' )
            v = np.linspace(gate.min,gate.max,gate.divs)
            idx = np.argmin(abs(v-Vrest))
            Pr[chan.name] = gate.tableA[idx]/gate.tableB[idx] #Since we don't care about channels other than h, KM, and KD
calcPr(Vrest=-0.07, Carest=0.13e-3)

def set_initialGk(Vrest=-0.07, Carest=0.13e-3):
    for chan in moose.wildcardFind('/model/elec/soma/#[TYPE=ZombieHHChannel]'):
        chan.Gk = chan.Gbar*chan.X**chan.Xpower*chan.Y**chan.Ypower*chan.Z**chan.Zpower
    for chan in moose.wildcardFind('/model/elec/soma/#[TYPE=HHChannel2D]'):
        chan.Gk = chan.Gbar*chan.X**chan.Xpower*chan.Y**chan.Ypower*chan.Z**chan.Zpower

rdes.buildModel()
K_DR_chan = moose.element('/model/elec/soma/K_DR_chan')
plot_K_DR_chanX = moose.Table('/model/graphs/plot_K_DR_chanX')
moose.connect(plot_K_DR_chanX, 'requestOut', K_DR_chan, 'getX')
clk = moose.element('/clock')
plot_t = moose.Table('/model/graphs/plott')
moose.connect(plot_t, 'requestOut', clk, 'getCurrentTime')

try:
    moose.element( '/model/elec/soma/vclamp' ).gain = 0.01*sm_area/0.00005
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*0.00005
    moose.element( '/model/elec/soma/vclamp' ).ti = 0.00005
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass

moose.reinit()
set_initialGk(Vrest=-0.07, Carest=0.13e-3)
print(K_DR_chan.Gk)
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
