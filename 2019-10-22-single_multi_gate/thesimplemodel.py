# exec(open('thesimplemodel.py').read())

import os
import sys
from neo.io import AxonIO
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import numpy as np
import warnings
import moose
import pickle
import rdesigneur as rd
import time

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    moose.delete('/model')
    moose.delete('/library')
except:
    pass

rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    cellProto = [
        ['somaProto', 'soma', sm_diam, sm_len],
    ],
    chanProto = [
        ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
        ['Na_Chan_(Migliore2018).Na_Chan()', 'Na_chan'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(1/gl), 'CM', str(CM), 'initVm', str(Vrest), 'Em', str(Vrest)],
    ],
    chanDistrib = [
        ['K_DR_chan', 'soma', 'Gbar', K_DR_changbar],
        ['Na_chan', 'soma', 'Gbar', Na_changbar],
    ],
    stimList = [
        ['soma', '1', '.', 'vclamp', str(Vrest) ],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
        # ['soma', '1', ',', 'inject', 'Injected current MOOSE'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
        # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
        # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
    ],
)

rdes.buildModel()

#Changing vclamp parameters
try:
    moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass

#Setting Ca_conc B value
try:
    moose.element('/model/elec/soma/Ca_conc').B = 1000e3/(2*F*depth*np.pi*sm_diam*sm_len*2)
    # moose.element('/model/elec/soma/Ca_conc').B *= 2
    # moose.element('/model/elec/soma/Ca_conc').B = 0
except:
    pass

moose.reinit()
moose.start(0.4)
# rdes.display()
