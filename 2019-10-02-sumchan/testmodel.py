#exec(open('testmodel.py').read())

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

# vfinal = -0.030

# Define constants not to be changed
F = 96485.3329
elecPlotDt = 50e-6
elecDt = 10e-6
sm_len = 60e-6
sm_diam = 60e-6
sm_area = np.pi*sm_len*sm_diam
CM = 0.01

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    moose.delete('/model')
    moose.delete('/library')
except:
    pass

rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    elecDt = elecDt,
    cellProto = [
        ['somaProto', 'soma', sm_diam, sm_len],
    ],
    chanProto = [
        ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
        ['K_A_Chan_(Migliore2018).K_A_Chan()', 'K_A_chan'],
    ],
    passiveDistrib = [
        ['soma', 'RM', str(1), 'CM', str(CM), 'initVm', str(-0.07), 'Em', str(-0.07)],
    ],
    chanDistrib = [
        ['K_DR_chan', 'soma', 'Gbar', '50'],
        ['K_A_chan', 'soma', 'Gbar', '50'],
    ],
    stimList = [
         ['soma', '1', '.', 'vclamp', f'-0.100 + (t>0.5 && t<1) * ({vfinal}-(-0.100))' ]
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
        # ['soma', '1', 'K_DR_chan', 'Ik', 'K_DR Channel current'],
        # ['soma', '1', 'K_A_chan', 'Ik', 'K_A Channel current'],
        # ['soma', '1', ',', 'inject', 'Injected current MOOSE'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
        # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
        # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
    ],
)

rdes.buildModel()

clk = moose.element('/clock')
plott = moose.Table('/model/graphs/plott')
moose.connect(plott, 'requestOut', clk, 'getCurrentTime')

#Changing vclamp parameters
try:
    moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecDt
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecDt
    moose.element( '/model/elec/soma/vclamp' ).ti = elecDt
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass

moose.reinit()
moose.start(1.5)
# rdes.display()

tvec = plott.vector
I = moose.element('/model/graphs/plot1')
plt.plot(tvec, I.vector, label=f'{vfinal:.2f}')

# import numpy as np
# import matplotlib.pyplot as plt
# plt.cla()
# ax = plt.subplot(111)
# for vfinal in np.arange(-0.080,0.040,0.010):
#     exec(open('testmodel.py').read())
#
# plt.title('Voltage clamp just KA and KDR')
# plt.ylabel('Holding current (A)')
# plt.xlabel('time (s)')
# plt.legend()
# pickle.dump(ax, open('justKADR_Vclamp.pickle', 'wb'))
#
# import numpy as np
# import matplotlib.pyplot as plt
# import pickle
# ax = pickle.load(open('justKADR_Vclamp.pickle', 'rb'))
# plt.show()
