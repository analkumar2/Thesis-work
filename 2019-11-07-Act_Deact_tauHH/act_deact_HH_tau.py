#exec(open('act_deact_HH_tau.py').read())

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
from scipy import signal
from scipy.optimize import minimize
from scipy.optimize import curve_fit

#Activation
try:
    [moose.delete(x) for x in ['/model', '/library']]
    # moose.delete('/model')
except:
    pass


Vlevel = 0
rdes = rd.rdesigneur(
    elecDt = 0.00005,
    chanProto = [['make_HH_K()', 'K']],
    chanDistrib = [['K', 'soma', 'Gbar', '360' ]],
    stimList = [['soma', '1', '.', 'vclamp', f'-0.080 + (t>0.3 && t<0.5) * {Vlevel+0.080}' ]],
    plotList = [
        # ['soma', '1', '.', 'Vm', 'Membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
    ]
)
rdes.buildModel()
#Changing vclamp parameters
try:
    moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecDt
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecDt
    moose.element( '/model/elec/soma/vclamp' ).ti = elecDt
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass

moose.reinit()
moose.start(0.6)
# rdes.display()

actdata = moose.element('/model/graphs/plot0').vector


#deactivation
try:
    [moose.delete(x) for x in ['/model', '/library']]
    # moose.delete('/model')
except:
    pass


Vlevel = 0
rdes = rd.rdesigneur(
    elecDt = 0.00005,
    chanProto = [['make_HH_K()', 'K']],
    chanDistrib = [['K', 'soma', 'Gbar', '360' ]],
    stimList = [['soma', '1', '.', 'vclamp', f'0.070 + (t>0.3 && t<0.5) * {Vlevel-0.070}' ]],
    plotList = [
        # ['soma', '1', '.', 'Vm', 'Membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
    ]
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
moose.start(0.6)
# rdes.display()
deactdata = moose.element('/model/graphs/plot0').vector

##Plotting
plt.plot(plott.vector, actdata, label='Clamped from -80mV to 0mV')
plt.plot(plott.vector, deactdata, label='Clamped from 70mV to 0mV')
plt.xlim(0.250,0.350)
plt.legend()
plt.show()

##Fitting
def actfunc(t, act_tau):
    return (1 - np.exp(-t/act_tau))*1e6

def deactfunc(t, C, A, deact_tau):
    return (C + A*np.exp(-t/deact_tau))*1e6

t = plott.vector
actIpeak_idx = (actdata[3010:4000]).argmax()
fittedparam_act, cov = curve_fit(actfunc, np.linspace(0,t[actIpeak_idx],actIpeak_idx), actdata[3010:3010+actIpeak_idx]*1e6, p0=[0.010], bounds=([0],[1]))

t = np.linspace(0,t[actIpeak_idx],actIpeak_idx)
plt.plot(t, actfunc(t, fittedparam_act[0])*1e-6)
plt.plot(t, actdata[3010:3010+actIpeak_idx])
plt.show()
