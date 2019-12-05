#exec(open('act_deact_HH_tau_simpcheck.py').read())

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


def m(t,minff,minfi,mtauf):
    return minff + (minfi-minff)*np.exp(-t/mtauf)

t=np.linspace(0,0.400,1000)
plt.plot(t,m(t,0.6,0,0.005)*m(t,0.3,1,0.100), label='Activation')
plt.plot(t,m(t,0.6,1,0.005)*m(t,0.3,0,0.100), label='Deactivation')
plt.legend()
plt.title('An artificial channel subjected to activation and deactivation protocols')
plt.xlabel('Time (s)')
plt.ylabel('Normalized conductance')
plt.show()

##Fitting
def actfunc(t, inf, act_tau):
    return inf*(1 - np.exp(-t/act_tau))

def deactfunc(t, C, A, deact_tau):
    return (C + A*np.exp(-t/deact_tau))

def inactfunc(t, C, A, inact_tau):
    return (C + A*np.exp(-t/inact_tau))

t=np.linspace(0,0.400,1000)
actdata = m(t,0.6,0,0.005)*m(t,0.3,1,0.100)
deactdata = m(t,0.6,1,0.005)*m(t,0.3,0,0.100)
actIpeak_idx = actdata.argmax()
fittedparam_act, cov = curve_fit(actfunc, t[:actIpeak_idx], actdata[:actIpeak_idx], p0=[0.6,0.010], bounds=([0,0],[1,1]))
fittedparam_deact, cov = curve_fit(deactfunc, t, deactdata, p0=[0.2,-0.2,0.010], bounds=([-1,-1,0],[1,1,1]))
fittedparam_inact, cov = curve_fit(inactfunc, t[actIpeak_idx:], actdata[actIpeak_idx:], p0=[0.2,-0.2,0.010], bounds=([-1,-1,0],[1,1,5]))

print(fittedparam_act)
print(fittedparam_deact)
print(fittedparam_inact)
# t = t[:actIpeak_idx]
plt.plot(t, actfunc(t, *fittedparam_act), label='fitted')
plt.plot(t, actdata, label='original')
plt.legend()
plt.title('Activation')
plt.xlabel('Time (s)')
plt.ylabel('Normalized conductance')
plt.show()

plt.plot(t, deactfunc(t, *fittedparam_deact), label='fitted')
# plt.plot(t, deactfunc(t,*[0,0.2,0.010]))
plt.plot(t, deactdata, label='original')
plt.legend()
plt.title('Deactivation')
plt.xlabel('Time (s)')
plt.ylabel('Normalized conductance')
plt.show()

plt.plot(t, inactfunc(t, *fittedparam_inact), label='fitted')
# plt.plot(t, inactfunc(t,*[0,0.2,0.010]))
plt.plot(t, actdata, label='original')
plt.legend()
plt.title('Inactivation')
plt.xlabel('Time (s)')
plt.ylabel('Normalized conductance')
plt.show()

plt.plot(t, m(t,fittedparam_act[0],0,fittedparam_act[1])*m(t,fittedparam_inact[0],1,fittedparam_inact[2]), label='fitted')
# plt.plot(t, inactfunc(t,*[0,0.2,0.010]))
plt.plot(t, actdata, label='original')
plt.legend()
plt.title('Inactivation')
plt.xlabel('Time (s)')
plt.ylabel('Normalized conductance')
plt.show()
