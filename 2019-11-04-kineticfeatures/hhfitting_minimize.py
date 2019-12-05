#exec(open('hhfitting_minimize.py').read())
#Using fmin
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from scipy.optimize import minimize
from scipy.optimize import curve_fit
from scipy.optimize import Bounds


# We do it for cell rCell10070.nwb
# Take from 993 to 5987 idx which is 100.3ms to 599.7ms
# First we extract the trace, then fit individual trace and then in a loop collect for all inf and tau values.
# In another file we check if the model also fits to the deactivation traces.

reader = h5py.File('../../Raw_data/Channelpedia/Kv1.1/DataKv1.1RatCHO/rCell10070.nwb')
data = reader['acquisition']['timeseries']['Activation']['repetitions']['repetition2']['data']
leak = data[993:5987,4]
datamleak = np.transpose(data[993:5987,4:15]) - np.transpose(leak)
t = np.linspace(0,0.4994, len(datamleak[1]))
# plt.figure(1)
# plt.plot(np.transpose(np.transpose(data[993:5987,4:15]) - np.transpose(leak)))
# plt.figure(2)
# plt.plot(data[993:5987,4:15])
# plt.show()

def kineticfunc_array(Gbar, minfV, mtauV, hinfV, htauV, min, hin, mpow, hpow):
    #Assuming that at time=0, the channel is at steady state at -80mV.
    t = np.linspace(0,0.4994, len(datamleak[1]))
    m = minfV + (min-minfV)*np.exp(-t/mtauV)
    h = hinfV + (hin-hinfV)*np.exp(-t/htauV)
    return Gbar*m**mpow*h**hpow

def error(x):
    Gbar, minfV, mtauV, hinfV, htauV, min, hin, mpow, hpow = x
    return np.sum((kineticfunc_array(Gbar, minfV, mtauV, hinfV, htauV, min, hin, mpow, hpow) - datamleak[10])**2)

bounds = Bounds([0,0,0.0001,0,0.0001,0,0,0,0], [10,1,1,1,1,1,1,4,4])
# minimum = minimize(error, [3,1,0.005,0.05,0.050,0,1,1,1], method='L-BFGS-B', bounds=bounds)
# minimum = minimize(error, [3,1,0.005,0.05,0.050,0,1,1,1], method='TNC', bounds=bounds)
# minimum = minimize(error, [3,1,0.005,0.05,0.050,0,1,1,1], method='SLSQP', bounds=bounds)
minimum = minimize(error, [3,1,0.005,0.05,0.050,0,1,1,1], method='trust-constr', bounds=bounds)


print(minimum.x)
plt.plot(t,datamleak[10])
plt.plot(t,kineticfunc_array(*minimum.x))
plt.show()
