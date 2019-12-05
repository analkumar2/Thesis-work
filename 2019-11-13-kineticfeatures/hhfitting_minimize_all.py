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

levelnum = 13
reader = h5py.File('../../Raw_data/Channelpedia/Kv1.1/DataKv1.1RatCHO/rCell10070.nwb')
data = reader['acquisition']['timeseries']['Activation']['repetitions']['repetition2']['data']
leak = np.mean(data[993:5987,4]*1e-9)
G = (np.transpose(data[993:5987,:])*1e-9 - np.transpose(leak))/(np.arange(-0.090,0.090,0.010)+0.0962)[:,None]
Gnorm = G/np.max(G)
t = np.arange(0,len(Gnorm[levelnum])/10000,1/10000)
# plt.figure(1)
# plt.plot(np.transpose(G))
# plt.figure(2)
# plt.plot(data[993:5987,4:15])
# plt.show()

def kineticfunc_array(minfV, mtauV, hinfV, htauV, min, hin, mpow, hpow):
    #Assuming that at time=0, the channel is at steady state at -80mV.
    m = minfV + (min-minfV)*np.exp(-t/mtauV)
    h = hinfV + (hin-hinfV)*np.exp(-t/htauV)
    return m**mpow*h**hpow

def error(x):
    minfV, mtauV, hinfV, htauV = x
    min, hin, mpow, hpow = 0,1,1,1
    return np.sum((kineticfunc_array(minfV, mtauV, hinfV, htauV, min, hin, mpow, hpow) - Gnorm[levelnum])**2)

bounds = Bounds([0,0.00001,0,0.0001],[1,1,1,1])
#bb = [(0,100e-8),(0.1),(1e-5,1),(0,1),(1e-4,1)]
# minimum = minimize(error, [3,1,0.005,0.05,0.050,0,1,1,1], method='L-BFGS-B', bounds=bounds)
# minimum = minimize(error, [3,1,0.005,0.05,0.050,0,1,1,1], method='TNC', bounds=bounds)
# minimum = minimize(error, [1,0.005,0.05,0.050,0,1], method='Nelder-Mead', bounds=bounds)
# minimum = minimize(error, [5e-8,1,0.0005,0.25,0.1], method='trust-constr', bounds=bounds)


for level in np.arange(-0.040,0.060,0.010):
    levelnum = int((level+0.060)/0.010)
    minimum = minimize(error, [1,0.005,0.05,0.050], method='Nelder-Mead', bounds=bounds)
    print(minimum.x)
    plt.plot(t,Gnorm[levelnum])
    plt.plot(t,kineticfunc_array(*minimum.x, *[0,1,1,1]))
    # plt.plot(t,kineticfunc_array(*[10e-8,1,0.0005,0.25,0.1], *[0,1,1,1]))
plt.show()
