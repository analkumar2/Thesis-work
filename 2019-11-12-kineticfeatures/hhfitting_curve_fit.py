#exec(open('hhfitting_curve_fit.py').read())
#Using curve_fit
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from scipy.optimize import minimize
from scipy.optimize import curve_fit


# We do it for cell rCell10070.nwb
# Take from 993 to 5987 idx which is 100.3ms to 599.7ms
# First we extract the trace, then fit individual trace and then in a loop collect for all inf and tau values.
# In another file we check if the model also fits to the deactivation traces.

reader = h5py.File('../../Raw_data/Channelpedia/Kv1.1/DataKv1.1RatCHO/rCell10070.nwb')
data = reader['acquisition']['timeseries']['Activation']['repetitions']['repetition2']['data']
leak = data[993:5987,4]*1e-9
G = (np.transpose(data[993:5987,4:15])*1e-9 - np.transpose(leak))/(np.arange(-0.050,0.060,0.010)+0.0962)[:,None]
# plt.figure(1)
# plt.plot(np.transpose(G))
# plt.figure(2)
# plt.plot(data[993:5987,4:15])
# plt.show()

def kineticfunc(t, Gbar, minfV, mtauV, hinfV, htauV):
    #Assuming that at time=0, the channel is at steady state at -80mV.
    min, hin = 0,1
    m = minfV + (min-minfV)*np.exp(-t/mtauV)
    h = hinfV + (hin-hinfV)*np.exp(-t/htauV)
    return Gbar*m**1*h**1*1e8

levelnum = 10
fittedparam, cov = curve_fit(kineticfunc, np.arange(0,len(G[levelnum])/10000,1/10000), G[levelnum]*1e8, p0=[2,1,0.005,0.25,0.046], bounds=([0,0,0.0001,0,0.0001],[100,1,1,1,1]))

print(fittedparam)
plt.plot(np.arange(0,len(G[levelnum])/10000,1/10000), G[levelnum])
plt.plot(np.arange(0,len(G[levelnum])/10000,1/10000), kineticfunc(np.arange(0,len(G[levelnum])/10000,1/10000), *fittedparam)*1e-8)
plt.show()
