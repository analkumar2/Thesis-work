#exec(open('hhfitting_brute_curvefit.py').read())
#Using brute_curvefit
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import brute_curvefit


# We do it for cell rCell10070.nwb
# Take from 993 to 5987 idx which is 100.3ms to 599.7ms
# First we extract the trace, then fit individual trace and then in a loop collect for all inf and tau values.
# In another file we check if the model also fits to the deactivation traces.

# levelnum = 13
reader = h5py.File('../../Raw_data/Channelpedia/Kv1.1/DataKv1.1RatCHO/rCell2251.nwb', 'r')
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

mhin = Gnorm[levelnum][0]

def kineticfunc(t, minfV, mtauV, hinfV, htauV, min):
    #Assuming that at time=0, the channel is at steady state at -80mV.
    # min, hin = 0,1
    hin = mhin/min
    m = minfV + (min-minfV)*np.exp(-t/mtauV)
    h = hinfV + (hin-hinfV)*np.exp(-t/htauV)
    return m**1*h**1

fittedparam, error = brute_curvefit.brute_then_scipy(kineticfunc, t, Gnorm[levelnum],  restrict=[[0,0,0,0,0],[1,1,1,1,1]], ntol=100000)

print(fittedparam)
plt.plot(t, Gnorm[levelnum])
plt.plot(t, kineticfunc(t, *fittedparam))
plt.show()

# # Just plotting the brute force version of plots
# plt.plot(t, Gnorm[levelnum])
# paramsfitted,errors = brute_curvefit.bruteforce(kineticfunc, t, Gnorm[levelnum],  restrict=[[0,0,0,0,0],[1,1,1,1,1]], ntol=100, returnnfactor = 1)
# for param in paramsfitted:
#     plt.plot(t,kineticfunc(t, *param), label='fitted')
# plt.legend()
# plt.show()

## Just plotting for scipy version
# fittedparam, error = brute_curvefit.scipy_fit(kineticfunc, t, Gnorm[levelnum],  restrict=[[0,0,0,0,0],[1,1,1,1,1]], p0list=[[1,0.010,0,0.1,0.1]])
#
# print(fittedparam)
# plt.plot(t, Gnorm[levelnum])
# plt.plot(t, kineticfunc(t, *fittedparam))
# plt.show()
