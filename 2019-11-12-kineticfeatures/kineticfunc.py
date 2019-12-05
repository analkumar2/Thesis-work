#exec(open('kineticfunc.py').read())

import numpy as np
import matplotlib.pyplot as plt

def kineticfunc1(Gbar, minfV, mtauV, hinfV, htauV, minfm80, hinfm80):
    #Assuming that at time=0, the channel is at steady state at -80mV.

    dt = 0.00001
    tlist = np.arange(0,0.5,dt)
    mlist = [minfm80 + dt*(minfV-minfm80)/mtauV]
    hlist = [hinfm80 + dt*(hinfV-hinfm80)/htauV]
    for t in tlist[1:]:
        mlist.append(mlist[-1]+dt*(minfV-mlist[-1])/mtauV)
        hlist.append(hlist[-1]+dt*(hinfV-hlist[-1])/htauV)
    G = Gbar*np.array(mlist)**1*np.array(hlist)
    return G

at0 = kineticfunc1(1, 0.5,0.005, 0.5,0.050, 0,1)
plt.figure(1)
plt.plot(at0)

def kineticfunc2(t, Gbar, minfV, mtauV, hinfV, htauV, min, hin):
    #Assuming that at time=0, the channel is at steady state at -80mV.
    m = minfV + (min-minfV)*np.exp(-t/mtauV)
    h = hinfV + (hin-hinfV)*np.exp(-t/htauV)
    return Gbar*m*h

at0 = [kineticfunc2(t, 1, 0.5,0.005, 0.5,0.050, 0,1) for t in np.arange(0,0.5,1e-5)]
plt.figure(2)
plt.plot(at0)
plt.show()
