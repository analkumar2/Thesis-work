# exec(open('Passive properties/featurecalc.py').read())
# only for RC circuit type passive models. Only positive currents allowed


import numpy as np
import matplotlib.pyplot as plt

def getRin(plotobject, currAmp):
    E_rest = min(plotobject.vector)
    Vmax = max(plotobject.vector)
    return abs((Vmax-E_rest)/currAmp)

def getCm(plotobject, currAmp,stimstart_time):
    V = plotobject.vector[:int(len(plotobject.vector)*0.75)]
    E_rest = min(plotobject.vector)
    Vmax = max(plotobject.vector)
    V63 = abs(Vmax-E_rest)*0.63+E_rest
    t63 = plotobject.dt*(np.abs(V63-V)).argmin()
    tau = t63-stimstart_time
    return tau/abs((Vmax-E_rest)/currAmp)

print(getRin(moose.element('/model/graphs/plot0'),10e-12))
print(getCm(moose.element('/model/graphs/plot0'),10e-12,1))
