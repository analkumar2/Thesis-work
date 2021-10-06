import moose
import rdesigneur as rd
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import brute_curvefit
from pprint import pprint
import MOOSEModel_17_somamulti as mm
from Combined100models import Models as Models

stim_start = 1

# rdes = rd.rdesigneur(
# 	cellProto=[['Compartments_v7.swc', 'elec'],],
# 	passiveDistrib=[["#", "RM", "1", "CM", "0.01", "RA", "1", "Em", "-0.065"],],
# 	stimList=[["soma", "1", ".", "inject", f"(t>=1 && t<=1.5) ? {-25e-12} : 0"],],
# 	plotList=[["soma", "1", ".", "Vm", "Soma membrane potential MOOSE"],],
# 	)

# rdes.buildModel()

# # Setup clock table to record time
# clk = moose.element("/clock")
# moose.Neutral("Graphs")
# plott = moose.Table("/Graphs/plott")
# moose.connect(plott, "requestOut", clk, "getCurrentTime")

# moose.reinit()
# moose.start(2)

# Vmvec = moose.element("/model/graphs/plot0").vector
# tvec = moose.element("/Graphs/plott").vector

tvec, Vmvec, Cavec = mm.runModel(Models["Model4"], -25e-12)


##### Calculating Rin and Cin ################
Erest = np.median(Vmvec[(tvec<1) & (tvec>0.5)])
print(f"{Erest = }")

def chargingm25(t, R1, R2, tau1, tau2):
        return Erest - R1 * 25e-12 * (1 - np.exp(-t / tau1)) - R2 * 25e-12 * (1 - np.exp(-t / tau2))

tempv = Vmvec[(tvec > stim_start) & (tvec < stim_start + 0.1)]
RCfitted_chm25, errorm25 = brute_curvefit.brute_scifit(
    chargingm25,
    np.linspace(0, 0.1, len(tempv)),
    tempv,
    restrict=[[5e6, 5e6, 0, 0], [1000e6, 1000e6, 0.1, 0.1]],
    ntol=1000,
    printerrors=False,
)
Rin = RCfitted_chm25[0]+RCfitted_chm25[1]
if RCfitted_chm25[2]>RCfitted_chm25[3]:
    Cin = RCfitted_chm25[2]/RCfitted_chm25[0]
else:
    Cin = RCfitted_chm25[3]/RCfitted_chm25[1]

print(f"{Rin = }, {Cin = }")
################################################

#### Calculating surface area #########
Sarea = 0
for compt in moose.wildcardFind("/model/elec/#[CLASS==ZombieCompartment]"):
    SA = compt.length*np.pi*compt.diameter
    print(compt,SA)
    Sarea = Sarea+SA
print(f'Total Surface area = {Sarea}')
#########################################

plt.plot(tvec, Vmvec)
plt.show()