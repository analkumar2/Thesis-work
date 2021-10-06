import numpy as np
import matplotlib.pyplot as plt
import moose
import MOOSEModel_17_somamulti as mm
from pprint import pprint
from copy import deepcopy
from Combined100models import Models as Models

modeldict = Models["Model4"]

hold = np.arange(-0.055, 0.045, 0.010)
exec(open('PlotexpSKcurr.py').read()) #import exp apamin protocol data
for j in WT_old.keys():
    WT_old[j] = np.array(WT_old[j])*1e-12
Wt_old10 = np.percentile(np.array(list(WT_old.values())), 10, 0)
Wt_old90 = np.percentile(np.array(list(WT_old.values())), 90, 0)
# Wt_old10 = np.min(np.array(list(WT_old.values())), 0)
# Wt_old90 = np.max(np.array(list(WT_old.values())), 0)

apacurrent_list = []
K_SKcurrent_list = []

for holdV in hold:
    print(holdV, modeldict["Parameters"]["Channels"]["K_SK_Chan"]["gbar"])
    tbAp, IbAp, Ca = mm.runModel(
        modeldict, vClamp=f"-0.055 + (t>1 && t<1.1)*-0.010 + (t>1.1 && t<1.9)*{holdV+0.055}", refreshKin=False,
    )
    modeldict_temp = deepcopy(modeldict)
    modeldict_temp["Parameters"]["Channels"]["K_SK_Chan"]["gbar"] = 0
    taAp, IaAp, Ca = mm.runModel(
        modeldict_temp,
        vClamp=f"-0.055 + (t>1 && t<1.1)*-0.010 + (t>1.1 && t<1.9)*{holdV+0.055}", refreshKin=False,
    )
    apacurrent_list.append(
        IbAp[np.argmin(np.abs(1.925 - tbAp))]
        - IaAp[np.argmin(np.abs(1.925 - taAp))]
    )

print(apacurrent_list)
plt.plot(hold, apacurrent_list, color='red')
plt.plot(hold, Wt_old10, label='10th and 90th percentile of experimental data', color='black')
plt.plot(hold, Wt_old90, color='black')
plt.legend()
plt.title('apamin Vclamp protocol current')
plt.xlabel('Preholding potential (V)')
plt.ylabel('apamin sensitive current (A)')
plt.show()