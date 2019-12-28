# exec(open('offset_representation.py').read())

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
import csv
import pandas as pd
import plotexp
import scipy.signal as scs
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

curramp = 150e-12
t,v = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+'Cell 6 of 171117.abf', curramp)

plt.plot(t,v, label = 'Cell 6 of 171117 150e-12')
plt.axvline(x=0.08135+0.400, color='r', ls='--')
plt.axvline(x=0.08135+0.499, color='r', ls='--')
off = np.nanmin(v[(t<=0.08135+0.499) & (t>=0.08135+0.400)])
offidx = np.argmin(v[(t<=0.08135+0.499) & (t>=0.08135+0.400)])
offidx = offidx+np.argmin(np.abs(t-0.08135-0.400))
Erest = np.nanmin(v[t<=0.0813])
plt.axhline(y=off, color='g', ls='--')
plt.axhline(y=Erest, color='g', ls='--')
plt.arrow(t[offidx], Erest, 0, -Erest+off, length_includes_head=True, head_width=0.02, head_length=0.004)
plt.arrow(t[offidx], off, 0, Erest-off, length_includes_head=True, head_width=0.02, head_length=0.004)
plt.text(0.53, -0.060, 'offset')
# plt.annotate(s='offset', xy=(t[offidx],off), xycoords='data', xytext=(t[offidx],Erest), textcoords='data', arrowprops=dict(arrowstyle='<->'))
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.show()