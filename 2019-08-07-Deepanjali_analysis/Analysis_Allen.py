# exec(open('Deepanjali data/Analysis_Allen.py').read())

# Extracts single feature using allensdk package

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
from allensdk.core.cell_types_cache import CellTypesCache

filename = 'Cell 6 of 171117.abf'
seg_no = 17 #0 is -100pA, 4 is 0pA, 20 is 400pA. Now extrapolate
stim_start = 81.4e-3 #in s
stim_stop = 581.4e-3 #in s

Actualstim_start = seg_no*2+stim_start
Actualstim_stop = seg_no*2+stim_stop
Inputcurr = seg_no*25 - 100 #in pA

reader  = AxonIO(filename='Deepanjali data/WT step input cells/'+filename)
Vtrace = reader.read_block().segments[seg_no].analogsignals[0]

i = np.zeros(int((Vtrace.t_stop - Vtrace.t_start)*Vtrace.sampling_rate))
i[int(stim_start*Vtrace.sampling_rate):int(stim_stop*Vtrace.sampling_rate)] = Inputcurr
i = np.array(i)
v = np.array([float(V) for V in Vtrace])
t = np.linspace(0,float(Vtrace.t_stop - Vtrace.t_start), int((Vtrace.t_stop - Vtrace.t_start)*Vtrace.sampling_rate))
t = np.array(t)

sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, filter=5)
sweep_ext.process_spikes()

sweep_ext.spike_feature_keys() #Lists all the features that can be extracted
sweep_ext.spike_feature("width") #Extracts AP width
