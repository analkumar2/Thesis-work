## exec(open('features.py').read())

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
from plotexp
import os
import csv
import pandas as pd
import brute_curvefit

from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

def expfeatures(cellpath,stim_start,stim_end):
    reader  = AxonIO(filename=cellpath)
    Vtracem25 = reader.read_block().segments[3].analogsignals[0]
    Vtrace25 = reader.read_block().segments[5].analogsignals[0]
    Vtrace150 = reader.read_block().segments[10].analogsignals[0]
    Vtrace300 = reader.read_block().segments[16].analogsignals[0]
    features = {}
    features['Cell name'] = cellpath.split('/')[-1]
    features['Sampling rate'] = Vtracem25.sampling_rate
    features['stim_start'] = stim_start
    features['stim_end'] = stim_end

    Inputcurr = 3*25e-12 - 100e-12 #in A
    i = np.zeros(int((Vtracem25.t_stop - Vtracem25.t_start)*Vtracem25.sampling_rate))
    i[int(stim_start*Vtracem25.sampling_rate):int(stim_end*Vtracem25.sampling_rate)] = Inputcurr
    i = np.array(i) #in A
    v = np.array([float(V) for V in Vtracem25])*1e-3 #in V
    t = np.linspace(0,float(Vtracem25.t_stop - Vtracem25.t_start), int((Vtracem25.t_stop - Vtracem25.t_start)*Vtracem25.sampling_rate))
    t = np.array(t) # in s

    sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, filter=float(Vtracem25.sampling_rate)/2500)
    sweep_ext.process_spikes()

    features[f'E_rest'] = float(np.mean(reader.read_block().segments[4].analogsignals[0]))*1e-3

    def charging25(t, R,C):
        return features[f'E_rest'] + R*25e-12*(1+np.exp(-t/R/C))

    def chargingm25(t, R,C):
        return features[f'E_rest'] - R*25e-12*(1+np.exp(-t/R/C))

    def discharging25(t, R,C):
        return features[f'E_rest'] + R*25e-12*(np.exp(-t/R/C))

    def dischargingm25(t, R,C):
        return features[f'E_rest'] + R*25e-12*(np.exp(-t/R/C))

    plt.plot(t, Vtrace25)
    plt.plot(t[(t>stim_start) & (t<stim_start+0.1)], np.array(np.ravel(Vtrace25))[(t>stim_start) & (t<stim_start+0.1)])
    plt.show()

    tempv = np.array(np.ravel(Vtrace25))[(t>stim_start) & (t<stim_start+0.1)]
    print(tempv)
    RCfitted_c25, error = brute_curvefit.brute_scifit(charging25, np.linspace(0,0.1,len(tempv)), tempv, restrict=[[50e6,50e-10],[250e6,250e-10]], )