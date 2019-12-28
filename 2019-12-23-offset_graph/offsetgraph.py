# exec(open('offsetgraph.py').read())

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

# //////////////////////////////////////////////////////////////////
foldername = '../../Raw_data/Deepanjali_data/WT step input cells' #Name of the folder where files are stored
fskip = ['Cell 2 of 21_3_2017.abf'] #List of files that need to be skipped
seg_nol = [10,16] #list of segment numbers to evaluate the features for. Note that 0 refers to -100pA injection. 10 to 150pA, 16 to 300pA, 20 to 400pA.
# //////////////////////////////////////////////////////////////////

def somefeatures(currAmp,Vtrace,stim_start,stim_end,sampling_rate):
    features = {}
    # i = np.zeros(len(Vtrace))
    # i[int(stim_start*sampling_rate):int(stim_end*sampling_rate)] = currAmp
    # i = np.array(i) #in A
    # v = np.array([float(V) for V in Vtrace]) #in V
    t = np.linspace(0,len(Vtrace)/sampling_rate, len(Vtrace))
    t = np.array(t) # in s
    # sweep_ext = EphysSweepFeatureExtractor(t=t, v=v*1e3, i=i*1e12, filter=4.9)
    # sweep_ext.process_spikes()
    features[f'E_rest'] = np.nanmean(Vtrace[(t<=stim_start) & (t>=stim_start-0.2)])
    # features[f'freq'] = len(sweep_ext.spike_feature("peak_t"))/(stim_end - stim_start)
    features[f'freq'] = len(scs.find_peaks(Vtrace, height=-0.03, prominence=0.01)[0])/(stim_end - stim_start)
    features[f'offset'] = np.nanmin(Vtrace[(t<=stim_end-0.001) & (t>=stim_start+0.4)]) - features[f'E_rest']
    # features[f'find_peaks'] = scs.find_peaks(Vtrace, height=-0.03)
    
    return features

offset = {}

for filename in os.listdir(foldername):
    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf', 'Cell 2 of 19_10_2016', 'Cell 1 of 27_10_2016.abf', 'Cell 1 of 14_10_2016.abf', 'Cell 4 of 7_10_2016.abf', 'Cell 6 of 12_10_2016.abf', 'Cell 7 of 12_10_2016.abf']
    if filename in stim1391:
        stim_start = 139.1e-3
        stim_end = 639.1e-3
    else:
        stim_start = 81.4e-3 #in s
        stim_end = 581.4e-3 #in s

    if filename in fskip:
        print(f'{filename} skipped')
        continue

    if filename[-4:] == '.abf':
	    offset[filename] = []
	    for curramp in np.arange(-100e-12, 425e-12, 25e-12):
	    	t,v = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, curramp)
	    	features = somefeatures(curramp,v,stim_start,stim_end,len(t)/t[-1])
	    	offset[filename].append(features['offset'])
    else:
    	continue

offsetl = []
for off in offset.keys():
	offsetl.append(offset[off])
medianoffset = np.median(offsetl,0)
avgoffset = np.nanmean(offsetl,0)


I = np.arange(-100e-12, 425e-12, 25e-12)
plt.plot(I, avgoffset)
plt.plot(I, medianoffset)
plt.show()

plt.plot(I, np.transpose(offsetl))
plt.title('Offset of all the expereimental cells')
plt.xlabel('Injected current (A)')
plt.ylabel('Offset (V)')
plt.show()
