# exec(open('VisualInspectionModels.py').read())

import os
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib as mpl
# mpl.rcParams["savefig.directory"] = os.chdir(os.path.dirname(__file__))
import scipy.signal as scs
import MOOSEModel_2
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

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
    features[f'offset'] = np.nanmin(Vtrace[(t<=stim_end) & (t>=stim_start)]) - features[f'E_rest']

    return features

for filee in os.listdir('Modelparameters_temp'):
    print(filee)
    exec(open('Modelparameters_temp/'+filee).read())
    for model in Models.keys():
        print(model)
        print(Models[model]['Error'])
        tvec, Vmvec = MOOSEModel_2.runModel(Models[model], 150e-12)
        features = somefeatures(150e-12,Vmvec,1,1.5,len(tvec)/tvec[-1])
        print(features[f'freq'])
        print('############################################################')
        print('')
        if features[f'freq']>=6:
            plt.plot(tvec, Vmvec)
            plt.title(filee.split('_')[-1].split('.')[0] + '_' + model)
            plt.xlabel('Time (s)')
            plt.ylabel('Membrane potential (V)')
            plt.savefig('Output/'+filee.split('_')[-1].split('.')[0]+'_'+model)
            plt.clf()
            # plt.show()
