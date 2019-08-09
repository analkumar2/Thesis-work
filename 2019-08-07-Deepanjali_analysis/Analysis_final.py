# exec(open('Deepanjali data/Analysis_final.py').read())

# Extracts all the needed features from all the files in a folder and stores it in a csv file. Uses allensdk package.

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
import csv
import pandas as pd
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

# //////////////////////////////////////////////////////////////////
foldername = 'Deepanjali data/WT step input cells' #Name of the folder where files are stored
fskip = ['Cell 2 of 21_3_2017.abf'] #List of files that need to be skipped
seg_nol = [10,16] #list of segment numbers to evaluate the features for. Note that 0 refers to -100pA injection. 10 to 150pA, 16 to 300pA, 20 to 400pA.
# //////////////////////////////////////////////////////////////////

def singlefeature(seg_nol,reader,filename,stim_start,stim_end):
    '''This function is used to extract features from a neo.io.axonio.AxonIO object and extracts features from it.
    seg_nol takes in a list of segment numbers. reader is the objecct. filename is the name of the neuron. stim_start and stim_end are when the current injection was started and ended.
    Outputs a pandas dataframe with all the features'''
    features = {}
    for seg_no in seg_nol:
        Vtrace = reader.read_block().segments[seg_no].analogsignals[0]
        #Extra
        features['Cell name'] = filename
        features['Sampling rate'] = Vtrace.sampling_rate
        features['stim_start'] = stim_start
        features['stim_end'] = stim_end
        Inputcurr = seg_no*25e-12 - 100e-12 #in A

        i = np.zeros(int((Vtrace.t_stop - Vtrace.t_start)*Vtrace.sampling_rate))
        i[int(stim_start*Vtrace.sampling_rate):int(stim_end*Vtrace.sampling_rate)] = Inputcurr*1e12
        i = np.array(i) #in pA
        v = np.array([float(V) for V in Vtrace]) #in mV
        t = np.linspace(0,float(Vtrace.t_stop - Vtrace.t_start), int((Vtrace.t_stop - Vtrace.t_start)*Vtrace.sampling_rate))
        t = np.array(t) # in s

        # plt.plot(t,v, label = f'{filename}')
        # plt.legend()
        # plt.show()

        sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, filter=float(Vtrace.sampling_rate)/2500)
        sweep_ext.process_spikes()


        # E_rest
        features[f'E_rest_{Inputcurr*1e12}'] = np.nanmean(v[t<=stim_start])*1e-3

        # Input resistance #steady-state V is average of last 100ms of the current clamp duration
        Vtracep25 = reader.read_block().segments[5].analogsignals[0]
        Vtracep25 = np.array([float(V) for V in Vtracep25]) #in mV
        Vtracen25 = reader.read_block().segments[3].analogsignals[0]
        Vtracen25 = np.array([float(V) for V in Vtracen25]) #in mV
        str_ind = (np.abs(t-stim_end+100e-3)).argmin()
        stp_ind = (np.abs(t-stim_end)).argmin()
        features[f'Rinput'] = (-1*np.nanmean(Vtracen25[str_ind:stp_ind]) + np.nanmean(Vtracep25[str_ind:stp_ind]))*1e-3/2/25e-12

        # Total capacitance
        Vtracep25_choppped = Vtracep25[:stp_ind]
        Vtracen25_choppped = Vtracen25[:stp_ind]
        vp63 = (np.nanmean(Vtracep25[str_ind:stp_ind]) - features[f'E_rest_{Inputcurr*1e12}'])*0.63 + features[f'E_rest_{Inputcurr*1e12}']
        vn63 = (np.nanmean(Vtracen25[str_ind:stp_ind]) - features[f'E_rest_{Inputcurr*1e12}'])*0.63 + features[f'E_rest_{Inputcurr*1e12}']
        tau = (t[(np.abs(Vtracep25_choppped-vp63)).argmin()] - stim_start + t[(np.abs(Vtracep25_choppped-vp63)).argmin()] - stim_start)/2
        tauinv = (t[len(Vtracep25_choppped)-(np.abs(Vtracep25_choppped[stp_ind::-1]-vp63)).argmin()] - stim_start + t[len(Vtracep25_choppped)-(np.abs(Vtracep25_choppped[::-1]-vp63)).argmin()] - stim_start)/2
        features[f'Cm'] = (tau+tauinv)/features[f'Rinput']/2

        # AP1_amp
        features[f'AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("peak_v")[0]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #APp_amp
        features[f'APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("peak_v")[-2]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #AP1_width
        features[f'AP1_width_{Inputcurr*1e12}'] = sweep_ext.spike_feature("width")[0]

        #APp_width
        features[f'APp_width_{Inputcurr*1e12}'] = sweep_ext.spike_feature("width")[-2]

        #AP1_thresh
        features[f'AP1_thresh_{Inputcurr*1e12}'] = sweep_ext.spike_feature("threshold_v")[0]*1e-3

        #APp_thresh
        features[f'APp_thresh_{Inputcurr*1e12}'] = sweep_ext.spike_feature("threshold_v")[-2]*1e-3

        #AP1_lat
        features[f'AP1_lat_{Inputcurr*1e12}'] = sweep_ext.spike_feature("threshold_t")[0] - stim_start

        #ISI1
        features[f'ISI1_{Inputcurr*1e12}'] = sweep_ext.spike_feature("peak_t")[1] - sweep_ext.spike_feature("peak_t")[0]

        #ISIl
        features[f'ISIl_{Inputcurr*1e12}'] = sweep_ext.spike_feature("peak_t")[-1] - sweep_ext.spike_feature("peak_t")[-2]

        #ISIavg
        pt = sweep_ext.spike_feature("peak_t")
        features[f'ISIavg_{Inputcurr*1e12}'] = np.nanmean([s-f for s,f in zip(pt[1:],pt[:-1])])

        #freq
        features[f'freq_{Inputcurr*1e12}'] = len(sweep_ext.spike_feature("peak_t"))/(stim_end - stim_start)

        #Adptn_id = 1-ISI1/ISIl
        features[f'Adptn_id_{Inputcurr*1e12}'] = 1 - features[f'ISI1_{Inputcurr*1e12}']/features[f'ISIl_{Inputcurr*1e12}']

        #fAHP_AP1_amp
        features[f'fAHP_AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("fast_trough_v")[0]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #fAHP_APp_amp
        features[f'fAHP_APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("fast_trough_v")[-2]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #mAHP_AP1_amp
        features[f'mAHP_AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("slow_trough_v")[0]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #mAHP_APp_amp
        features[f'mAHP_APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("slow_trough_v")[-2]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #mAHP_AP1_dur
        features[f'mAHP_AP1_dur_{Inputcurr*1e12}'] = (sweep_ext.spike_feature("slow_trough_t")[0] - sweep_ext.spike_feature("peak_t")[0])/features[f'ISI1_{Inputcurr*1e12}']

        #mAHP_APp_dur = mAHP of second last spike (penultimate)
        features[f'mAHP_APp_dur_{Inputcurr*1e12}'] = (sweep_ext.spike_feature("slow_trough_t")[-2] - sweep_ext.spike_feature("peak_t")[-2])/features[f'ISIl_{Inputcurr*1e12}']

        #ADP_AP1_amp
        features[f'ADP_AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("adp_v")[0]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #ADP_APp_amp
        features[f'ADP_APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("adp_v")[-2]*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #mAHP_stimend_amp = within 50ms
        features[f'mAHP_stimend_amp_{Inputcurr*1e12}'] = np.min(v[int((stim_end)*Vtrace.sampling_rate):int((stim_end + 50e-3)*Vtrace.sampling_rate)])*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

        #sAHP_stimend_amp = within 200ms
        features[f'sAHP_stimend_amp_{Inputcurr*1e12}'] = np.min(v[int((stim_end)*Vtrace.sampling_rate):int((stim_end + 200e-3)*Vtrace.sampling_rate)])*1e-3 - features[f'E_rest_{Inputcurr*1e12}']

    return features


Sno = 0
features_pd = pd.DataFrame()
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
        reader  = AxonIO(filename='Deepanjali data/WT step input cells/'+filename)
    else:
        continue

    features = singlefeature(seg_nol,reader,filename,stim_start,stim_end)
    features_pd = features_pd.append(pd.DataFrame(features,index = [Sno]))
    Sno = Sno +1
    print(Sno)
    print(features['Cm'])
    print(features['Rinput'])

# features_pd.to_csv('Deepanjali data/WT step input cells/features.csv', sep='\t')
