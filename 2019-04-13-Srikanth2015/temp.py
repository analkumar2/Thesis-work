# exec(open('Srikanth2015/temp.py').read())

# First do nrniv -python in the terminal with cwd in the folder where mod files are there.
# Then execute this file
# exec(open('Combe2018/Analysis_Ae.py').read())
import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
import csv
import pandas as pd
from neuron import h,gui
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor



def get_features(t,v, nameit):
    #v is in V
    #t is in s
    stim_start = 1000e-3
    stim_end = 1500e-3
    features = {}
    #Extra
    features['Cell name'] = nameit
    features['stim_start'] = stim_start
    features['stim_end'] = stim_end
    Inputcurr = 150e-12 #in A

    v = v*1e3 #in mV
    t = t #in s
    start_idx = (np.abs(t - stim_start)).argmin()
    end_idx = (np.abs(t - stim_end)).argmin()
    i = np.zeros(len(t))
    i[start_idx:end_idx] = Inputcurr*1e12 #in pA

    sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, filter = 4.9)
    sweep_ext.process_spikes()

    # E_rest
    features[f'E_rest_{Inputcurr*1e12}'] = np.nanmean(v[i==0])

    # AP1_amp
    features[f'AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("peak_v")[0] - features[f'E_rest_{Inputcurr*1e12}']

    #APp_amp
    features[f'APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("peak_v")[-2] - features[f'E_rest_{Inputcurr*1e12}']

    #AP1_width
    features[f'AP1_width_{Inputcurr*1e12}'] = sweep_ext.spike_feature("width")[0]

    #APp_width
    features[f'APp_width_{Inputcurr*1e12}'] = sweep_ext.spike_feature("width")[-2]

    #AP1_thresh
    features[f'AP1_thresh_{Inputcurr*1e12}'] = sweep_ext.spike_feature("threshold_v")[0]

    #APp_thresh
    features[f'APp_thresh_{Inputcurr*1e12}'] = sweep_ext.spike_feature("threshold_v")[-2]

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
    features[f'fAHP_AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("fast_trough_v")[0] - features[f'E_rest_{Inputcurr*1e12}']

    #fAHP_APp_amp
    features[f'fAHP_APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("fast_trough_v")[-2] - features[f'E_rest_{Inputcurr*1e12}']

    #mAHP_AP1_amp
    features[f'mAHP_AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("slow_trough_v")[0] - features[f'E_rest_{Inputcurr*1e12}']

    #mAHP_APp_amp
    features[f'mAHP_APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("slow_trough_v")[-2] - features[f'E_rest_{Inputcurr*1e12}']

    #mAHP_AP1_dur
    features[f'mAHP_AP1_dur_{Inputcurr*1e12}'] = (sweep_ext.spike_feature("slow_trough_t")[0] - sweep_ext.spike_feature("peak_t")[0])/features[f'ISI1_{Inputcurr*1e12}']

    #mAHP_APp_dur = mAHP of second last spike (penultimate)
    features[f'mAHP_APp_dur_{Inputcurr*1e12}'] = (sweep_ext.spike_feature("slow_trough_t")[-2] - sweep_ext.spike_feature("peak_t")[-2])/features[f'ISIl_{Inputcurr*1e12}']

    #ADP_AP1_amp
    features[f'ADP_AP1_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("adp_v")[0] - features[f'E_rest_{Inputcurr*1e12}']

    #ADP_APp_amp
    features[f'ADP_APp_amp_{Inputcurr*1e12}'] = sweep_ext.spike_feature("adp_v")[-2] - features[f'E_rest_{Inputcurr*1e12}']

    #mAHP_stimend_amp = within 50ms
    end50_idx = (np.abs(t - stim_end - 50e-3)).argmin()
    print(end50_idx)
    features[f'mAHP_stimend_amp_{Inputcurr*1e12}'] = np.min(v[end_idx:end50_idx]) - features[f'E_rest_{Inputcurr*1e12}']

    #sAHP_stimend_amp = within 200ms
    end200_idx = (np.abs(t - stim_end - 200e-3)).argmin()
    features[f'sAHP_stimend_amp_{Inputcurr*1e12}'] = np.min(v[end_idx:end200_idx]) - features[f'E_rest_{Inputcurr*1e12}']

    return [features, t, v]
