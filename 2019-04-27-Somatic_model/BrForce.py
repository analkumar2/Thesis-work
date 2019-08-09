#exec(open('Somatic model/BrForce.py').read())

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import csv
import random
import sys
import time


exec(open('Somatic model/feature_dict.py').read())

#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

def features(Vtrace):
    stim_start = preStimTime
    stim_end = preStimTime + injectTime
    features_df = feature_range_df.copy()
    Inputcurr = 150e-12 #in A

    v = np.array(Vtrace) #in mV
    t = np.linspace(0,runtime, len(Vtrace)) #in s
    start_idx = (np.abs(t - stim_start)).argmin()
    end_idx = (np.abs(t - stim_end)).argmin()
    i = np.zeros(len(t))
    i[start_idx:end_idx] = Inputcurr*1e12 #in pA

    sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, filter = 9.9)
    sweep_ext.process_spikes()

    features_df['raw'] = ''
    # E_rest
    features_df.loc[f'E_rest_{Inputcurr*1e12}','raw'] = np.nanmean(v[i==0])
    # AP1_amp
    features_df.loc[f'AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #APp_amp
    features_df.loc[f'APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #AP1_width
    features_df.loc[f'AP1_width_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("width")[0]
    #APp_width
    features_df.loc[f'APp_width_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("width")[-2]
    #AP1_thresh
    features_df.loc[f'AP1_thresh_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("threshold_v")[0]
    #APp_thresh
    features_df.loc[f'APp_thresh_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("threshold_v")[-2]
    #AP1_lat
    features_df.loc[f'AP1_lat_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("threshold_t")[0] - stim_start
    #ISI1
    features_df.loc[f'ISI1_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_t")[1] - sweep_ext.spike_feature("peak_t")[0]
    #ISIl
    features_df.loc[f'ISIl_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_t")[-1] - sweep_ext.spike_feature("peak_t")[-2]
    # #ISIavg
    # pt = sweep_ext.spike_feature("peak_t")
    # features_df.loc[f'ISIavg_{Inputcurr*1e12}','raw'] = np.nanmean([s-f for s,f in zip(pt[1:],pt[:-1])])
    #freq
    features_df.loc[f'freq_{Inputcurr*1e12}','raw'] = len(sweep_ext.spike_feature("peak_t"))/(stim_end - stim_start)
    #Adptn_id = 1-ISI1/ISIl
    features_df.loc[f'Adptn_id_{Inputcurr*1e12}','raw'] = 1 - features_df.loc[f'ISI1_{Inputcurr*1e12}','raw']/features_df.loc[f'ISIl_{Inputcurr*1e12}','raw']
    #fAHP_AP1_amp
    features_df.loc[f'fAHP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("fast_trough_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #fAHP_APp_amp
    features_df.loc[f'fAHP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("fast_trough_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    # #mAHP_AP1_amp
    # features_df.loc[f'mAHP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("slow_trough_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #mAHP_APp_amp
    features_df.loc[f'mAHP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("slow_trough_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    # #mAHP_AP1_dur
    # features_df.loc[f'mAHP_AP1_dur_{Inputcurr*1e12}','raw'] = (sweep_ext.spike_feature("slow_trough_t")[0] - sweep_ext.spike_feature("peak_t")[0])/features_df.loc[f'ISI1_{Inputcurr*1e12}','raw']
    #mAHP_APp_dur = mAHP of second last spike (penultimate)
    features_df.loc[f'mAHP_APp_dur_{Inputcurr*1e12}','raw'] = (sweep_ext.spike_feature("slow_trough_t")[-2] - sweep_ext.spike_feature("peak_t")[-2])/features_df.loc[f'ISIl_{Inputcurr*1e12}','raw']
    # #ADP_AP1_amp
    # features_df.loc[f'ADP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("adp_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    # #ADP_APp_amp
    # features_df.loc[f'ADP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("adp_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #mAHP_stimend_amp = within 50ms
    end50_idx = (np.abs(t - stim_end - 50e-3)).argmin()
    print(end50_idx)
    features_df.loc[f'mAHP_stimend_amp_{Inputcurr*1e12}','raw'] = np.min(v[end_idx:end50_idx]) - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #sAHP_stimend_amp = within 200ms
    end200_idx = (np.abs(t - stim_end - 200e-3)).argmin()
    features_df.loc[f'sAHP_stimend_amp_{Inputcurr*1e12}','raw'] = np.min(v[end_idx:end200_idx]) - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']

    features_df = features_df.replace('', np.nan)

    features_df['rescaled'] = ''
    for i in features_df.index:
        if features_df.loc[i,'raw'] !='' and ~np.isnan(features_df.loc[i,'raw']):
            Max = features_df.loc[i,'Max']
            Min = features_df.loc[i,'Min']
            features_df.loc[i,'rescaled'] = 2/(Max-Min)*(features_df.loc[i,'raw'] - Min) -1
    features_df = features_df.replace('', np.nan)

    features_df['cost'] = ''
    for i in features_df.index:
        if features_df.loc[i,'rescaled'] !='' and ~np.isnan(features_df.loc[i,'rescaled']):
            if features_df.loc[i,'rescaled'] >1 or features_df.loc[i,'rescaled'] <-1:
                features_df.loc[i,'cost'] = np.abs(features_df.loc[i,'rescaled'])
            else:
                features_df.loc[i,'cost'] = 0
    features_df = features_df.replace('', np.nan)

    return features_df
