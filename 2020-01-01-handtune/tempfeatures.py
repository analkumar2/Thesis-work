## exec(open('tempfeatures.py').read())

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
import csv
import scipy.signal as scs
import pandas as pd
import brute_curvefit
from pprint import pprint
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor


def modelfeatures(stim_start=0.5,stim_end=1):
    Vtrace0 = np.load('rep1_0pA.npy')
    tt = np.linspace(0,1.2, len(Vtrace0))
    t0, Vtrace0 = tt, np.load('rep1_0pA.npy')
    t25, Vtrace25 = tt, np.load('rep1_25pA.npy')
    t150, Vtrace150 = tt, np.load('rep1_150pA.npy')
    t300, Vtrace300 = tt, np.load('rep1_300pA.npy')

    features = {}
    # features['Model name'] = cellpath.split('/')[-1]
    features['Sampling rate'] = len(t25)/t25[-1]
    features['stim_start'] = stim_start
    features['stim_end'] = stim_end

    ## Erest
    features[f'E_rest_0'] = np.median(Vtrace0)
    features[f'E_rest_25'] = np.median(Vtrace25[t25<stim_start])
    features[f'E_rest_150'] = np.median(Vtrace150[t150<stim_start])
    features[f'E_rest_300'] = np.median(Vtrace300[t300<stim_start])

    ## Input resistance and Cell capacitance
    def charging25(t, R,C):
        return features[f'E_rest_25'] + R*25e-12*(1-np.exp(-t/R/C))

    tempv = Vtrace25[(t25>stim_start) & (t25<stim_start+0.1)]
    RCfitted_ch25, error25 = brute_curvefit.brute_scifit(charging25, np.linspace(0,0.1,len(tempv)), tempv, restrict=[[50e6,50e-12],[250e6,250e-6]], ntol = 10000, printerrors=False)
    features[f'Input resistance'], features[f'Cell capacitance'] = RCfitted_ch25

    for I in [150e-12, 300e-12]:
        if I==150e-12:
            tt = t150
            vv = Vtrace150
            Erest = features[f'E_rest_150']
        elif I==300e-12:
            tt = t300
            vv = Vtrace300
            Erest = features[f'E_rest_300']
        ii = np.zeros(len(tt))
        ii[(tt>=stim_start) & (tt<=stim_end)] = I
        sweep_ext = EphysSweepFeatureExtractor(t=tt, v=vv*1e3, i=ii*1e12, filter=len(tt)/tt[-1]/2500, min_height=20)
        sweep_ext.process_spikes()
        if len(sweep_ext.spike_feature("width"))<=3:
            return False

        # AP1_amp
        features[f'AP1_amp_{I}'] = sweep_ext.spike_feature("peak_v")[0]*1e-3 - Erest

        #APp_amp
        features[f'APp_amp_{I}'] = sweep_ext.spike_feature("peak_v")[-2]*1e-3 - Erest

        #AP1_width
        features[f'AP1_width_{I}'] = sweep_ext.spike_feature("width")[0]

        #APp_width
        features[f'APp_width_{I}'] = sweep_ext.spike_feature("width")[-2]

        #AP1_thresh
        features[f'AP1_thresh_{I}'] = sweep_ext.spike_feature("threshold_v")[0]*1e-3

        #APp_thresh
        features[f'APp_thresh_{I}'] = sweep_ext.spike_feature("threshold_v")[-2]*1e-3

        #AP1_lat
        features[f'AP1_lat_{I}'] = sweep_ext.spike_feature("threshold_t")[0] - stim_start

        #ISI1
        features[f'ISI1_{I}'] = sweep_ext.spike_feature("peak_t")[1] - sweep_ext.spike_feature("peak_t")[0]

        #ISIl
        features[f'ISIl_{I}'] = sweep_ext.spike_feature("peak_t")[-1] - sweep_ext.spike_feature("peak_t")[-2]

        #ISIavg
        pt = sweep_ext.spike_feature("peak_t")
        features[f'ISIavg_{I}'] = np.nanmean([s-f for s,f in zip(pt[1:],pt[:-1])])

        #freq
        features[f'freq_{I}'] = len(sweep_ext.spike_feature("peak_t"))/(stim_end - stim_start)

        #Adptn_id = 1-ISI1/ISIl
        features[f'Adptn_id_{I}'] = 1 - features[f'ISI1_{I}']/features[f'ISIl_{I}']

        #fAHP_AP1_amp
        features[f'fAHP_AP1_amp_{I}'] = sweep_ext.spike_feature("fast_trough_v")[0]*1e-3 - Erest

        #fAHP_APp_amp
        features[f'fAHP_APp_amp_{I}'] = sweep_ext.spike_feature("fast_trough_v")[-2]*1e-3 - Erest

        #mAHP_AP1_amp
        features[f'mAHP_AP1_amp_{I}'] = sweep_ext.spike_feature("slow_trough_v")[0]*1e-3 - Erest

        #mAHP_APp_amp
        features[f'mAHP_APp_amp_{I}'] = sweep_ext.spike_feature("slow_trough_v")[-2]*1e-3 - Erest

        #mAHP_AP1_time
        features[f'mAHP_AP1_time_{I}'] = sweep_ext.spike_feature("slow_trough_t")[0] - sweep_ext.spike_feature("peak_t")[0]

        #mAHP_APp_time = mAHP of second last spike (penultimate)
        features[f'mAHP_APp_time_{I}'] = sweep_ext.spike_feature("slow_trough_t")[-2] - sweep_ext.spike_feature("peak_t")[-2]

        #ADP_AP1_amp
        features[f'ADP_AP1_amp_{I}'] = sweep_ext.spike_feature("adp_v")[0]*1e-3 - Erest

        #ADP_APp_amp
        features[f'ADP_APp_amp_{I}'] = sweep_ext.spike_feature("adp_v")[-2]*1e-3 - Erest

        #mAHP_stimend_amp = within 50ms
        features[f'mAHP_stimend_amp_{I}'] = np.min(vv[(tt>stim_end)&(tt<stim_end+0.050)]) - Erest

        #sAHP_stimend_amp = within 200ms
        features[f'sAHP_stimend_amp_{I}'] = np.min(vv[(tt>stim_end)&(tt<stim_end+0.200)]) - Erest

        #AHP_AP1_amp = Minimum between the first two spikes
        features[f'AHP_AP1_amp_{I}'] = np.nanmin(vv[(tt>sweep_ext.spike_feature("peak_t")[0]) & (tt<sweep_ext.spike_feature("peak_t")[1])]) - Erest

        #AHP_APp_amp = Minimum between the last two spikes
        features[f'AHP_APp_amp_{I}'] = np.nanmin(vv[(tt>sweep_ext.spike_feature("peak_t")[-2]) & (tt<sweep_ext.spike_feature("peak_t")[-1])]) - Erest

        #AHP_AP1_time = Minimum between the first two spikes
        features[f'AHP_AP1_time_{I}'] = np.argmin(vv[(tt>=sweep_ext.spike_feature("peak_t")[0]) & (tt<sweep_ext.spike_feature("peak_t")[1])])/features['Sampling rate']

        #AHP_APp_time = Minimum between the last two spikes
        features[f'AHP_APp_time_{I}'] = np.argmin(vv[(tt>=sweep_ext.spike_feature("peak_t")[-2]) & (tt<sweep_ext.spike_feature("peak_t")[-1])])/features['Sampling rate']

    return features

def modelscore(modeldict=None, meanfeature=meanfeatures, stdfeature=stdfeatures, modelfeature=modelfeatures(stim_start=0.5,stim_end=1)):
    if modelfeature==None:
        modelfeature = modelfeatures(modeldict,stim_start=1,stim_end=1.5)
    if modelfeature==False:
        return False
    def removestrcol(dictt):
        dicttk = dictt.keys()
        for keyy in dicttk:
            if isinstance(dictt[keyy],str):
                dictt.pop(keyy)
    [removestrcol(dictt) for dictt in [meanfeature, stdfeature, modelfeature]]

    score = {}
    for keyy in meanfeature.keys():
        score[keyy] = (modelfeature[keyy]-meanfeature[keyy])/stdfeature[keyy]

    return score


score = modelscore()
scoreabove = {}
for keyy in score.keys():
    if score[keyy]>2 or score[keyy]<-2:
        scoreabove[keyy] = score[keyy]

print(scoreabove)