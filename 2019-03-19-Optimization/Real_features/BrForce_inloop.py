# exec(open('Optimization/Custom/Real_features/BrForce_inloop.py').read())

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import itertools
import pandas as pd
import csv
import random as rnd
import sys
import time
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor


exec(open('Optimization/Custom/Real_features/feature_dict.py').read()) #imports the feature allowed range
F = 96485.3329
ChP = 'Optimization/Custom/Real_features/ChannelProtos_Sri2015_base' #CHannel Kinetics
sm_diam = 100e-6
sm_len = 100e-6
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
RA = 1.0
stim_start = 2.0
stim_end = 2.5
runtime = 3
Inputcurr = 150e-12
elecPlotDt = 0.00005
times = 5
works = pd.DataFrame()
prms = {}

#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

#Extracting features
def features(Vtrace,stim_start,stim_end,Inputcurr):
    features_df = feature_range_df.copy()

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
    features_df.loc[f'E_rest_{Inputcurr*1e12}','raw'] = np.nanmean(v[:int(stim_start*len(v)/runtime)])
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
    features_df.loc[f'ISIavg_{Inputcurr*1e12}','raw'] = 'skip'
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
    features_df.loc[f'mAHP_AP1_amp_{Inputcurr*1e12}','raw'] = 'skip'
    #mAHP_APp_amp
    features_df.loc[f'mAHP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("slow_trough_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    # #mAHP_AP1_dur
    # features_df.loc[f'mAHP_AP1_dur_{Inputcurr*1e12}','raw'] = (sweep_ext.spike_feature("slow_trough_t")[0] - sweep_ext.spike_feature("peak_t")[0])/features_df.loc[f'ISI1_{Inputcurr*1e12}','raw']
    features_df.loc[f'mAHP_AP1_dur_{Inputcurr*1e12}','raw'] = 'skip'
    #mAHP_APp_dur = mAHP of second last spike (penultimate)
    features_df.loc[f'mAHP_APp_dur_{Inputcurr*1e12}','raw'] = (sweep_ext.spike_feature("slow_trough_t")[-2] - sweep_ext.spike_feature("peak_t")[-2])/features_df.loc[f'ISIl_{Inputcurr*1e12}','raw']
    # #ADP_AP1_amp
    # features_df.loc[f'ADP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("adp_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    features_df.loc[f'ADP_AP1_amp_{Inputcurr*1e12}','raw'] = 'skip'
    # #ADP_APp_amp
    # features_df.loc[f'ADP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("adp_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    features_df.loc[f'ADP_APp_amp_{Inputcurr*1e12}','raw'] = 'skip'
    #mAHP_stimend_amp = within 50ms
    end50_idx = (np.abs(t - stim_end - 50e-3)).argmin()
    features_df.loc[f'mAHP_stimend_amp_{Inputcurr*1e12}','raw'] = np.min(v[end_idx:end50_idx]) - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #sAHP_stimend_amp = within 200ms
    end200_idx = (np.abs(t - stim_end - 200e-3)).argmin()
    features_df.loc[f'sAHP_stimend_amp_{Inputcurr*1e12}','raw'] = np.min(v[end_idx:end200_idx]) - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']

    features_df = features_df.replace('', np.nan)

    features_df['rescaled'] = ''
    for i in features_df.index:
        if features_df.loc[i,'raw'] == 'skip':
            features_df.loc[i,'rescaled'] = 0
        elif np.isnan(features_df.loc[i,'raw']):
            features_df.loc[i,'rescaled'] = 10 #Penalty for not having the feature altogether
        else:
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

def makeModel(stim_start = 2.0, stim_end = 2.5, curr = 150e-12):
    rdes = rd.rdesigneur(
        elecPlotDt = elecPlotDt,
        stealCellFromLibrary = False,
        verbose = False,
        # cellProto = [['somaProto', 'soma', 12.76e-6, 0.01e-6]],
        cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
        chanProto = [
            [ChP+'.Na_Chan()', 'Na_chan'],
            [ChP+'.KDR_Chan()', 'K_DR_chan'],
            [ChP+'.KA_Chan()', 'K_A_chan'],
            [ChP+'.KM_Chan()', 'K_M_chan'],
            [ChP+'.h_Chan()', 'h_chan'],
            [ChP+'.CaT_Chan()', 'Ca_T_chan'],
            [ChP+'.CaR_Chan()', 'Ca_R_chan'],
            [ChP+'.CaL_Chan()', 'Ca_L_chan'],
            [ChP+'.CaN_Chan()', 'Ca_N_chan'],
            [ChP+'.KBK_Chan()', 'K_BK_chan'],
            [ChP+'.KSK_Chan()', 'K_SK_chan'],
            [ChP+'.Ca_Conc()', 'Ca_conc'],
        ],
        passiveDistrib = [
            ['soma', 'RM', str(RM), 'RA', '1.5', 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
        ],
        chanDistrib = [
            ['Na_chan', 'soma', 'Gbar', str(Na_Gbar)],
            ['K_DR_chan', 'soma', 'Gbar', str(K_DR_Gbar)],
            ['K_A_chan', 'soma', 'Gbar', str(K_A_Gbar)],
            ['K_M_chan', 'soma', 'Gbar', str(K_M_Gbar)],
            ['h_chan', 'soma', 'Gbar', str(h_Gbar)],
            ['Ca_T_chan', 'soma', 'Gbar', str(Ca_T_Gbar)],
            ['Ca_R_chan', 'soma', 'Gbar', str(Ca_R_Gbar)],
            ['Ca_L_chan', 'soma', 'Gbar', str(Ca_L_Gbar)],
            ['Ca_N_chan', 'soma', 'Gbar', str(Ca_N_Gbar)],
            ['K_SK_chan', 'soma', 'Gbar', str(K_SK_Gbar)],
            ['K_BK_chan', 'soma', 'Gbar', str(K_BK_Gbar)],
            ['Ca_conc', 'soma', 'thick', '177.9e-6'],
        ],
        stimList = [
            # ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<6.5) * 0.050' ],
            ['soma', '1', '.', 'inject', f'(t>={stim_start} && t<={stim_end}) ? {curr} : 0'],
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
            # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
            # ['soma', '1', 'Na_Schan', 'Gk', 'Soma Sodium conductance'],
            # ['soma', '1', 'KDR_Schan', 'Gk', 'Soma Kdr conductance'],
            # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
            # ['soma', '1', 'KBK_chan', 'Gk', 'Soma KBK conductance'],
            # ['soma', '1', 'KSK_chan', 'Gk', 'Soma KSK conductance'],
            # ['soma', '1', 'CaT_chan', 'Gk', 'Soma CaT conductance'],
            # ['soma', '1', 'CaR_Schan', 'Gk', 'Soma CaR_S conductance'],
            # ['soma', '1', 'CaL_Schan', 'Gk', 'Soma CaL_S conductance'],

            # ['soma', '1', 'Na_Schan', 'Ik', 'Soma Sodium current'],
            # ['soma', '1', 'KDR_Schan', 'Ik', 'Soma Kdr current'],
            # ['soma', '1', 'KA_Schan', 'Ik', 'Soma KA current'],
            # ['soma', '1', 'KM_chan', 'Ik', 'Soma KM current'],
            # ['soma', '1', 'h_chan', 'Ik', 'Soma h current'],
            # ['soma', '1', 'KSK_chan', 'Ik', 'Soma KSK current'],
            # ['soma', '1', 'KBK_chan', 'Ik', 'Soma KBK current'],
            # ['soma', '1', 'CaT_chan', 'Ik', 'Soma CaT current'],
            # ['soma', '1', 'CaR_Schan', 'Ik', 'Soma CaR_S current'],
            # ['soma', '1', 'CaL_Schan', 'Ik', 'Soma CaL_S current'],
            # ['soma', '1', 'CaN_Schan', 'Ik', 'Soma CaN_S current'],
        ],
    )
    return rdes

def parameters():
    soma = moose.element('/model/elec/soma')
    soma_area = np.pi*soma.diameter*soma.length
    global prms
    prms['RM'] = soma.Rm*soma_area
    prms['CM'] = soma.Cm/soma_area
    prms['RA'] = soma.Cm*soma_area/soma.length
    prms['Em'] = soma.Em
    print(f"RM={prms['RM']}; CM={prms['CM']}; RA={prms['RA']}; Em={prms['Em']}")
    for i in moose.wildcardFind( '/model/elec/soma/#[ISA=ChanBase]' ):
        prms[f'{i.name} Gbar'] = i.Gbar/soma_area
        print(f"{i.name} Gbar = {prms[f'{i.name} Gbar']}")
    Caconcc = moose.element('/model/elec/soma/Ca_conc')
    prms['Ca_Basal'] = Caconcc.CaBasal
    prms['Ca_B'] = Caconcc.B
    prms['Ca_tau'] = Caconcc.tau
    print(f"Ca_Basal = {prms['Ca_Basal']}; Ca_B = {prms['Ca_B']}; Ca_tau = {prms['Ca_tau']}")

st = time.time()
for i in range(times):
    try:
        del(features_df)
    except:
        pass
    exec(open('Optimization/Custom/Real_features/BrForce.py').read())
    if tcost<1001:
        parameters()
        prms['tcost'] = tcost
        works = works.append(prms,ignore_index=True)
        works.to_csv('Optimization/Custom/Real_features/works.csv')
    print(f'{i} --> {tcost}')
    plt.plot(Vtrace)
    plt.show()
tt = time.time() - st
print(tt)
