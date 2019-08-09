#exec(open('Somatic model/channel_slider_myown_SK.py').read())

import io
import os
import sys
from neo.io import AxonIO
import quantities as pq
import csv
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import numpy as np
import warnings
import moose
import pickle
import rdesigneur as rd
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

####################### The experimental trace########################
def exp_tracef(currno=10):
    global flnme
    global exp_trace
    global exp_sampdur
    global exp_samprate
    global exp_samppoints
    global exp_trace_injend
    global exp_trace_injstart
    # stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf']
    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf']
    flnme = 'Cell 3 of 10717.abf'
    exp_tracefile = f'Deepanjali data/WT step input cells/{flnme}'
    reader = AxonIO(filename=exp_tracefile)
    seg = reader.read_block().segments[currno] # 10 means 15pA current
    exp_trace = seg.analogsignals[0]
    exp_samprate = float(exp_trace.sampling_rate)
    exp_sampdur = float(exp_trace.t_stop) - float(exp_trace.t_start)
    exp_samppoints = int(exp_samprate*exp_sampdur)
    if flnme in stim1391:
        exp_trace_injstart = 139.1e-3
        exp_trace_injend = 639.1e-3
    else:
        exp_trace_injstart = 81.4e-3
        exp_trace_injend = 581.4e-3
    exp_trace = np.array(exp_trace).flatten()


# print(seg) # To confirm the sampling rate and duration of 2 sec
# If your choosen file is one of stim1391, the current injection starts at 139.1e-3 and ends at 639.1e-3
# Otherwise, from 81.4e-3 till 581.4e-3

#####################################################################

################## Initial defining of parameters##################
#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

#Myown
ChPo = 'Somatic model/ChannelProtos_Sri2015_base'
ChP = 'Somatic model/ChannelProtos_Sri2015_base'
# ChP = 'Somatic model/ChannelProtos_Combe2018'
exec(open('Somatic model/feature_dict.py').read()) #Get the feature_range
F = 96485.3329
sm_diam = 100e-6
sm_len = 100e-6
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
RA = 1.0

with open('Somatic model/best_params.pkl', 'rb') as f:
    Em, RM, RA, CM, sm_diam, sm_len, \
    Ca_Basal, Ca_tau, Ca_B, Na_Gbar, K_DR_Gbar, K_A_Gbar, \
    K_M_Gbar, h_Gbar, Ca_T_Gbar, Ca_R_Gbar, \
    Ca_L_Gbar, Ca_N_Gbar, K_SK_Gbar, K_BK_Gbar = pickle.load(f)

# with open('Somatic model/cost0.pkl', 'rb') as f:
#     Em, RM, RA, CM, sm_diam, sm_len, \
#     Ca_Basal, Ca_tau, Ca_B, Na_Gbar, K_DR_Gbar, K_A_Gbar, \
#     K_M_Gbar, h_Gbar, Ca_T_Gbar, Ca_R_Gbar, \
#     Ca_L_Gbar, Ca_N_Gbar, K_SK_Gbar, K_BK_Gbar = pickle.load(f)

# #Srikanth2015 params aftertuning
# Em = -0.070
# RM = 2.8
# CM = 0.01
# Ca_Basal = 5e-5
# Ca_tau = 0.029
# Ca_B = 1/(sm_vol*F*2)
# Na_Gbar = 320
# K_DR_Gbar = 100
# K_A_Gbar = 300
# K_M_Gbar = 150
# h_Gbar =  6
# Ca_T_Gbar = 1
# Ca_R_Gbar = 30
# Ca_L_Gbar = 5
# Ca_N_Gbar = 1.76
# K_SK_Gbar = 100
# K_BK_Gbar = 100

# #Srikanth2015 params default
# Em = -0.070
# RM = 3.5
# CM = 0.01
# Ca_Basal = 5e-5
# Ca_tau = 0.029
# Ca_B = 1/(sm_vol*F*2)
# Na_Gbar = 70
# K_DR_Gbar = 30
# K_A_Gbar = 80
# K_M_Gbar = 0.01
# h_Gbar =  0.8
# Ca_T_Gbar = 1
# Ca_R_Gbar = 1
# Ca_L_Gbar = 1
# Ca_N_Gbar = 1
# K_SK_Gbar = 0.01
# K_BK_Gbar = 0.01

# # Slider params
# Em = -0.070
# RM = 2
# CM = 0.015
# Ca_Basal = 100e-6
# Ca_tau = 0.020
# Ca_B = 575792.7
# Na_Gbar = 1.2*0.035*1e4
# K_DR_Gbar = 2.2*0.015*1e4
# K_A_Gbar = 7*0.0005*1e4
# K_M_Gbar = 0.001*1e4
# h_Gbar =  1.8e-6*1e4
# Ca_T_Gbar = 0.0003*1e4
# Ca_R_Gbar = 0.00008*1e4
# Ca_L_Gbar = 0.00006*0.1*1e4
# Ca_N_Gbar = 6.415
# K_SK_Gbar = 0.7*4.5*0.0001*0.5*1e4
# K_BK_Gbar = 0.9*1.5*0.03*5.5*1e4

##########
lines = []
tplot = []
axes_list = []
sliders_list = []
fields = []
prms = {}

diameter = sm_diam
elecPlotDt = 0.00005
preStimTime = 1
injectTime = 0.1
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 150e-12
freeprm_dict = {'Em':Em, 'RM':RM, 'CM':CM, 'Ca_Basal':Ca_Basal, 'Ca_tau':Ca_tau, 'Ca_B':Ca_B, 'Na_Gbar':Na_Gbar, 'K_DR_Gbar':K_DR_Gbar, 'K_A_Gbar':K_A_Gbar, 'K_M_Gbar':K_M_Gbar, 'h_Gbar':h_Gbar, 'Ca_T_Gbar':Ca_T_Gbar, 'Ca_R_Gbar':Ca_R_Gbar, 'Ca_L_Gbar':Ca_L_Gbar, 'Ca_N_Gbar':Ca_N_Gbar, 'K_SK_Gbar':K_SK_Gbar, 'K_BK_Gbar':K_BK_Gbar}
freeprm_list = list(freeprm_dict)
##############################################################

############Building the initial model#######################
def makeModel(stim_start = preStimTime, stim_end = preStimTime+injectTime, curr = Injectcurr):
    rdes = rd.rdesigneur(
        elecPlotDt = 0.00005,
        # cellProto = [['somaProto', 'soma', 12.76e-6, 0.01e-6]],
        cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
        chanProto = [
            [ChP+'.Na_Chan()', 'Na_chan'],
            [ChP+'.K_DR_Chan()', 'K_DR_chan'],
            [ChP+'.K_A_Chan()', 'K_A_chan'],
            [ChP+'.K_M_Chan()', 'K_M_chan'],
            [ChP+'.h_Chan()', 'h_chan'],
            [ChP+'.Ca_T_Chan()', 'Ca_T_chan'],
            [ChP+'.Ca_R_Chan()', 'Ca_R_chan'],
            [ChP+'.Ca_L_Chan()', 'Ca_L_chan'],
            [ChPo+'.Ca_N_Chan()', 'Ca_N_chan'],
            [ChP+'.K_BK_Chan()', 'K_BK_chan'],
            [ChP+'.K_SK_Chan()', 'K_SK_chan'],
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

def SKproto():
    exec(open('Somatic model/myown_SKproto.py').read())
    StSI_list = []
    for Ii in Itrace_col:
        StSI_list.append(np.mean(Ii[33376:37376]))
    return StSI_list

def features(Vtrace,stim_start=preStimTime,stim_end=preStimTime+injectTime,truntime=runtime,Inputcurr=Injectcurr):
    features_df = feature_range_df.copy()

    v = np.array(Vtrace) #in mV
    t = np.linspace(0,truntime, len(Vtrace)) #in s
    start_idx = (np.abs(t - stim_start)).argmin()
    end_idx = (np.abs(t - stim_end)).argmin()
    i = np.zeros(len(t))
    i[start_idx:end_idx] = Inputcurr*1e12 #in pA

    sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, start=stim_start, end=stim_end, filter = 9.9)
    sweep_ext.process_spikes()

    features_df['raw'] = ''
    # E_rest
    if stim_start>0.5:
        features_df.loc[f'E_rest_{Inputcurr*1e12}','raw'] = np.nanmean(v[int(0.5*len(v)/truntime):int(stim_start*len(v)/truntime)])
    else:
        features_df.loc[f'E_rest_{Inputcurr*1e12}','raw'] = np.nanmean(v[-int(0.5*len(v)/truntime):])
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
    pt = sweep_ext.spike_feature("peak_t")
    features_df.loc[f'ISIavg_{Inputcurr*1e12}','raw'] = np.nanmean([s-f for s,f in zip(pt[1:],pt[:-1])])
    # features_df.loc[f'ISIavg_{Inputcurr*1e12}','raw'] = 'skip'
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
                features_df.loc[i,'cost'] = np.abs(features_df.loc[i,'rescaled'])*features_df.loc[i,'weight']
            else:
                features_df.loc[i,'cost'] = 0
    features_df = features_df.replace('', np.nan)

    return features_df

def parameters():
    soma = moose.element('/model/elec/soma')
    soma_area = np.pi*soma.diameter*soma.length
    global prms
    prms['RM'] = soma.Rm*soma_area
    prms['CM'] = soma.Cm/soma_area
    prms['RA'] = soma.Rm*soma_area/soma.length
    prms['Em'] = soma.Em
    print(f"RM={prms['RM']}; CM={prms['CM']}; RA={prms['RA']}; Em={prms['Em']}")
    for i in moose.wildcardFind( '/model/elec/soma/#[ISA=ChanBase]' ):
        prms[f'{i.name} Gbar'] = i.Gbar/soma_area
        print(f"{i.name} Gbar = {prms[f'{i.name} Gbar']}")
    Caconcc = moose.element('/model/elec/soma/Ca_conc')
    prms['Ca_Basal'] = Caconcc.CaBasal
    prms['Ca_tau'] = Caconcc.tau
    prms['Ca_B'] = Caconcc.B
    print(f"Ca_Basal = {prms['Ca_Basal']}; Ca_tau = {prms['Ca_tau']}; Ca_B = {prms['Ca_B']}")

def update(val):
    global Em
    global RM
    global CM
    global Ca_Basal
    global Ca_tau
    global Ca_B
    global Na_Gbar
    global K_DR_Gbar
    global K_A_Gbar
    global K_M_Gbar
    global h_Gbar
    global Ca_T_Gbar
    global Ca_R_Gbar
    global Ca_L_Gbar
    global Ca_N_Gbar
    global K_SK_Gbar
    global K_BK_Gbar

    Em = sliders_list[0].val
    RM = sliders_list[1].val
    CM = sliders_list[2].val
    Ca_Basal = sliders_list[3].val
    Ca_tau = sliders_list[4].val
    Ca_B = sliders_list[5].val
    Na_Gbar = sliders_list[6].val
    K_DR_Gbar = sliders_list[7].val
    K_A_Gbar = sliders_list[8].val
    K_M_Gbar = sliders_list[9].val
    h_Gbar = sliders_list[10].val
    Ca_T_Gbar = sliders_list[11].val
    Ca_R_Gbar = sliders_list[12].val
    Ca_L_Gbar = sliders_list[13].val
    Ca_N_Gbar = sliders_list[14].val
    K_SK_Gbar = sliders_list[15].val
    K_BK_Gbar = sliders_list[16].val

    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        [moose.delete(x) for x in ['/model']]
    except:
        pass
    rdes = makeModel()
    text_trap = io.StringIO() # JUst to suppress the inbuilt terminal output
    sys.stdout = text_trap # JUst to suppress the inbuilt terminal output
    rdes.buildModel()
    sys.stdout = sys.__stdout__ # restoring terminal output
    try:
        moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).td = 0
    except:
        pass
    moose.element('/model/elec/soma/Ca_conc').CaBasal = Ca_Basal
    moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
    moose.element('/model/elec/soma/Ca_conc').B = Ca_B
    moose.reinit()
    moose.start( runtime )
    Vtrace = moose.element('/model/graphs/plot0' ).vector
    v = np.array(Vtrace)

    try:
        del(features_df)
    except:
        pass
    try:
        features_df = features(v*1e3,preStimTime,preStimTime+injectTime,runtime,Injectcurr)
        features_df.insert(3,'raw_exp',features_exp['raw'])
        print(features_df)
    except:
        pass
    try:
        tcost = np.nansum(features_df['cost'])
        # print(f"Total cost = {tcost}")
    except:
        tcost = 1000
        # print(f'Sorry, No spike')
    parameters()
    prms['tcost'] = tcost

    print(f'Total cost --> {tcost}')
    print('###############################################')

    l.set_ydata(v)
    fig.canvas.draw_idle()

exp_tracef(10)
features_exp = features(exp_trace,exp_trace_injstart,exp_trace_injend,exp_sampdur,Injectcurr)

rdes = makeModel()
text_trap = io.StringIO() # JUst to suppress the inbuilt terminal output
sys.stdout = text_trap # JUst to suppress the inbuilt terminal output
rdes.buildModel()
sys.stdout = sys.__stdout__ # restoring terminal output
try:
    moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass
moose.element('/model/elec/soma/Ca_conc').CaBasal = Ca_Basal
moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
moose.element('/model/elec/soma/Ca_conc').B = Ca_B
moose.reinit()
moose.start( runtime )
Vtrace = moose.element('/model/graphs/plot0' ).vector
###################################################################

########### Initial setting up of axes############################
fig, ax = plt.subplots()
plt.subplots_adjust(top = 0.90, bottom=0.60)
t = np.linspace(0,runtime, len(Vtrace)) #in s
v = np.array(Vtrace) #in V
l, = plt.plot(t, v, lw=2, color='red',label='in-silico')
exp, = plt.plot(np.linspace(1-exp_trace_injstart,1+exp_sampdur-exp_trace_injstart,exp_samppoints), exp_trace*1e-3, label=flnme)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title(f'Injected current = {Injectcurr}A')
plt.legend()
plt.axis([0, runtime, -0.090, 0.060])

Ii = SKproto()
print(Ii - [44,116,206,322,583,1035,1574,2191,2715,3265])


parameters()

axcolor = 'lightgoldenrodyellow'
for x in np.linspace(0.5,0.1,len(freeprm_list)):
    axes_list.append(plt.axes([0.1, x, 0.8, 0.02], facecolor=axcolor))

for i in range(len(freeprm_dict)):
    if i == 0:
        sliders_list.append(Slider(axes_list[0], 'Em', Em-0.010, Em+0.010, valinit=Em))
    else:
        sliders_list.append(Slider(axes_list[i], freeprm_list[i], freeprm_dict[freeprm_list[i]]*0.05, freeprm_dict[freeprm_list[i]]*20, valinit=freeprm_dict[freeprm_list[i]]))
    sliders_list[-1].on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
resetbut = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

saveax = plt.axes([0.01, 0.025, 0.2, 0.04])
savebut = Button(saveax, 'Save params', color=axcolor, hovercolor='0.975')

currax = plt.axes([0.5, 0.025, 0.2, 0.04])
currbut = TextBox(currax, 'Input current', color=axcolor, hovercolor='0.975', initial = str(Injectcurr))

def reset(event):
    for sldr in sliders_list:
        sldr.reset()
    print('parameters reset')
resetbut.on_clicked(reset)

def save_params(event):
    with open('Somatic model/best_params.pkl', 'wb') as f:
        pickle.dump([Em, RM, RA, CM, sm_diam, sm_len, \
            Ca_Basal, Ca_tau, Ca_B, Na_Gbar, K_DR_Gbar, K_A_Gbar, \
            K_M_Gbar, h_Gbar, Ca_T_Gbar, Ca_R_Gbar, \
            Ca_L_Gbar, Ca_N_Gbar, K_SK_Gbar, K_BK_Gbar], f)
    print('Parameters saved')
savebut.on_clicked(save_params)

def curr(text):
    Injectcurr = eval(text)
    if Injectcurr==150e-12:
        currno=10
    elif Injectcurr==100e-12:
        currno=8
    elif Injectcurr==0:
        currno=4
    elif Injectcurr==300e-12:
        currno=16
    # exp_tracef(currno)
    exp.set_ydata(exp_trace*1e-3)


    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        [moose.delete(x) for x in ['/model']]
    except:
        pass
    rdes = makeModel(curr = Injectcurr)
    text_trap = io.StringIO() # JUst to suppress the inbuilt terminal output
    sys.stdout = text_trap # JUst to suppress the inbuilt terminal output
    rdes.buildModel()
    sys.stdout = sys.__stdout__ # restoring terminal output
    try:
        moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).td = 0
    except:
        pass
    moose.element('/model/elec/soma/Ca_conc').CaBasal = Ca_Basal
    moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
    moose.element('/model/elec/soma/Ca_conc').B = Ca_B
    moose.reinit()
    moose.start( runtime )
    Vtrace = moose.element('/model/graphs/plot0' ).vector
    v = np.array(Vtrace)
    l.set_ydata(v)
    fig.canvas.draw_idle()
    print(Injectcurr)
    try:
        features_df = features(v*1e3,preStimTime,preStimTime+injectTime,runtime,Injectcurr)
        features_df.insert(3,'raw_exp',features_exp['raw'])
        print(features_df)
    except:
        pass
currbut.on_submit(curr)

plt.show(block=False)