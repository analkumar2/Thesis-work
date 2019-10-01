#exec(open('ChannelSlider.py').read())

import os
import sys
from neo.io import AxonIO
import csv
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import numpy as np
import warnings
import moose
import pickle
import rdesigneur as rd
import time
import quantities as pq
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

foldername=os.path.basename(os.getcwd()) #Where parameter values will be stored

# Define constants not to be changed
F = 96485.3329
elecPlotDt = 0.00005

# Define model parameters
sm_diam=60e-6; sm_len=60e-6 #To be set so that Cm is as needed when CM=0.01
CM = 0.01
depth = 0.1 # No units. For ca_conc B

# Define characteristics we need in the model or is used to better the model
Gin=7.35e-9
mingl = 0.15
Pr={'K_D_chan':98.6e-2, 'h_chan':20.2e-2, 'K_M_chan':4.74e-2, 'Na_P_chan': 0.128} #Default at -70mV
Vrest = -0.065

# Calculate secondary parameters
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
maxgl = Gin/sm_area

# Stimulus
preStimTime = 0.5
injectTime = 0.5
postStimTime = 0.2
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 150e-12
Vleveli = Vrest
Vlevelf = -0.0

# Plotting stuff
lines = []
tplot = []
axes_list = []
sliders_list = []
cid = []

################Read experiment file#######################
def exp_tracef(Injectcurr=150e-12):
    global flnme
    global exp_sampdur
    global exp_samprate
    global exp_samppoints
    global exp_trace_injend
    global exp_trace_injstart
    global Vrest
    global sm_diam
    global sm_len
    global Gin
    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf']
    # flnme = 'Cell 3 of 10717.abf'
    flnme = 'cell 4 of 61016.abf'
    exp_tracefile = f'../../Raw_data/Deepanjali_data/WT step input cells/{flnme}'
    reader = AxonIO(filename=exp_tracefile)
    currno = int(Injectcurr*1e12/25+4)
    seg = reader.read_block().segments[currno] # 10 means 150pA current
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

    exp_trace = np.array(exp_trace).flatten()*1e-3
    Vrest = np.mean(np.array(reader.read_block().segments[4].analogsignals[0]).flatten())*1e-3
    Rinp25 = np.abs(np.max(np.array(reader.read_block().segments[5].analogsignals[0]).flatten())*1e-3 - Vrest)/25e-12
    Rinn25 = np.abs(np.min(np.array(reader.read_block().segments[3].analogsignals[0]).flatten())*1e-3 - Vrest)/25e-12
    Gin = 2/(Rinp25+Rinn25)
    print(Gin)

    t = np.linspace(0,exp_sampdur, int(exp_sampdur*exp_samprate))
    Vtracep25 = np.array(reader.read_block().segments[5].analogsignals[0]).flatten()*1e-3
    Vtracen25 = np.array(reader.read_block().segments[3].analogsignals[0]).flatten()*1e-3
    str_ind = (np.abs(t-exp_trace_injend+100e-3)).argmin()
    stp_ind = (np.abs(t-exp_trace_injend)).argmin()
    Vtracep25_choppped = Vtracep25[:stp_ind]
    Vtracen25_choppped = Vtracen25[:stp_ind]
    vp63 = np.abs(np.max(np.array(reader.read_block().segments[5].analogsignals[0]).flatten())*1e-3 - Vrest)*0.63 + Vrest
    vn63 = -np.abs(np.min(np.array(reader.read_block().segments[3].analogsignals[0]).flatten())*1e-3 - Vrest)*0.63 + Vrest
    tau = (t[(np.abs(Vtracep25_choppped-vp63)).argmin()] - exp_trace_injstart + t[(np.abs(Vtracen25_choppped-vn63)).argmin()] - exp_trace_injstart)/2
    tauinv = (t[len(Vtracep25_choppped)-(np.abs(Vtracep25_choppped[stp_ind::-1]-vp63)).argmin()] - exp_trace_injstart + t[len(Vtracep25_choppped)-(np.abs(Vtracep25_choppped[::-1]-vp63)).argmin()] - exp_trace_injstart)/2
    Cm = (tau+tauinv)*Gin/2
    print(Cm)
    sm_len = np.sqrt(Cm/CM/np.pi)
    sm_diam = sm_len

    return exp_trace

exp_trace150 = exp_tracef(Injectcurr=150e-12)
exp_trace300 = exp_tracef(Injectcurr=300e-12)
Vleveli = Vrest #Reassigning as Vrest was changed in the exp_trace function


def calcPr(Vrest=-0.07, Carest=0.13e-3): #Calculates the probability of the ion cchannels open at a particular V and Ca (Currently only for KD, h, KM, NaP)
    for chan in moose.wildcardFind('/library/#[TYPE=HHChannel]'):
        if chan.Xpower>0:
            gate = moose.element( chan.path + '/gateX' )
            v = np.linspace(gate.min,gate.max,gate.divs)
            idx = np.argmin(abs(v-Vrest))
            Pr[chan.name] = gate.tableA[idx]/gate.tableB[idx] #Since we don't care about channels other than h, NaP, KM, and KD

exec(open('Modelfunc.py').read())

#Random parameter values
Parameters = {}
Parameters['sm_diam'] = sm_diam
Parameters['sm_len'] = sm_len
Parameters['RM'] = 2.053773343
Parameters['CM'] = CM
Parameters['Em'] = Vrest
Parameters['Ca_concCaBasal'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.05e-3
Parameters['Ca_conctau'] = (0.2 - 0.01)*np.random.random(1)[0] + 0.01 #0.1
Parameters['Ca_concB'] = (5e9 - 1e9)*np.random.random(1)[0] + 1e9 #2291006575
Parameters['Ca_L_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #5.051749119
Parameters['Ca_N_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #5.533811207
Parameters['Ca_T_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #2.910963364
Parameters['h_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.771319723
Parameters['K_A_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.435938559
Parameters['K_BK_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.010677492
Parameters['K_D_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.000219033
Parameters['K_DR_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #4.477080862
Parameters['K_M_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.02933468
Parameters['K_SK_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.025938828
Parameters['Na_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #51.0502324
Parameters['Na_P_changbar'] = (0.1e-3 - 0.01e-3)*np.random.random(1)[0] + 0.01e-3 #0.510502324e-1

[Parameters, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=runtime, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'], Ca_concCaBasal=str(Parameters['Ca_concCaBasal']), Ca_conctau=str(Parameters['Ca_conctau']), Ca_L_changbar=str(Parameters['Ca_L_changbar']), Ca_N_changbar=str(Parameters['Ca_N_changbar']), Ca_T_changbar=str(Parameters['Ca_T_changbar']), h_changbar=str(Parameters['h_changbar']), K_A_changbar=str(Parameters['K_A_changbar']), K_BK_changbar=str(Parameters['K_BK_changbar']), K_D_changbar=str(Parameters['K_D_changbar']), K_DR_changbar=str(Parameters['K_DR_changbar']), K_M_changbar=str(Parameters['K_M_changbar']), K_SK_changbar=str(Parameters['K_SK_changbar']), Na_changbar=str(Parameters['Na_changbar']), Na_P_changbar=str(Parameters['Na_P_changbar']))
print(Parameters)
