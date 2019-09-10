#exec(open('plotcsv.py').read())

import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

invidx = 4

foldername = os.path.basename(os.getcwd())
Pl = pd.read_csv(f'../../Output/{foldername}/Parametersdf.csv').tail(invidx).iloc[0]
Parameters = {key:Pl[key] for key in Pl.keys()}

Vtrace = list(pd.read_csv(f'../../Output/{foldername}/Vmvecdf.csv').tail(invidx).iloc[0])
ttrace = list(pd.read_csv(f'../../Output/{foldername}/tvecdf.csv').tail(invidx).iloc[0])

def exp_tracef(Injectcurr=150e-12):
    global flnme
    global exp_sampdur
    global exp_samprate
    global exp_samppoints
    global exp_trace_injend
    global exp_trace_injstart
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
    exp_trace = np.array(exp_trace).flatten()
    return exp_trace
exp_trace = exp_tracef(Injectcurr=Injectcurr)


plt.plot(np.linspace(preStimTime-exp_trace_injstart,preStimTime+exp_sampdur-exp_trace_injstart,exp_samppoints), exp_trace*1e-3, label=flnme)
plt.plot(ttrace,Vtrace, label='KA as ghk')
plt.axis([0, 1.2, -0.100, 0.060])
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (mV)')
print(Parameters)
plt.show()
