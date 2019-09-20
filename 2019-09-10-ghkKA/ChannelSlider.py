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

foldername=os.path.basename(os.getcwd())
#Parameters on sliders: Em, Gl, Ca_concCaBasal,

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
Pr={'K_D_chan':98.6e-2, 'h_chan':20.2e-2, 'K_M_chan':4.74e-2} #Default at -70mV
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

####################################################################################
# def set_initialgatevalues(Vrest=-0.07, Carest=0.13e-3):
#     for chan in moose.wildcardFind('/library/#[TYPE=HHChannel]'):
#         if chan.Xpower>0:
#             gate = moose.element( chan.path + '/gateX' )
#             v = np.linspace(gate.min,gate.max,gate.divs)
#             idx = np.argmin(abs(v-Vrest))
#             chan.X = gate.tableA[idx]/gate.tableB[idx]
#             Pr[chan.name] = chan.X #Since we don't care about channels other than h, KM, and KD
#         if chan.Ypower>0:
#             gate = moose.element( chan.path + '/gateY' )
#             v = np.linspace(gate.min,gate.max,gate.divs)
#             idx = np.argmin(abs(v-Vrest))
#             chan.Y = gate.tableA[idx]/gate.tableB[idx]
#         if chan.Zpower>0:
#             gate = moose.element( chan.path + '/gateZ' )
#             vorca = np.linspace(gate.min,gate.max,gate.divs)
#             if chan.useConcentration == 1:
#                 idx = np.argmin(abs(vorca-Carest))
#             else:
#                 idx = np.argmin(abs(vorca-Vrest))
#             chan.Z = gate.tableA[idx]/gate.tableB[idx]
#     for chan in moose.wildcardFind('/library/#[TYPE=HHChannel2D]'):
#         if chan.Xpower>0:
#             gate = moose.element( chan.path + '/gateX' )
#             v = np.linspace(gate.xminA,gate.xmaxA,gate.xdivsA)
#             ca = np.linspace(gate.yminA,gate.ymaxA,gate.ydivsA)
#             idxv = np.argmin(abs(v-Vrest))
#             if len(ca)==0:
#                 idxca = 0
#             else:
#                 idxca = np.argmin(abs(ca-Carest))
#             chan.X = gate.tableA[idx][idxca]/gate.tableB[idx][idxca]
#         if chan.Ypower>0:
#             gate = moose.element( chan.path + '/gateY' )
#             v = np.linspace(gate.xminA,gate.xmaxA,gate.xdivsA)
#             ca = np.linspace(gate.yminA,gate.ymaxA,gate.ydivsA)
#             idxv = np.argmin(abs(v-Vrest))
#             if len(ca)==0:
#                 idxca = 0
#             else:
#                 idxca = np.argmin(abs(ca-Carest))
#             chan.Y = gate.tableA[idx][idxca]/gate.tableB[idx][idxca]
#         if chan.Zpower>0:
#             gate = moose.element( chan.path + '/gateZ' )
#             v = np.linspace(gate.xminA,gate.xmaxA,gate.xdivsA)
#             ca = np.linspace(gate.yminA,gate.ymaxA,gate.ydivsA)
#             idxv = np.argmin(abs(v-Vrest))
#             if len(ca)==0:
#                 idxca = 0
#             else:
#                 idxca = np.argmin(abs(ca-Carest))
#             chan.Z = gate.tableA[idx][idxca]/gate.tableB[idx][idxca]
#     moose.element('/library/Ca_conc').Ca = Carest

def calcPr(Vrest=-0.07, Carest=0.13e-3):
    for chan in moose.wildcardFind('/library/#[TYPE=HHChannel]'):
        if chan.Xpower>0:
            gate = moose.element( chan.path + '/gateX' )
            v = np.linspace(gate.min,gate.max,gate.divs)
            idx = np.argmin(abs(v-Vrest))
            Pr[chan.name] = gate.tableA[idx]/gate.tableB[idx] #Since we don't care about channels other than h, KM, and KD

# def set_initialGk(Vrest=-0.07, Carest=0.13e-3):
#     for chan in moose.wildcardFind('/model/elec/soma/#[TYPE=HHChannel]'):
#         chan.Gk = chan.Gbar*chan.X**chan.Xpower*chan.Y**chan.Ypower*chan.Z**chan.Zpower
#     for chan in moose.wildcardFind('/model/elec/soma/#[TYPE=HHChannel2D]'):
#         chan.Gk = chan.Gbar*chan.X**chan.Xpower*chan.Y**chan.Ypower*chan.Z**chan.Zpower

exec(open('Modelfunc.py').read())
# tinitial = time.time()
if os.path.exists(f'../../Output/{foldername}/Parametersdf.csv'):
    prevrun = 1
    invidx = 4
    Pl = pd.read_csv(f'../../Output/{foldername}/Parametersdf.csv').tail(invidx).iloc[0]
    Parameters = {key:Pl[key] for key in Pl.keys()}
else:
    prevrun = 0
    Parameters = {}
    Parameters['sm_diam'] = sm_diam
    Parameters['sm_len'] = sm_len
    Parameters['RM'] = 2.053773343
    Parameters['CM'] = CM
    Parameters['Em'] = Vrest
    Parameters['Ca_concCaBasal'] = 0.05e-3
    Parameters['Ca_conctau'] = 0.1
    Parameters['Ca_concB'] = 2291006575
    Parameters['Ca_L_changbar'] = 5.051749119
    Parameters['Ca_N_changbar'] = 5.533811207
    Parameters['Ca_T_changbar'] = 2.910963364
    Parameters['h_changbar'] = 0.771319723
    Parameters['K_A_changbar'] = 0.435938559
    Parameters['K_BK_changbar'] = 0.010677492
    Parameters['K_D_changbar'] = 0.000219033
    Parameters['K_DR_changbar'] = 4.477080862
    Parameters['K_M_changbar'] = 0.02933468
    Parameters['K_SK_changbar'] = 0.025938828
    Parameters['Na_changbar'] = 51.0502324

[Parameters, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=runtime, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'], Ca_concCaBasal=str(Parameters['Ca_concCaBasal']), Ca_conctau=str(Parameters['Ca_conctau']), Ca_L_changbar=str(Parameters['Ca_L_changbar']), Ca_N_changbar=str(Parameters['Ca_N_changbar']), Ca_T_changbar=str(Parameters['Ca_T_changbar']), h_changbar=str(Parameters['h_changbar']), K_A_changbar=str(Parameters['K_A_changbar']), K_BK_changbar=str(Parameters['K_BK_changbar']), K_D_changbar=str(Parameters['K_D_changbar']), K_DR_changbar=str(Parameters['K_DR_changbar']), K_M_changbar=str(Parameters['K_M_changbar']), K_SK_changbar=str(Parameters['K_SK_changbar']), Na_changbar=str(Parameters['Na_changbar']))
print(Parameters)
# set_initialgatevalues(Vrest=Vrest, Carest=0.13e-3)
# calcPr(Vrest=Vrest, Carest=0.13e-3)
########### Defining some Matplotlib functions####################
def donothing(val):
    pass
def update(val):
    tme=time.time()
    global Parameters
    global sliders_list
    global characteristics
    global Vmvec
    global Ivec
    global Cavec
    global tvec
    Parameters['sm_diam'] = sm_diam
    Parameters['sm_len'] = sm_len
    Parameters['RM'] = 1/sliders_list[1].val
    Parameters['CM'] = CM
    Parameters['Em'] = sliders_list[0].val
    Parameters['Ca_concCaBasal'] = sliders_list[2].val
    Parameters['Ca_conctau'] = sliders_list[3].val
    Parameters['Ca_concB'] = sliders_list[4].val
    Parameters['Ca_L_changbar'] = sliders_list[5].val
    Parameters['Ca_N_changbar'] = sliders_list[6].val
    Parameters['Ca_T_changbar'] = sliders_list[7].val
    Parameters['h_changbar'] = sliders_list[8].val
    Parameters['K_A_changbar'] = sliders_list[9].val
    Parameters['K_BK_changbar'] = sliders_list[10].val
    Parameters['K_D_changbar'] = sliders_list[11].val
    Parameters['K_DR_changbar'] = sliders_list[12].val
    Parameters['K_M_changbar'] = sliders_list[13].val
    Parameters['K_SK_changbar'] = sliders_list[14].val
    Parameters['Na_changbar'] = sliders_list[15].val

    print('#####################################')
    print(Parameters)
    [Parameters, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=runtime, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'], Ca_concCaBasal=str(Parameters['Ca_concCaBasal']), Ca_L_changbar=str(Parameters['Ca_L_changbar']), Ca_N_changbar=str(Parameters['Ca_N_changbar']), Ca_T_changbar=str(Parameters['Ca_T_changbar']), h_changbar=str(Parameters['h_changbar']), K_A_changbar=str(Parameters['K_A_changbar']), K_BK_changbar=str(Parameters['K_BK_changbar']), K_D_changbar=str(Parameters['K_D_changbar']), K_DR_changbar=str(Parameters['K_DR_changbar']), K_M_changbar=str(Parameters['K_M_changbar']), K_SK_changbar=str(Parameters['K_SK_changbar']), Na_changbar=str(Parameters['Na_changbar']))

    for s in sliders_list:
        s.observers[0] = donothing
    sliders_list[0].set_val(Parameters['Em'])
    sliders_list[1].set_val(1/Parameters['RM'])
    sliders_list[8].set_val(Parameters['h_changbar'])
    sliders_list[10].set_val(Parameters['K_BK_changbar'])
    sliders_list[11].set_val(Parameters['K_D_changbar'])
    sliders_list[13].set_val(Parameters['K_M_changbar'])
    for s in sliders_list:
        s.observers[0] = update

    l.set_ydata(Vmvec)
    fig.canvas.draw_idle()
    print(time.time()-tme)

########### Initial setting up of axes############################
fig, ax = plt.subplots()
plt.subplots_adjust(top = 0.90, bottom=0.60)
l, = plt.plot(tvec, Vmvec, lw=2, color='red',label='in-silico')
exp, = plt.plot(np.linspace(preStimTime-exp_trace_injstart,preStimTime+exp_sampdur-exp_trace_injstart,exp_samppoints), exp_trace*1e-3, label=flnme)
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title(f'Injected current = {Injectcurr}A')
plt.legend()
plt.axis([0, runtime, -0.100, 0.060])

axcolor = 'lightgoldenrodyellow'
for x in np.linspace(0.5,0.1,16):
    axes_list.append(plt.axes([0.1, x, 0.8, 0.02], facecolor=axcolor))

sliders_list.append(Slider(axes_list[0], 'Em', Parameters['Em']-0.010, Parameters['Em']+0.010, valinit=Parameters['Em'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[1], 'gl', mingl, Gin/sm_area, valinit=1/Parameters['RM'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[2], 'Ca_concCaBasal', Parameters['Ca_concCaBasal']*0.1, Parameters['Ca_concCaBasal']*10, valinit=Parameters['Ca_concCaBasal'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[3], 'Ca_conctau', Parameters['Ca_conctau']*0.1, Parameters['Ca_conctau']*10, valinit=Parameters['Ca_conctau'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[4], 'Ca_concB', Parameters['Ca_concB']*0.1, Parameters['Ca_concB']*10, valinit=Parameters['Ca_concB'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[5], 'Ca_L_changbar', Parameters['Ca_L_changbar']*0.1, Parameters['Ca_L_changbar']*10, valinit=Parameters['Ca_L_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[6], 'Ca_N_changbar', Parameters['Ca_N_changbar']*0.1, Parameters['Ca_N_changbar']*10, valinit=Parameters['Ca_N_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[7], 'Ca_T_changbar', Parameters['Ca_T_changbar']*0.1, Parameters['Ca_T_changbar']*10, valinit=Parameters['Ca_T_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[8], 'h_changbar', Parameters['h_changbar']*0.1, Parameters['h_changbar']*10, valinit=Parameters['h_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[9], 'K_A_changbar', Parameters['K_A_changbar']*0.1, Parameters['K_A_changbar']*10, valinit=Parameters['K_A_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[10], 'K_BK_changbar', Parameters['K_BK_changbar']*0.1, Parameters['K_BK_changbar']*10, valinit=Parameters['K_BK_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[11], 'K_D_changbar', Parameters['K_D_changbar']*0.1, Parameters['K_D_changbar']*10, valinit=Parameters['K_D_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[12], 'K_DR_changbar', Parameters['K_DR_changbar']*0.1, Parameters['K_DR_changbar']*10, valinit=Parameters['K_DR_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[13], 'K_M_changbar', Parameters['K_M_changbar']*0.1, Parameters['K_M_changbar']*10, valinit=Parameters['K_M_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[14], 'K_SK_changbar', Parameters['K_SK_changbar']*0.1, Parameters['K_SK_changbar']*10, valinit=Parameters['K_SK_changbar'], valfmt="%1.2e"))
sliders_list.append(Slider(axes_list[15], 'Na_changbar', Parameters['Na_changbar']*0.1, Parameters['Na_changbar']*10, valinit=Parameters['Na_changbar'], valfmt="%1.2e"))
for s in sliders_list:
    cid.append(s.on_changed(update))

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
    l.set_ydata(Vmvec)
    fig.canvas.draw_idle()
resetbut.on_clicked(reset)

def save_params(event):
    Parametersdf = pd.DataFrame([Parameters])
    characteristicsdf = pd.DataFrame([characteristics])
    Vmvecdf = pd.DataFrame([Vmvec])
    Ivecdf = pd.DataFrame([Ivec])
    Cavecdf = pd.DataFrame([Cavec])
    tvecdf = pd.DataFrame([tvec])
    if os.path.exists(f'../../Output/{foldername}/Parametersdf.csv'):
        Parametersdf.to_csv(f'../../Output/{foldername}/Parametersdf.csv', index=None, mode='a', header=False)
    else:
        Parametersdf.to_csv(f'../../Output/{foldername}/Parametersdf.csv', index=None, header=True)
    if os.path.exists(f'../../Output/{foldername}/characteristicsdf.csv'):
        characteristicsdf.to_csv(f'../../Output/{foldername}/characteristicsdf.csv', index=None, mode='a', header=False)
    else:
        characteristicsdf.to_csv(f'../../Output/{foldername}/characteristicsdf.csv', index=None, header=True)
    if os.path.exists(f'../../Output/{foldername}/Vmvecdf.csv'):
        Vmvecdf.to_csv(f'../../Output/{foldername}/Vmvecdf.csv', index=None, mode='a', header=False)
    else:
        Vmvecdf.to_csv(f'../../Output/{foldername}/Vmvecdf.csv', index=None, header=True)
    if os.path.exists(f'../../Output/{foldername}/Ivecdf.csv'):
        Ivecdf.to_csv(f'../../Output/{foldername}/Ivecdf.csv', index=None, mode='a', header=False)
    else:
        Ivecdf.to_csv(f'../../Output/{foldername}/Ivecdf.csv', index=None, header=True)
    if os.path.exists(f'../../Output/{foldername}/Cavecdf.csv'):
        Cavecdf.to_csv(f'../../Output/{foldername}/Cavecdf.csv', index=None, mode='a', header=False)
    else:
        Cavecdf.to_csv(f'../../Output/{foldername}/Cavecdf.csv', index=None, header=True)
    if os.path.exists(f'../../Output/{foldername}/tvecdf.csv'):
        tvecdf.to_csv(f'../../Output/{foldername}/tvecdf.csv', index=None, mode='a', header=False)
    else:
        tvecdf.to_csv(f'../../Output/{foldername}/tvecdf.csv', index=None, header=True)
    print('Parameters saved')
savebut.on_clicked(save_params)

def curr(text):
    global Injectcurr
    global Parameters
    global characteristics
    global Vmvec
    global Ivec
    global Cavec
    global tvec
    Injectcurr = eval(text)
    [Parameters, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=runtime, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'], Ca_concCaBasal=str(Parameters['Ca_concCaBasal']), Ca_L_changbar=str(Parameters['Ca_L_changbar']), Ca_N_changbar=str(Parameters['Ca_N_changbar']), Ca_T_changbar=str(Parameters['Ca_T_changbar']), h_changbar=str(Parameters['h_changbar']), K_A_changbar=str(Parameters['K_A_changbar']), K_BK_changbar=str(Parameters['K_BK_changbar']), K_D_changbar=str(Parameters['K_D_changbar']), K_DR_changbar=str(Parameters['K_DR_changbar']), K_M_changbar=str(Parameters['K_M_changbar']), K_SK_changbar=str(Parameters['K_SK_changbar']), Na_changbar=str(Parameters['Na_changbar']))
    exp_trace = exp_tracef(Injectcurr=Injectcurr)
    print(Injectcurr)
    l.set_ydata(Vmvec)
    exp.set_ydata(exp_trace*1e-3)
    fig.canvas.draw_idle()
currbut.on_submit(curr)

plt.show(block=False)
