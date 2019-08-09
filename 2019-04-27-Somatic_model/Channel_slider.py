#exec(open('Somatic model/Channel_slider.py').read())
# No features are calculated in this version

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

Exp_tracedir = 'Deepanjali data/WT step input cells/'
flnme = 'Cell 1 of 201217.abf'
currno = 10

def exp_tracef(currno=10):
    global flnme
    global exp_trace
    global exp_sampdur
    global exp_samprate
    global exp_samppoints
    global exp_trace_injend
    global exp_trace_injstart

    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf']
    exp_tracefile = Exp_tracedir+flnme
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
ChPo = 'Somatic model/ChannelProtos_Sri2015_base_revchanged'
ChP = 'Somatic model/ChannelProtos_Sri2015_base_revchanged'


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

Ca_B = 1000/(2*96485*0.1*18*np.pi*sm_diam*sm_len)

lines = []
tplot = []
axes_list = []
sliders_list = []
fields = []
prms = {}

diameter = sm_diam
elecPlotDt = 0.00005
preStimTime = 1
injectTime = 0.5
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

def parameters():
    soma = moose.element('/model/elec/soma')
    soma_area = np.pi*soma.diameter*soma.length
    global prms
    prms['RM'] = soma.Rm*soma_area
    prms['CM'] = soma.Cm/soma_area
    prms['RA'] = soma.Rm*np.pi*sm_diam**2/sm_len/4
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

    parameters()
    print('###############################################')

    l.set_ydata(v)
    fig.canvas.draw_idle()

exp_tracef(currno)

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
parameters()

axcolor = 'lightgoldenrodyellow'
for x in np.linspace(0.5,0.1,len(freeprm_list)):
    axes_list.append(plt.axes([0.1, x, 0.8, 0.02], facecolor=axcolor))

for i in range(len(freeprm_dict)):
    if i == 0:
        sliders_list.append(Slider(axes_list[0], 'Em', Em-0.010, Em+0.010, valinit=Em))
    else:
        sliders_list.append(Slider(axes_list[i], freeprm_list[i], freeprm_dict[freeprm_list[i]]*0.5, freeprm_dict[freeprm_list[i]]*20, valinit=freeprm_dict[freeprm_list[i]]))
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
    currno = int(4+Injectcurr/25e-12)
    print(Injectcurr)
    exp_tracef(currno)
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
currbut.on_submit(curr)

plt.show(block=False)
