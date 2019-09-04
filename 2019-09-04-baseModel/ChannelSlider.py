#exec(open('ChannelSlider.py').read())

import os
import sys
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
P_K_D_chan_m70 = 98.6e-2
P_h_chan_m70 = 20.2e-2
P_K_M_chan_m70 = 4.74e-2
Vrest = -0.070

# Calculate secondary parameters
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
maxgl = Gin/sm_area

# Stimulus
preStimTime = 1
injectTime = 0.5
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 150e-12
Vleveli = -0.070
Vlevelf = -0.0

# Plotting stuff
lines = []
tplot = []
axes_list = []
sliders_list = []
fields = []


exec(open('Modelfunc.py').read())
# tinitial = time.time()
if os.path.exists(f'../../Output/{foldername}/Parametersdf.csv'):
    prevrun = 1
    Parameters = pd.read_csv(f'../../Output/{foldername}/Parametersdf.csv').tail(1)
else:
    prevrun = 0
    Parametersout = {}
    sm_area = sm_area
    Parameters['sm_diam'] = sm_diam
    Parameters['sm_len'] = sm_len
    Parameters['RM'] = 4
    Parameters['CM'] = CM
    Parameters['Em'] = Vrest
    Parameters['Ca_concCaBasal'] = 0.05e-3
    Parameters['Ca_conctau'] = 0.1
    Parameters['Ca_concB'] = 2291006575
    Parameters['Ca_L_changbar'] = 3
    Parameters['Ca_N_changbar'] = 3
    Parameters['Ca_T_changbar'] = 3
    Parameters['h_changbar'] = 0.25
    Parameters['K_A_changbar'] = 30
    Parameters['K_BK_changbar'] = 8
    Parameters['K_D_changbar'] = 0.05
    Parameters['K_DR_changbar'] = 3
    Parameters['K_M_changbar'] = 1.1
    Parameters['K_SK_changbar'] = 1
    Parameters['Na_changbar'] = 1000

[Parametersout, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=2, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'][0], Ca_concCaBasal=str(Parameters['Ca_concCaBasal'][0]), Ca_L_changbar=str(Parameters['Ca_L_changbar'][0]), Ca_N_changbar=str(Parameters['Ca_N_changbar'][0]), Ca_T_changbar=str(Parameters['Ca_T_changbar'][0]), h_changbar=str(Parameters['h_changbar'][0]), K_A_changbar=str(Parameters['K_A_changbar'][0]), K_BK_changbar=str(Parameters['K_BK_changbar'][0]), K_D_changbar=str(Parameters['K_D_changbar'][0]), K_DR_changbar=str(Parameters['K_DR_changbar'][0]), K_M_changbar=str(Parameters['K_M_changbar'][0]), K_SK_changbar=str(Parameters['K_SK_changbar'][0]), Na_changbar=str(Parameters['Na_changbar'][0]))

########### Defining some Matplotlib functions####################
def update(val):
    global Parameters
    global fig
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

    [Parametersout, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=2, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'][0], Ca_concCaBasal=str(Parameters['Ca_concCaBasal'][0]), Ca_L_changbar=str(Parameters['Ca_L_changbar'][0]), Ca_N_changbar=str(Parameters['Ca_N_changbar'][0]), Ca_T_changbar=str(Parameters['Ca_T_changbar'][0]), h_changbar=str(Parameters['h_changbar'][0]), K_A_changbar=str(Parameters['K_A_changbar'][0]), K_BK_changbar=str(Parameters['K_BK_changbar'][0]), K_D_changbar=str(Parameters['K_D_changbar'][0]), K_DR_changbar=str(Parameters['K_DR_changbar'][0]), K_M_changbar=str(Parameters['K_M_changbar'][0]), K_SK_changbar=str(Parameters['K_SK_changbar'][0]), Na_changbar=str(Parameters['Na_changbar'][0]))

    # sliders_list[0].set_val(Parametersout['Em'])
    # sliders_list[1].set_val(1/Parametersout['RM'])
    # sliders_list[8].set_val(Parametersout['h_changbar'])
    # sliders_list[8].set_val(1)
    # sliders_list[8].set_val(1)

    l.set_ydata(Vmvec)
    fig.canvas.draw_idle()

########### Initial setting up of axes############################
fig, ax = plt.subplots()
plt.subplots_adjust(top = 0.90, bottom=0.60)
l, = plt.plot(tvec, Vmvec, lw=2, color='red',label='in-silico')
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.title(f'Injected current = {Injectcurr}A')
plt.legend()
plt.axis([0, runtime, -0.090, 0.060])

axcolor = 'lightgoldenrodyellow'
for x in np.linspace(0.5,0.1,16):
    axes_list.append(plt.axes([0.1, x, 0.8, 0.02], facecolor=axcolor))

sliders_list.append(Slider(axes_list[0], 'Em', Parametersout['Em']-0.010, Parametersout['Em']+0.010, valinit=Parametersout['Em']))
sliders_list.append(Slider(axes_list[1], 'gl', mingl, Gin/sm_area, valinit=1/Parametersout['RM']))
sliders_list.append(Slider(axes_list[2], 'Ca_concCaBasal', Parametersout['Ca_concCaBasal']*0.5, Parametersout['Ca_concCaBasal']*2, valinit=Parametersout['Ca_concCaBasal']))
sliders_list.append(Slider(axes_list[3], 'Ca_conctau', Parametersout['Ca_conctau']*0.5, Parametersout['Ca_conctau']*2, valinit=Parametersout['Ca_conctau']))
sliders_list.append(Slider(axes_list[4], 'Ca_concB', Parametersout['Ca_concB']*0.5, Parametersout['Ca_concB']*2, valinit=Parametersout['Ca_concB']))
sliders_list.append(Slider(axes_list[5], 'Ca_L_changbar', Parametersout['Ca_L_changbar']*0.5, Parametersout['Ca_L_changbar']*2, valinit=Parametersout['Ca_L_changbar']))
sliders_list.append(Slider(axes_list[6], 'Ca_N_changbar', Parametersout['Ca_N_changbar']*0.5, Parametersout['Ca_N_changbar']*2, valinit=Parametersout['Ca_N_changbar']))
sliders_list.append(Slider(axes_list[7], 'Ca_T_changbar', Parametersout['Ca_T_changbar']*0.5, Parametersout['Ca_T_changbar']*2, valinit=Parametersout['Ca_T_changbar']))
sliders_list.append(Slider(axes_list[8], 'h_changbar', Parametersout['h_changbar']*0.5, Parametersout['h_changbar']*2, valinit=Parametersout['h_changbar']))
sliders_list.append(Slider(axes_list[9], 'K_A_changbar', Parametersout['K_A_changbar']*0.5, Parametersout['K_A_changbar']*2, valinit=Parametersout['K_A_changbar']))
sliders_list.append(Slider(axes_list[10], 'K_BK_changbar', Parametersout['K_BK_changbar']*0.5, Parametersout['K_BK_changbar']*2, valinit=Parametersout['K_BK_changbar']))
sliders_list.append(Slider(axes_list[11], 'K_D_changbar', Parametersout['K_D_changbar']*0.5, Parametersout['K_D_changbar']*2, valinit=Parametersout['K_D_changbar']))
sliders_list.append(Slider(axes_list[12], 'K_DR_changbar', Parametersout['K_DR_changbar']*0.5, Parametersout['K_DR_changbar']*2, valinit=Parametersout['K_DR_changbar']))
sliders_list.append(Slider(axes_list[13], 'K_M_changbar', Parametersout['K_M_changbar']*0.5, Parametersout['K_M_changbar']*2, valinit=Parametersout['K_M_changbar']))
sliders_list.append(Slider(axes_list[14], 'K_SK_changbar', Parametersout['K_SK_changbar']*0.5, Parametersout['K_SK_changbar']*2, valinit=Parametersout['K_SK_changbar']))
sliders_list.append(Slider(axes_list[15], 'Na_changbar', Parametersout['Na_changbar']*0.5, Parametersout['Na_changbar']*2, valinit=Parametersout['Na_changbar']))
for s in sliders_list:
    s.on_changed(update)

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
    Parametersdf = pd.DataFrame([Parametersout])
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
    Injectcurr = eval(text)
    [Parametersout, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=2, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'][0], Ca_concCaBasal=str(Parameters['Ca_concCaBasal'][0]), Ca_L_changbar=str(Parameters['Ca_L_changbar'][0]), Ca_N_changbar=str(Parameters['Ca_N_changbar'][0]), Ca_T_changbar=str(Parameters['Ca_T_changbar'][0]), h_changbar=str(Parameters['h_changbar'][0]), K_A_changbar=str(Parameters['K_A_changbar'][0]), K_BK_changbar=str(Parameters['K_BK_changbar'][0]), K_D_changbar=str(Parameters['K_D_changbar'][0]), K_DR_changbar=str(Parameters['K_DR_changbar'][0]), K_M_changbar=str(Parameters['K_M_changbar'][0]), K_SK_changbar=str(Parameters['K_SK_changbar'][0]), Na_changbar=str(Parameters['Na_changbar'][0]))
    print(Injectcurr)
currbut.on_submit(curr)

plt.show(block=False)
