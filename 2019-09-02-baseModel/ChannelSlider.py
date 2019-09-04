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

# Define constants not to be changed
F = 96485.3329
elecPlotDt = 0.00005

# Define model parameters
sm_diam=60e-6; sm_len=60e-6
CM = 0.01
Em = -0.070
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
Vleveli = -0.070
Vlevelf = -0.0


exec(open('Modelfunc.py').read())
# tinitial = time.time()
[Parameters, characteristics, Vmvec, Ivec, Cavec, tvec] = Modelfunc(runfor=2, stimul='Vclamp', Injectcurr=150e-12, gl = 0.25, Ca_concCaBasal='0.05e-3', Ca_L_changbar='3', Ca_N_changbar='3', Ca_T_changbar='3', h_changbar='0.25', K_A_changbar='30', K_BK_changbar='8', K_D_changbar='0.05', K_DR_changbar='3', K_M_changbar='1.1', K_SK_changbar='1', Na_changbar='1000')
# print(time.time()-tinitial)

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
