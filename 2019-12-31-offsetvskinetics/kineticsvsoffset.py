## exec(open('kineticsvsoffset.py').read())

import moose
import MOOSEModel_4 as mm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('poster')


def somefeatures(currAmp,Vtrace,stim_start,stim_end,sampling_rate):
    features = {}
    t = np.linspace(0,len(Vtrace)/sampling_rate, len(Vtrace))
    t = np.array(t) # in s

    # features[f'E_rest'] = np.median(Vtrace[t<=stim_start])
    # features[f'freq'] = len(scs.find_peaks(Vtrace[(t>=stim_start)&(t<=stim_end)], height=-0.03, prominence=0.01)[0])/(stim_end - stim_start)
    features[f'offset'] = np.nanmin(Vtrace[(t<=stim_end-0.001) & (t>=stim_start+0.400)])
    return features

Na_Chan_X_vhalf, K_DR_Chan_X_slope, K_DR_Chan_X_F = -0.03626557246957454, 0.004293268274831886, 0.016008788932111953
initModel = {'Error': 5.537655149381235, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.8904569830116381e-10, 'Rm': 787836790.3806579, 'Em': -0.07084106424846002}, 'Channels': {'Na_Chan': {'Gbar': 3.771544601682249e-07, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.08021304143433096, 'gateX': [Na_Chan_X_vhalf, 0.004367504181119175, -0.03137093478722176, 0.01953703125568758, 0.06340763317450275, 0.06855864702943262, 0.018880962188814493, 0.0009427141127737687], 'gateY': [-0.032027246727940814, -0.0024429964267569876, -0.045644301575264404, 0.006059039958509024, 0.09933568933519707, 0.06112850733311163, 0.009885697510926968, 0.02501976667692509]}, 'K_DR_Chan': {'Gbar': 1.5159047773904547e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.08030796790883893, 'gateX': [0.028026663300130376, K_DR_Chan_X_slope, -0.015514359010036357, 0.02756184788269313, 0.09468903432281536, 0.028297116277257783, 0.03852341077851185, K_DR_Chan_X_F]}}}}

# mm.plotModel(initModel, 300e-12)

offset_l = []
Na_Chan_X_vhalf_l = np.arange(-0.036,-0.023, 0.001)
for Na_Chan_X_vhalf in Na_Chan_X_vhalf_l:
    K_DR_Chan_X_slope, K_DR_Chan_X_F = 0.004293268274831886, 0.016008788932111953
    initModel = {'Error': 5.537655149381235, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.8904569830116381e-10, 'Rm': 787836790.3806579, 'Em': -0.07084106424846002}, 'Channels': {'Na_Chan': {'Gbar': 3.771544601682249e-07, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.08021304143433096, 'gateX': [Na_Chan_X_vhalf, 0.004367504181119175, -0.03137093478722176, 0.01953703125568758, 0.06340763317450275, 0.06855864702943262, 0.018880962188814493, 0.0009427141127737687], 'gateY': [-0.032027246727940814, -0.0024429964267569876, -0.045644301575264404, 0.006059039958509024, 0.09933568933519707, 0.06112850733311163, 0.009885697510926968, 0.02501976667692509]}, 'K_DR_Chan': {'Gbar': 1.5159047773904547e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.08030796790883893, 'gateX': [0.028026663300130376, K_DR_Chan_X_slope, -0.015514359010036357, 0.02756184788269313, 0.09468903432281536, 0.028297116277257783, 0.03852341077851185, K_DR_Chan_X_F]}}}}
    ttrace,Vtrace = mm.runModel(initModel, 150e-12)
    features = somefeatures(150e-12,Vtrace,stim_start=1,stim_end=1.5,sampling_rate=len(ttrace)/ttrace[-1])
    offset_l.append(features[f'offset'])
    # plt.plot(ttrace,Vtrace, label=Na_Chan_X_vhalf)

plt.plot(Na_Chan_X_vhalf_l, offset_l, label='Sodium channel activation gate half activation potential')

offset_l = []
Na_Chan_Y_vhalf_l = np.arange(-0.036,-0.023, 0.001)
for Na_Chan_Y_vhalf in Na_Chan_Y_vhalf_l:
    K_DR_Chan_X_slope, K_DR_Chan_X_F = 0.004293268274831886, 0.016008788932111953
    initModel = {'Error': 5.537655149381235, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.8904569830116381e-10, 'Rm': 787836790.3806579, 'Em': -0.07084106424846002}, 'Channels': {'Na_Chan': {'Gbar': 3.771544601682249e-07, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.08021304143433096, 'gateX': [-0.03626557246957454, 0.004367504181119175, -0.03137093478722176, 0.01953703125568758, 0.06340763317450275, 0.06855864702943262, 0.018880962188814493, 0.0009427141127737687], 'gateY': [Na_Chan_Y_vhalf, -0.0024429964267569876, -0.045644301575264404, 0.006059039958509024, 0.09933568933519707, 0.06112850733311163, 0.009885697510926968, 0.02501976667692509]}, 'K_DR_Chan': {'Gbar': 1.5159047773904547e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.08030796790883893, 'gateX': [0.028026663300130376, K_DR_Chan_X_slope, -0.015514359010036357, 0.02756184788269313, 0.09468903432281536, 0.028297116277257783, 0.03852341077851185, K_DR_Chan_X_F]}}}}
    ttrace,Vtrace = mm.runModel(initModel, 150e-12)
    features = somefeatures(150e-12,Vtrace,stim_start=1,stim_end=1.5,sampling_rate=len(ttrace)/ttrace[-1])
    offset_l.append(features[f'offset'])
    # plt.plot(ttrace,Vtrace, label=Na_Chan_Y_vhalf)

plt.plot(Na_Chan_Y_vhalf_l, offset_l, label='Sodium channel inactivation gate half activation potential')

offset_l = []
K_DR_Chan_X_vhalf_l = np.arange(0.020,0.040, 0.001)
for K_DR_Chan_X_vhalf in K_DR_Chan_X_vhalf_l:
    K_DR_Chan_X_slope, K_DR_Chan_X_F = 0.004293268274831886, 0.016008788932111953
    initModel = {'Error': 5.537655149381235, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.8904569830116381e-10, 'Rm': 787836790.3806579, 'Em': -0.07084106424846002}, 'Channels': {'Na_Chan': {'Gbar': 3.771544601682249e-07, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.08021304143433096, 'gateX': [-0.03626557246957454, 0.004367504181119175, -0.03137093478722176, 0.01953703125568758, 0.06340763317450275, 0.06855864702943262, 0.018880962188814493, 0.0009427141127737687], 'gateY': [-0.032027246727940814, -0.0024429964267569876, -0.045644301575264404, 0.006059039958509024, 0.09933568933519707, 0.06112850733311163, 0.009885697510926968, 0.02501976667692509]}, 'K_DR_Chan': {'Gbar': 1.5159047773904547e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.08030796790883893, 'gateX': [K_DR_Chan_X_vhalf, K_DR_Chan_X_slope, -0.015514359010036357, 0.02756184788269313, 0.09468903432281536, 0.028297116277257783, 0.03852341077851185, K_DR_Chan_X_F]}}}}
    ttrace,Vtrace = mm.runModel(initModel, 150e-12)
    features = somefeatures(150e-12,Vtrace,stim_start=1,stim_end=1.5,sampling_rate=len(ttrace)/ttrace[-1])
    offset_l.append(features[f'offset'])
    # plt.plot(ttrace,Vtrace, label=K_DR_Chan_X_vhalf)

plt.plot(K_DR_Chan_X_vhalf_l, offset_l, label='Delayed rectifier potassium channel activation gate half activation potential')


plt.title('The after-hyperpolarization peak potential')
plt.xlabel('half activation Voltage (V)')
plt.ylabel('after-hyperpolarization peak potential (V)')
plt.legend()
plt.show()
