## exec(open('CVvsMean.py').read())

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import xmltodict
import sys
import os
import io
import importlib
import MOOSEModel_17_somamulti as mm
import pickle
# import featuresv26_nonallen_uniform as fts
import moose
# import plotexpv2 as pex
from copy import deepcopy
import numpy.random as nr
from multiprocessing import Pool
import time
import argparse
from pprint import pprint
import pickle
from copy import deepcopy
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.linear_model import LinearRegression
import gc

# from sklearn.linear_model import LinearRegression

elecid_ori = None
elecPlotDt = 0.0001
elecDt = 0.0001

stim_start = 1
stim_end = 1.9
totalsec = 2.5
samprate = 20000

nan = 10000
# model = {'Error': 470.08426631177923, 'Parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.19e-10, 'Rm': 198905511.35013187, 'Em': -0.05608805404962253}, 'Channels': {'Na_Chan': {'Gbar': 4.932428436178534e-06, 'Erev': 0.06, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_Custom4', 'KineticVars': {'m_vhalf_inf': -0.0359019320487981, 'm_slope_inf': 0.00788319407034278, 'm_A': -0.044359772011543515, 'm_B': 0.02, 'm_C': 0.0161, 'm_D': 0.0547, 'm_E': 0.0311, 'm_F': 0.00064, 'h_vhalf_inf': -0.050316776512570736, 'h_slope_inf': -0.005219098993311621, 'h_A': -0.04763346729779447, 'h_B': 0.003464, 'h_C': 0.0, 'h_D': 0.0262, 'h_E': 0.00854, 'h_F': 0.3069864140957819, 's_vhalf_inf': -0.04492406530182888, 's_slope_inf': -0.010911310412463028, 's_A': 1, 's_B': 0.001, 's_C': 0, 's_D': 0.6152461928009552, 's_E': 0.001, 's_F': 1}}, 'K_DR_Chan': {'Gbar': 7.658644948611434e-07, 'Erev': -0.09, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_Custom3', 'KineticVars': {'n_vhalf_inf': 0.01731374470366132, 'n_slope_inf': 0.01471839705516013, 'n_A': 0.014985238813908366, 'n_E': 0.020993521832043036, 'n_F': 0.01380360065197889}}, 'K_A_Chan': {'Gbar': 1.8215306367504489e-07, 'Erev': -0.09, 'Kinetics': '../../Compilations/Kinetics/K_A_Chan_Custom3', 'KineticVars': {'n_vhalf_inf': 0.02149776726410334, 'n_slope_inf': 0.012912009256498078, 'n_A': -0.024208971400308105, 'n_B': 0.08860901204015699, 'n_C': 0, 'n_D': 0, 'n_E': 0.0060496419089759845, 'n_F': 0.004765782060316019, 'l_vhalf_inf': -0.04389847493434692, 'l_slope_inf': -0.025778506595107165, 'l_min': 0.002, 'l_m': 0.3418388911022859, 'l_cm': 0.05}}, 'K_M_Chan': {'Gbar': 1.6808195467313064e-10, 'Erev': -0.09, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_Custom1', 'KineticVars': {'factor': 3.3e-05}}, 'K_SK_Chan': {'Gbar': 2.056391414204614e-15, 'Erev': -0.09, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_Custom3'}, 'Ca_T_Chan': {'Gbar': 1.6721696312505114e-15, 'Erev': 0.12, 'Kinetics': '../../Compilations/Kinetics/Ca_T_Chan_Custom1'}}, 'Ca_Conc': {'Ca_B': 1800000000.0, 'Ca_tau': 0.15, 'Ca_base': 5e-05, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Common)'}}, 'Scores': {}, 'Score': {'E_rest_0': 0, 'Input resistance': -1.9528889095844233, 'Cell capacitance': 0, 'AP1_amp_1.5e-10': 0, 'APp_amp_1.5e-10': 0, 'APavgpratio_amp_1.5e-10': -3.0984125901593194, 'AP1_width_1.5e-10': -3.6253310169109465, 'APp_width_1.5e-10': -2.4509513665268226, 'AP1_thresh_1.5e-10': 0, 'APp_thresh_1.5e-10': 0, 'ISI1_1.5e-10': 5.549977871150843, 'ISIl_1.5e-10': 0, 'ISIavg_1.5e-10': 0, 'freq_1.5e-10': 0, 'Adptn_id_1.5e-10': -6.337688639200551, 'fAHP_AP1_amp_1.5e-10': 0, 'fAHP_APp_amp_1.5e-10': -1.601930907299259, 'mAHP_APp_amp_1.5e-10': -2.7346819933400606, 'mAHP_APp_time_1.5e-10': 0, 'AHP_AP1_amp_1.5e-10': 0, 'AHP_APp_amp_1.5e-10': -2.57430297658575, 'AHP_APp_time_1.5e-10': 0, 'Upstroke_AP1_time_1.5e-10': 1.7639475967188118, 'Upstroke_APp_time_1.5e-10': 4.166816666667142, 'Upstroke_AP1_amp_1.5e-10': -1.103959350109573, 'Upstroke_APp_amp_1.5e-10': -2.960660499806755, 'Upstroke_AP1_value_1.5e-10': 1.4455993645508278, 'Upstroke_APp_value_1.5e-10': 2.9806177283258517, 'Downstroke_AP1_time_1.5e-10': -3.901371573204427, 'Downstroke_APp_time_1.5e-10': -1.166711666665328, 'Downstroke_AP1_amp_1.5e-10': 4.945860117472832, 'Downstroke_APp_amp_1.5e-10': 0, 'Downstroke_AP1_value_1.5e-10': -3.376310966584173, 'Downstroke_APp_value_1.5e-10': -2.2735722411322983, 'UpDn_AP1_ratio_1.5e-10': 0, 'UpThr_AP1_diff_1.5e-10': -2.028035043099396, 'UpThr_APp_diff_1.5e-10': -3.054102872346163, 'offset_1.5e-10': -2.0129256585474393, 'AP1_amp_3e-10': 0, 'APp_amp_3e-10': -3.2455115126952085, 'APavgpratio_amp_3e-10': 0, 'AP1_width_3e-10': -3.742272086154128, 'APp_width_3e-10': -1.368315651114673, 'AP1_thresh_3e-10': 0, 'APp_thresh_3e-10': 0, 'ISI1_3e-10': 6.997038752423039, 'ISIl_3e-10': 0, 'ISIavg_3e-10': 0, 'freq_3e-10': 0, 'Adptn_id_3e-10': -5.2066952295950335, 'fAHP_AP1_amp_3e-10': 0, 'fAHP_APp_amp_3e-10': -2.0616055775966013, 'mAHP_APp_amp_3e-10': -1.9164540887580153, 'mAHP_APp_time_3e-10': 0, 'AHP_AP1_amp_3e-10': 0, 'AHP_APp_amp_3e-10': -1.9423259205205055, 'AHP_APp_time_3e-10': 0, 'Upstroke_AP1_time_3e-10': 3.1454726305580336, 'Upstroke_APp_time_3e-10': 0, 'Upstroke_AP1_amp_3e-10': 0, 'Upstroke_APp_amp_3e-10': -6.279148420998243, 'Upstroke_AP1_value_3e-10': 1.2877783413072783, 'Upstroke_APp_value_3e-10': 0, 'Downstroke_AP1_time_3e-10': -4.151434269793974, 'Downstroke_APp_time_3e-10': 0, 'Downstroke_AP1_amp_3e-10': 5.5410446811235525, 'Downstroke_APp_amp_3e-10': 0, 'Downstroke_AP1_value_3e-10': -4.102321125163147, 'Downstroke_APp_value_3e-10': 0, 'UpDn_AP1_ratio_3e-10': 0, 'UpThr_AP1_diff_3e-10': 0, 'UpThr_APp_diff_3e-10': 0, 'offset_3e-10': -1.1159060878515004, 'freq300to150ratio': 1.4756342221364624}, 'Features': {'Sampling rate': 10000.4, 'stim_start': 1, 'stim_end': 1.5, 'E_rest_0': -0.06716766200787411, 'E_rest_m25': -0.06716847233440179, 'E_rest_150': -0.06716847233440967, 'E_rest_300': -0.0671684723344067, 'Input resistance': 99007520.80599822, 'Cell capacitance': 1.189516331831558e-10, 'AP1_amp_1.5e-10': 0.11489228726817852, 'APp_amp_1.5e-10': 0.10309974181042947, 'APavgpratio_amp_1.5e-10': 1.011636128160689, 'AP1_width_1.5e-10': 0.000600000000000156, 'APp_width_1.5e-10': 0.0008000000000001339, 'AP1_thresh_1.5e-10': -0.04792232083616102, 'APp_thresh_1.5e-10': -0.04408274251225343, 'ISI1_1.5e-10': 0.04059999999999997, 'ISIl_1.5e-10': 0.04730000000000012, 'ISIavg_1.5e-10': 0.04504000000000001, 'freq_1.5e-10': 22.0, 'Adptn_id_1.5e-10': 0.14164904862579564, 'fAHP_AP1_amp_1.5e-10': 0.013265121619563294, 'fAHP_APp_amp_1.5e-10': 0.01462801477249473, 'mAHP_APp_amp_1.5e-10': 0.010687823364556537, 'mAHP_APp_time_1.5e-10': 0.01760000000000006, 'AHP_AP1_amp_1.5e-10': 0.010698216494546685, 'AHP_APp_amp_1.5e-10': 0.010687823364556537, 'AHP_APp_time_1.5e-10': 0.017599296028158875, 'Upstroke_AP1_time_1.5e-10': -0.00019999999999997797, 'Upstroke_APp_time_1.5e-10': -0.00019999999999997797, 'Upstroke_AP1_amp_1.5e-10': 0.052720853802266814, 'Upstroke_APp_amp_1.5e-10': 0.053067399588552736, 'Upstroke_AP1_value_1.5e-10': 458.853813569696, 'Upstroke_APp_value_1.5e-10': 310.57488393468816, 'Downstroke_AP1_time_1.5e-10': 0.0, 'Downstroke_APp_time_1.5e-10': 0.000300000000000189, 'Downstroke_AP1_amp_1.5e-10': 0.11489228726817852, 'Downstroke_APp_amp_1.5e-10': 0.08283788538314664, 'Downstroke_AP1_value_1.5e-10': -150.80847062436462, 'Downstroke_APp_value_1.5e-10': -103.94659286355234, 'UpDn_AP1_ratio_1.5e-10': 2.9878313023917076, 'UpThr_AP1_diff_1.5e-10': 0.03347470230401816, 'UpThr_APp_diff_1.5e-10': 0.029981669766396494, 'offset_1.5e-10': 0.010699343217505, 'AP1_amp_3e-10': 0.11501530256693318, 'APp_amp_3e-10': 0.0771327904013039, 'APavgpratio_amp_3e-10': 1.0320725001695672, 'AP1_width_3e-10': 0.0005999999999999339, 'APp_width_3e-10': 0.0013000000000000789, 'AP1_thresh_3e-10': -0.05174071145709095, 'APp_thresh_3e-10': -0.0420327284835539, 'ISI1_3e-10': 0.01859999999999995, 'ISIl_3e-10': 0.026799999999999935, 'ISIavg_3e-10': 0.024055000000000003, 'freq_3e-10': 42.0, 'Adptn_id_3e-10': 0.30597014925373156, 'fAHP_AP1_amp_3e-10': 0.014689867767231538, 'fAHP_APp_amp_3e-10': 0.018489331086231908, 'mAHP_APp_amp_3e-10': 0.017472600616231806, 'mAHP_APp_time_3e-10': 0.011900000000000022, 'AHP_AP1_amp_3e-10': 0.014395128935768262, 'AHP_APp_amp_3e-10': 0.017472600616231806, 'AHP_APp_time_3e-10': 0.011899524019039239, 'Upstroke_AP1_time_3e-10': -9.999999999998899e-05, 'Upstroke_APp_time_3e-10': -0.000400000000000178, 'Upstroke_AP1_amp_3e-10': 0.05988333479158259, 'Upstroke_APp_amp_3e-10': 0.047241378190381905, 'Upstroke_AP1_value_3e-10': 434.8558465335387, 'Upstroke_APp_value_3e-10': 103.05596623503826, 'Downstroke_AP1_time_3e-10': 9.999999999998899e-05, 'Downstroke_APp_time_3e-10': 0.0004999999999999449, 'Downstroke_AP1_amp_3e-10': 0.11396195824118535, 'Downstroke_APp_amp_3e-10': 0.06324332216898552, 'Downstroke_AP1_value_3e-10': -159.11017309229854, 'Downstroke_APp_value_3e-10': -45.730396281957546, 'UpDn_AP1_ratio_3e-10': 2.253555066517058, 'UpThr_AP1_diff_3e-10': 0.044455573914266834, 'UpThr_APp_diff_3e-10': 0.0221056343395291, 'offset_3e-10': 0.01756057160002228, 'freq300to150ratio': 1.9090909090909092}}

# mm.plotModel(model)

# exec(open("Combined100models.py").read())
from Combined100models import Models
modelname = 'Model4'
fullModel = deepcopy(Models[modelname])
mm.plotModel(
            fullModel,
            CurrInjection=150e-12,
            vClamp=None,
            refreshKin=True,
            Truntime=0.01,
            syn=True, synwg=0.01, synfq=0.5
        )
plt.close('all')


def get_Vmvec(fullModel_tI_II):
    fullModel = deepcopy(fullModel_tI_II[0])
    tI = fullModel_tI_II[1]
    II = fullModel_tI_II[2]
    tempt, tempv, Ca = mm.runModel(
            fullModel,
            CurrInjection=150e-12,
            vClamp=None,
            refreshKin=False,
            Truntime=0.01,
            syn=True, synwg=0.01, synfq=0.43
        )
    moose.delete("/model/stims/stim0")
    stimtable = moose.StimulusTable("/model/stims/stim2")
    soma = moose.element("/model/elec/soma")
    moose.connect(stimtable, "output", soma, "setInject")

    stimtable.vector = II
    stimtable.stepSize = (
        0  # This forces use of current time as x value for interpolation
    )
    stimtable.stopTime = tI[-1]

    Tdur = tI[-1]
    moose.reinit()
    moose.start(tI[-1])
    Vmvec = moose.element("/model/graphs/plot0").vector
    tvec = moose.element("/Graphs/plott").vector

    spiket = processVmvec(tvec,Vmvec)

    return spiket

def processVmvec(tvec,Vmvec):
    tt = tvec
    vv = Vmvec
    I = 150e-12
    ii = np.zeros(len(tt))
    ii[(tt >= stim_start) & (tt <= stim_end)] = I
    sweep_ext = EphysSweepFeatureExtractor(
        t=tt,
        v=vv * 1e3,
        i=ii * 1e12,
        filter=len(tt) / tt[-1] / 2500,
        start=stim_start,
        end=stim_end,
    )
    try:
        sweep_ext.process_spikes()
    except ValueError:
        return []
    spiket = sweep_ext.spike_feature("peak_t")
    return spiket

def calcCV(ISIlist):
    return np.std(ISIlist)/np.mean(ISIlist)

def calcjit(spiket_list):
    minspikes=100
    for spiket in spiket_list:
        if len(spiket) < minspikes:
            minspikes = len(spiket)

    for i in range(len(spiket_list)):
        spiket_list[i] = spiket_list[i][:minspikes]

    jitter = np.nanstd(spiket_list, 0)
    spikemean = np.nanmean(spiket_list, 0)
    # print(jitter, spikemean)

    x = spikemean.reshape((-1, 1))
    y = jitter
    model = LinearRegression().fit(x, y)
    # print(model.score(x, y), model.intercept_, model.coef_)

    return model.coef_[0]

def calcCVjit(fullModel):
    f = open('CVvsMean_True_channels4_spiket.py', 'a+')
    f.write(f'fullModel = {fullModel} \n \n')
    tI_list = []
    II_list = []
    for i in range(1000):
        curr = np.zeros(int(samprate*totalsec))
        curr[int(1*samprate):int(1.9*samprate)] = 150e-12
        noise = np.random.normal(0,20e-12,int(samprate*totalsec))
        curr = curr + noise
        t = np.linspace(0,totalsec, int(totalsec*samprate))
        tI_list.append(t)
        II_list.append(curr)

    tempspiket = get_Vmvec([fullModel,tI_list[0],II_list[0]])
    if len(tempspiket)<2:
        print('Too few spikes 1')
        return [np.nan,np.nan]

    spiket_list = []
    ISI500list = []
    pool = Pool(processes=os.cpu_count()-10) #opening processes
    A = pool.map(get_Vmvec, zip(np.repeat(fullModel, 1000), tI_list, II_list))
    for a in A:
        # tvec,Vmvec = a
        # plt.plot(tvec, Vmvec)
        # plt.show()
        # spiket = processVmvec(tvec,Vmvec)
        spiket = a
        if len(spiket)<2:
            print('Too few spikes 2')
            continue
            # return [np.nan,np.nan]

        if len([i for i in spiket if i>(stim_start+0.5)])<1 or len([i for i in spiket if i<(stim_start+0.5)])<1:
            print('No spikes around 0.5s mark')
            continue
            # return [np.nan,np.nan]

        ISI500list.append(min([i for i in spiket if i>(stim_start+0.5)]) - max([i for i in spiket if i<(stim_start+0.5)]))
        spiket_list.append(spiket)
    # plt.show()
    pool.terminate()

    CV500 = calcCV(ISI500list)
    jitter = calcjit(spiket_list)

    f.write(f'spiketlist = {spiket_list} \n \n')
    f.close()

    if len(ISI500list)<500 or len(spiket_list)<500:
        return [np.nan,np.nan]
    else:
        return [CV500, jitter]


def main(Channame):
    newModel = deepcopy(fullModel)
    CV500, jit = calcCVjit(newModel)
    print(CV500, jit)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Chan', type=str)
    args = parser.parse_args()
    main(args.Chan)
    #fig = pickle.load(open('CVvsMean.pkl', 'rb'))
    #plt.show()
