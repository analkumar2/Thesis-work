## exec(open('CVvsMean.py').read())

# import moose
# import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
# import xmltodict
import sys
import os
import io
# import importlib
# # import MOOSEModel_17_somamulti as mm
# import pickle
# # import featuresv26_nonallen_uniform as fts
# import moose
# # import plotexpv2 as pex
# from copy import deepcopy
import numpy.random as nr
from multiprocessing import Pool
import time
# import argparse
from pprint import pprint
import pickle
from copy import deepcopy
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
# import warnings
# warnings.filterwarnings("ignore", category=RuntimeWarning)
# from sklearn.linear_model import LinearRegression
# import gc

# from sklearn.linear_model import LinearRegression

from neuron import h,gui

elecid_ori = None
elecPlotDt = 0.0001
elecDt = 0.0001

stim_start = 1
stim_end = 1.9
totalsec = 2.5
samprate = 20000

nan = 10000

def get_Vmvec(tI_II):
    tI = tI_II[0]
    II = tI_II[1]
    v_vec = h.Vector()             # Membrane potential vector
    t_vec = h.Vector()             # Time stamp vector
    v_vec.record(h.soma[0](0.5)._ref_v)
    t_vec.record(h._ref_t)

    stim = h.IClamp(h.soma[0](0.5))

    currh = h.Vector(tI*1e9)
    th = h.Vector(tI * 1e3)

    pprint(np.array(th))
    pprint(np.array(currh))

    stim.delay = 0
    stim.dur = 1e9
    currh.play(stim._ref_amp, th, 1)

    h.finitialize()
    h.tstop = 1000
    h.run()

    pprint(np.array(t_vec))
    pprint(np.array(v_vec))

    plt.plot(np.array(t_vec)*1e-3,np.array(v_vec)*1e-3)
    plt.show()

    # spiket = processVmvec(np.array(t_vec)*1e-3,np.array(v_vec)*1e-3)

    # return spiket
    return 1

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

def calcCVjit():
    tI_list = []
    II_list = []
    for i in range(3):
        curr = np.zeros(int(samprate*totalsec))
        curr[int(1*samprate):int(1.9*samprate)] = 150e-12
        noise = np.random.normal(0,20e-12,int(samprate*totalsec))
        curr = curr + noise
        t = np.linspace(0,totalsec, int(totalsec*samprate))
        tI_list.append(t)
        II_list.append(curr)

    # tempspiket = get_Vmvec([tI_list[0],II_list[0]])
    # if len(tempspiket)<2:
    #     print('Too few spikes 1')
    #     return [np.nan,np.nan]

    spiket_list = []
    ISI500list = []
    pool = Pool(processes=os.cpu_count()-4) #opening processes
    A = pool.map(get_Vmvec, zip(tI_list, II_list))
    # for a in A:
    #     spiket = a
    #     if len(spiket)<2:
    #         print('Too few spikes 2')
    #         continue
    #         # return [np.nan,np.nan]

    #     if len([i for i in spiket if i>(stim_start+0.5)])<1 or len([i for i in spiket if i<(stim_start+0.5)])<1:
    #         print('No spikes around 0.5s mark')
    #         continue
    #         # return [np.nan,np.nan]

    #     ISI500list.append(min([i for i in spiket if i>(stim_start+0.5)]) - max([i for i in spiket if i<(stim_start+0.5)]))
    #     spiket_list.append(spiket)
    # # plt.show()
    # pool.terminate()

    # CV500 = calcCV(ISI500list)
    # jitter = calcjit(spiket_list)

    # f.write(f'spiketlist = {spiket_list} \n \n')
    # f.close()

    # if len(ISI500list)<500 or len(spiket_list)<500:
    #     return [np.nan,np.nan]
    # else:
    #     return [CV500, jitter]
    return [1,2]


def main():
    h('load_file("CombeEtAl2018/simplestim.hoc")')
    CV500, jit = calcCVjit()
    print(CV500, jit)


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('Chan', type=str)
    # args = parser.parse_args()
    # main(args.Chan)
    #fig = pickle.load(open('CVvsMean.pkl', 'rb'))
    #plt.show()
    main()
