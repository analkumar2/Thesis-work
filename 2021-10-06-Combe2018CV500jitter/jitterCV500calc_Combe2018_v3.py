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
import numpy.random as nr
from multiprocessing import Pool
import time
# import argparse
from pprint import pprint
import pickle
from copy import deepcopy
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
from sklearn.linear_model import LinearRegression
import gc

from neuron import h,gui

# os.system('mpiexec -n 4 python3 spiket_Combe2018.py')

elecid_ori = None
elecPlotDt = 0.0001
elecDt = 0.0001

stim_start = 1
stim_end = 1.9
totalsec = 2.5
samprate = 20000

nan = 10000

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
    A = np.load('spiket.npy', allow_pickle=True)

    spiket_list = []
    ISI500list = []
    for a in A:
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

    CV500 = calcCV(ISI500list)
    jitter = calcjit(spiket_list)

    if len(ISI500list)<500 or len(spiket_list)<500:
        return [np.nan,np.nan]
    else:
        return [CV500, jitter]


def main():
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
