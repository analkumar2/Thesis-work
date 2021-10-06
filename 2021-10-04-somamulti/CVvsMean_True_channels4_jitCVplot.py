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
# import featuresv26_nonallen as fts
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
nan = np.nan
exec(open('CVvsMean_True_channels4_jitCV.py').read())

gbarratio = [0,0.1,0.2,0.5,0.75,0.9,1,1.1,1.5,2,3,5,10]

fig, axs = plt.subplots(1, 2)
axs[0].plot(gbarratio[:len(CV500_Na)], CV500_Na, label='Na', c='black')
# axs[0].plot(gbarratio[:len(CV500_K_DR)], CV500_K_DR, label='K_DR', c='yellow')
axs[0].plot(gbarratio[:len(CV500_K_A)], CV500_K_A, label='K_A', c='blue')
axs[0].plot(gbarratio[:len(CV500_K_M)], CV500_K_M, label='K_M', c='green')
# # axs[0].plot(gbarratio[:len(CV500_Na)], CV500_h, label='h', c='grey')
axs[0].plot(gbarratio[:len(CV500_K_SK)], CV500_K_SK, label='K_SK', c='red')
axs[0].plot(gbarratio[:len(CV500_Ca_L)], CV500_Ca_L, label='Ca_L', c='purple')

axs[1].plot(gbarratio[:len(jit_Na)], jit_Na, label='Na', c='black')
# axs[1].plot(gbarratio[:len(jit_K_DR)], jit_K_DR, label='K_DR', c='yellow')
axs[1].plot(gbarratio[:len(jit_K_A)], jit_K_A, label='K_A', c='blue')
axs[1].plot(gbarratio[:len(jit_K_M)], jit_K_M, label='K_M', c='green')
# # axs[1].plot(gbarratio[:len(jit_Na)], jit_h, label='h', c='grey')
axs[1].plot(gbarratio[:len(jit_K_SK)], jit_K_SK, label='K_SK', c='red')
axs[1].plot(gbarratio[:len(jit_Ca_L)], jit_Ca_L, label='Ca_L', c='purple')

axs[0].legend()
axs[1].legend()
axs[0].set_xlabel('gbar ratio')
axs[0].set_ylabel('CV500')
axs[1].set_xlabel('gbar ratio')
axs[1].set_ylabel('jitter slope')

fig.tight_layout()

plt.show()