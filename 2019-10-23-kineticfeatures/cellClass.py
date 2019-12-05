#exec(open('cellClass.py').read())

import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import signal

reader = h5py.File('../../Raw_data/Channelpedia/Kv1.1/DataKv1.1RatCHO/rCell2241.nwb')

class KineticsCell:
    '''
    Input a proper h5pyreader object downloaded from channelpedia. And the repetion number.
    Please confirm that the following protocols were followed:
    (https://channelpedia.epfl.ch/expdata/stimulus)
    Activation (#3)
    Deactivation (#8)
    Inactivation (#10)
    Ramp (#12)
    AP (#31)
    Recovery (#36)

    '''
    def __init__(self, h5pyreader, repetition='repetition1'):
        self.ID = h5pyreader['general']['cell_id'][:][0][0]
        self.species = reader['general']['channel_info']['species'][:][0]
        self.host_cell = reader['general']['channel_info']['species'][:][0]
        self.ion_channel = reader['general']['channel_info']['ion_channel'][:][0]
        self.repetition = repetition

    def act_params(self):
        samprate = 10000
        data = reader['acquisition']['timeseries']['Activation']['repetitions'][repetition]['data']
        baseline = np.mean(data[:,1][1000:5990])
        act_params_dict = {}
        levels = np.arange(-0.090,0.090,0.010)
        levelnum = 0
        for level in levels:
            max = np.max(data[:,levelnum][995:1990] - baseline) if level>=0 else -np.max(baseline - data[:,levelnum][995:1990]) #995 because of capacitive artifacts till 994
            max0 = max0+max
            steady50msavg = np.mean(data[:,levelnum][5490:5990] - baseline) if level>=0 else -np.mean(baseline - data[:,levelnum][5490:5990])
            maxtimeidx = data[:,levelnum][990:1990].argmax() if level>=0 else data[:,levelnum][990:1990].argmin()
            if maxtimeidx<=5:
                risetau = np.nan
                decaytau = np.nan
            else:
                risetau = np.abs((max)*0.865+baseline-data[:,levelnum][990:maxtimeidx+990]).argmin()/samprate/2    #2tau is used because near the beginning, there is capacitance artifact
                decaytau = np.abs((max-steady50msavg)*0.368+steady50msavg-data[:,levelnum][maxtimeidx+990:5990]).argmin()/samprate
            levelnum +=1
            act_params_dict[f'{level:.2f}'+'_max']= max
            act_params_dict[f'{level:.2f}'+'_risetau']= risetau
            act_params_dict[f'{level:.2f}'+'_decaytau']= decaytau
        return act_params_dict

    def deact_params(self):
        samprate = 10000
        data = reader['acquisition']['timeseries']['Deactivation']['repetitions'][repetition]['data']
        baseline = np.mean(data[:,1][3990:5990])
        act_params_dict = {}
        levels = np.arange(-0.080,0.040,0.010)
        levelnum = 0
        for level in levels:
            max = np.max(data[:,levelnum][995:1990] - baseline) if level>=0 else -np.max(baseline - data[:,levelnum][995:1990]) #995 because of capacitive artifacts till 994
            max0 = max0+max
            steady50msavg = np.mean(data[:,levelnum][5490:5990] - baseline) if level>=0 else -np.mean(baseline - data[:,levelnum][5490:5990])
            maxtimeidx = data[:,levelnum][990:1990].argmax() if level>=0 else data[:,levelnum][990:1990].argmin()
            if maxtimeidx<=5:
                risetau = np.nan
                decaytau = np.nan
            else:
                risetau = np.abs((max)*0.865+baseline-data[:,levelnum][990:maxtimeidx+990]).argmin()/samprate/2    #2tau is used because near the beginning, there is capacitance artifact
                decaytau = np.abs((max-steady50msavg)*0.368+steady50msavg-data[:,levelnum][maxtimeidx+990:5990]).argmin()/samprate
            levelnum +=1
            act_params_dict[f'{level:.2f}'+'_max']= max
            act_params_dict[f'{level:.2f}'+'_risetau']= risetau
            act_params_dict[f'{level:.2f}'+'_decaytau']= decaytau
        return act_params_dict

    def inact_params(self):
        print(self.ID)

    def ramp_params(self):
        print(self.ID)

    def AP_params(self):
        print(self.ID)

    def rec_params(self):
        print(self.ID)

first = KineticsCell(reader)
