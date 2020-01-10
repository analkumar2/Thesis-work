## exec(open('featuresv2.py').read())

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import plotexp
import os
import csv
import pandas as pd
import brute_curvefit

from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

def expfeatures(cellpath,stim_start,stim_end):
    # tm25, Vtracem25 = plotexp.expdata(cellpath, -25e-12)
    t0, Vtrace0 = plotexp.expdata(cellpath, 0e-12)
    t25, Vtrace25 = plotexp.expdata(cellpath, 25e-12)
    t150, Vtrace150 = plotexp.expdata(cellpath, 150e-12)
    t300, Vtrace300 = plotexp.expdata(cellpath, 300e-12)
    features = {}
    features['Cell name'] = cellpath.split('/')[-1]
    features['Sampling rate'] = len(tm25)/tm25[-1]
    features['stim_start'] = stim_start
    features['stim_end'] = stim_end

    # features[f'E_rest_m25'] = np.median(Vtracem25[tm25<stim_start])
    features[f'E_rest_0'] = np.median(Vtrace0)
    features[f'E_rest_25'] = np.median(Vtrace25[t25<stim_start])
    features[f'E_rest_150'] = np.median(Vtrace150[t150<stim_start])
    features[f'E_rest_300'] = np.median(Vtrace300[t300<stim_start])

    def charging25(t, R,C):
        # C2 = 100
        return features[f'E_rest_25'] + R*25e-12*(1-np.exp(-t/R/C))

    # def chargingm25(t, R,C):
    #     return features[f'E_rest_m25'] - R*25e-12*(1-np.exp(-t/R/C))

    # def discharging25(t, R,C):
    #     return features[f'E_rest_25'] + R*25e-12*(np.exp(-t/R/C))

    # def dischargingm25(t, R,C):
    #     return features[f'E_rest_m25'] - R*25e-12*(np.exp(-t/R/C))

    ## 25
    # plt.plot(t25, Vtrace25)
    # plt.plot(t25[(t25>stim_start) & (t25<stim_start+0.1)], Vtrace25[(t25>stim_start) & (t25<stim_start+0.1)])
    # plt.show()

    ## charging25
    tempv = Vtrace25[(t25>stim_start) & (t25<stim_start+0.1)]
    RCfitted_ch25, error25 = brute_curvefit.brute_scifit(charging25, np.linspace(0,0.1,len(tempv)), tempv, restrict=[[50e6,50e-12],[250e6,250e-6]], ntol = 10000)
    # plt.plot(t25[(t25>stim_start) & (t25<stim_start+0.1)], charging25(np.linspace(0,0.1,len(tempv)), *RCfitted_ch25))

    # ## discharging25
    # tempv = Vtrace25[(t25>stim_end) & (t25<stim_end+0.1)]
    # RCfitted_disch25, error = brute_curvefit.brute_scifit(discharging25, np.linspace(0,0.1,len(tempv)), tempv, restrict=[[50e6,50e-12],[250e6,250e-6]], ntol = 10000)
    # plt.plot(t25[(t25>stim_end) & (t25<stim_end+0.1)], discharging25(np.linspace(0,0.1,len(tempv)), *RCfitted_disch25))
    # plt.show()

    ## m25
    # plt.plot(tm25, Vtracem25)
    # plt.plot(tm25[(tm25>stim_start) & (tm25<stim_start+0.1)], Vtracem25[(tm25>stim_start) & (tm25<stim_start+0.1)])
    # plt.show()

    ## charging25
    # tempv = Vtracem25[(tm25>stim_start) & (tm25<stim_start+0.1)]
    # RCfitted_chm25, errorm25 = brute_curvefit.brute_scifit(chargingm25, np.linspace(0,0.1,len(tempv)), tempv, restrict=[[50e6,50e-12],[250e6,250e-6]], ntol = 10000)
    # plt.plot(tm25[(tm25>stim_start) & (tm25<stim_start+0.1)], chargingm25(np.linspace(0,0.1,len(tempv)), *RCfitted_chm25))

    # ## discharging25
    # tempv = Vtracem25[(tm25>stim_end) & (tm25<stim_end+0.1)]
    # RCfitted_dischm25, error = brute_curvefit.brute_scifit(dischargingm25, np.linspace(0,0.1,len(tempv)), tempv, restrict=[[50e6,50e-12],[250e6,250e-6]], ntol = 10000)
    # plt.plot(tm25[(tm25>stim_end) & (tm25<stim_end+0.1)], dischargingm25(np.linspace(0,0.1,len(tempv)), *RCfitted_dischm25))
    # plt.show()
    features[f'Input Resistance'], features[f'Cell capacitance'] = RCfitted_ch25
    return features




if __name__ == '__main__':
    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf', 'Cell 2 of 19_10_2016', 'Cell 1 of 27_10_2016.abf', 'Cell 1 of 14_10_2016.abf', 'Cell 4 of 7_10_2016.abf', 'Cell 6 of 12_10_2016.abf', 'Cell 7 of 12_10_2016.abf']
    filename = 'cell 4 of 61016.abf'
    if filename in stim1391:
        stim_start = 139.1e-3
        stim_end = 639.1e-3
    else:
        stim_start = 81.4e-3 #in s
        stim_end = 581.4e-3 #in s
    charchar = expfeatures('Experimental recordings/'+filename, stim_start, stim_end)
