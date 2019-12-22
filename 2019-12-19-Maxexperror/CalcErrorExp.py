#exec(open('CalcErrorExp.py').read())

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
import csv
import pandas as pd
import scipy.signal as scs
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
import plotexp


# //////////////////////////////////////////////////////////////////
foldername = '../../Raw_data/Deepanjali_data/WT step input cells' #Name of the folder where files are stored
fskip = ['Cell 2 of 21_3_2017.abf'] #List of files that need to be skipped
# //////////////////////////////////////////////////////////////////
mostdiffpair, error = [[0,0],0]
errors = []

numfiles = len(os.listdir(foldername))
for i in np.arange(numfiles):
    filename = os.listdir(foldername)[i]

    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf', 'Cell 2 of 19_10_2016', 'Cell 1 of 27_10_2016.abf', 'Cell 1 of 14_10_2016.abf', 'Cell 4 of 7_10_2016.abf', 'Cell 6 of 12_10_2016.abf', 'Cell 7 of 12_10_2016.abf']
    if filename in stim1391:
        stim_start = 139.1e-3
        stim_end = 639.1e-3
    else:
        stim_start = 81.4e-3 #in s
        stim_end = 581.4e-3 #in s

    if filename in fskip:
        print(f'{filename} skipped')
        continue

    if filename[-4:] == '.abf':
        curr_expData = {}
        tempT, tempV = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, 25e-12)
        curr_expData['25pA'] = [np.linspace(0, tempT[-1], int(tempT[-1]*10000)) , scs.resample(tempV, int(tempT[-1]*10000))]
        curr_expData['25pA'][0], curr_expData['25pA'][1] = curr_expData['25pA'][0][curr_expData['25pA'][0]<1] + (1-stim_start), curr_expData['25pA'][1][curr_expData['25pA'][0]<1]

        tempT, tempV = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, 50e-12)
        curr_expData['50pA'] = [np.linspace(0, tempT[-1], int(tempT[-1]*10000)) , scs.resample(tempV, int(tempT[-1]*10000))]
        curr_expData['50pA'][0], curr_expData['50pA'][1] = curr_expData['50pA'][0][curr_expData['50pA'][0]<1] + (1-stim_start), curr_expData['50pA'][1][curr_expData['50pA'][0]<1]

        tempT, tempV = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, 150e-12)
        curr_expData['150pA'] = [np.linspace(0, tempT[-1], int(tempT[-1]*10000)) , scs.resample(tempV, int(tempT[-1]*10000))]
        curr_expData['150pA'][0], curr_expData['150pA'][1] = curr_expData['150pA'][0][curr_expData['150pA'][0]<1] + (1-stim_start), curr_expData['150pA'][1][curr_expData['150pA'][0]<1]
    else:
        continue

    for k in np.arange(i+1, numfiles):
        filename = os.listdir(foldername)[k]

        stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf', 'Cell 2 of 19_10_2016', 'Cell 1 of 27_10_2016.abf', 'Cell 1 of 14_10_2016.abf', 'Cell 4 of 7_10_2016.abf', 'Cell 6 of 12_10_2016.abf', 'Cell 7 of 12_10_2016.abf']
        if filename in stim1391:
            stim_start = 139.1e-3
            stim_end = 639.1e-3
        else:
            stim_start = 81.4e-3 #in s
            stim_end = 581.4e-3 #in s

        if filename in fskip:
            print(f'{filename} skipped')
            continue

        if filename[-4:] == '.abf':
            next_expData = {}
            tempT, tempV = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, 25e-12)
            next_expData['25pA'] = [np.linspace(0, tempT[-1], int(tempT[-1]*10000)) , scs.resample(tempV, int(tempT[-1]*10000))]
            next_expData['25pA'][0], next_expData['25pA'][1] = next_expData['25pA'][0][next_expData['25pA'][0]<1] + (1-stim_start), next_expData['25pA'][1][next_expData['25pA'][0]<1]

            tempT, tempV = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, 50e-12)
            next_expData['50pA'] = [np.linspace(0, tempT[-1], int(tempT[-1]*10000)) , scs.resample(tempV, int(tempT[-1]*10000))]
            next_expData['50pA'][0], next_expData['50pA'][1] = next_expData['50pA'][0][next_expData['50pA'][0]<1] + (1-stim_start), next_expData['50pA'][1][next_expData['50pA'][0]<1]

            tempT, tempV = plotexp.expdata('../../Raw_data/Deepanjali_data/WT step input cells/'+filename, 150e-12)
            next_expData['150pA'] = [np.linspace(0, tempT[-1], int(tempT[-1]*10000)) , scs.resample(tempV, int(tempT[-1]*10000))]
            next_expData['150pA'][0], next_expData['150pA'][1] = next_expData['150pA'][0][next_expData['150pA'][0]<1] + (1-stim_start), next_expData['150pA'][1][next_expData['150pA'][0]<1]
        else:
            continue

        error25pA = np.sum((curr_expData['25pA'][1]-next_expData['25pA'][1])**2)
        error50pA = np.sum((curr_expData['50pA'][1]-next_expData['50pA'][1])**2)
        error150pA = np.sum((curr_expData['150pA'][1]-next_expData['150pA'][1])**2)
        currerror = error25pA+error50pA+error150pA
        print(i,',',k, '-->', currerror)
        errors.append(currerror)
        if  currerror > error:
            error  = currerror
            mostdiffpair = [i,k]
