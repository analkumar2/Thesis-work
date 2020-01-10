## exec(open('Deepanjalidata_Analysis.py').read())

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import featuresv3 as fv3
import plotexp as pex
import os
from sklearn.decomposition import PCA


# //////////////////////////////////////////////////////////////////
# foldername = '../../Raw_data/Deepanjali_data/WT step input cells/' #Name of the folder where files are stored
foldername = '../../Raw_data/Deepanjali_data/KO step input cells/' #Name of the folder where files are stored
fskip = ['Cell 2 of 21_3_2017.abf'] #List of files that need to be skipped
# //////////////////////////////////////////////////////////////////

Sno = 0
features_pd = pd.DataFrame()
for filename in os.listdir(foldername):
    stim1391 = ['Cell 3 of 181016.abf', 'cell 4 of 61016.abf', 'cell 4 of 111016.abf', 'cell 4 of 131016.abf', 'Cell 4 of 181016.abf', 'cell 5 of 61016.abf', 'Cell 5 of 181016.abf', 'Cell 2 of 19_10_2016', 'Cell 1 of 27_10_2016.abf', 'Cell 1 of 14_10_2016.abf', 'Cell 4 of 7_10_2016.abf', 'Cell 6 of 12_10_2016.abf', 'Cell 7 of 12_10_2016.abf', 'Cell 2 of 19_10_2016.abf']
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
        print(filename)
        features = fv3.expfeatures(foldername+filename, stim_start, stim_end)
        features_pd = features_pd.append(pd.DataFrame(features,index = [Sno]))
        Sno = Sno +1
        print(Sno)
        plt.axvline(x=stim_start)
        plt.axvline(x=stim_end)
        plt.plot(*pex.expdata(foldername+filename, 150e-12))
        plt.plot(*pex.expdata(foldername+filename, 300e-12))
        plt.show()
    else:
        continue

# features_pd.to_csv('features.csv', sep='\t')
features_ess_pd= features_pd.dropna(axis='columns')
features_ess_pd= features_ess_pd.drop(['Sampling rate','stim_start','stim_end','E_rest_25','E_rest_150','E_rest_300'], axis='columns')
features_ess_pd.to_csv('features_ess_KO.csv', sep='\t')

### features_ess_WT = pd.read_csv('features_ess_WT.csv',sep='\t', index_col=0) #To load a csv to dataframe
