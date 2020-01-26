## exec(open('bestexprec.py').read())

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import featuresv6 as fts
import plotexpv1 as pex
import os

# Importing the means and std of WT
features_ess_WT = pd.read_csv('features_ess_WT.csv',sep='\t', index_col=0)
meanfeatures_WT = features_ess_WT.mean()
stdfeatures_WT = features_ess_WT.std()

features_score_WT = features_ess_WT
for key in meanfeatures_WT.keys():
	features_score_WT.loc[:,key] = np.abs((features_ess_WT.loc[:,key] - meanfeatures_WT[key])/stdfeatures_WT[key])


print(np.argsort(features_score_WT.sum(axis='columns')))
foldername = '../../Raw_data/Deepanjali_data/WT step input cells/'
# cellname = 'Cell 4 of 201117.abf' #best
cellname = 'cell 5 of 61016.abf' #second best. Finally used this one
# cellname = 'cell 4 of 61016.abf' #third best
# cellname = 'Cell 1 of 25517.abf' #fourth best
plt.plot(*pex.expdata(foldername+cellname, 150e-12), label='150pA')
plt.plot(*pex.expdata(foldername+cellname, 300e-12), label='300pA')
plt.title(f'best experimental WT - {cellname}')
plt.legend()
plt.show()



# Importing the means and std of KO
features_ess_KO = pd.read_csv('features_ess_KO.csv',sep='\t', index_col=0)
meanfeatures_KO = features_ess_KO.mean()
stdfeatures_KO = features_ess_KO.std()

features_score_KO = features_ess_KO
for key in meanfeatures_KO.keys():
	features_score_KO.loc[:,key] = np.abs((features_ess_KO.loc[:,key] - meanfeatures_KO[key])/stdfeatures_KO[key])

print(np.argsort(features_score_KO.sum(axis='columns')))
foldername = '../../Raw_data/Deepanjali_data/KO step input cells/'
cellname = 'Cell 6 of 12_10_2016.abf' #best. Finally use this one
# cellname = 'Cell 3 of 21_3_2017.abf' #second best
# cellname = 'Cell 2 of 19_10_2016.abf' #third best
plt.plot(*pex.expdata(foldername+cellname, 150e-12), label='150pA')
plt.plot(*pex.expdata(foldername+cellname, 300e-12), label='300pA')
plt.title(f'best experimental KO - {cellname}')
plt.legend()
plt.show()


#####Score analysis####
features_score_WT = features_score_WT.set_index('Cell name')
features_score_KO = features_score_KO.set_index('Cell name')

print(features_score_WT[features_score_WT<1.1].dropna())
print(features_score_KO[features_score_KO<1.1].dropna())