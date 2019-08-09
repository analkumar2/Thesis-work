# exec(open('Deepanjali data/temp.py').read())

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
from neo.io import AxonIO
import os
import csv
import pandas as pd
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor


############ Was used to visualize each plot ##########
# Sno = 0
# for filename in os.listdir('Deepanjali data/WT step input cells'):
#     if Sno<19:
#         Sno = Sno +1
#         continue
#     if filename[-4:] == '.abf':
#         reader  = AxonIO(filename='Deepanjali data/WT step input cells/'+filename)
#     else:
#         continue
#     seg_no=16
#     Vtrace = reader.read_block().segments[seg_no].analogsignals[0]
#     v = np.array([float(V) for V in Vtrace])
#     t = np.linspace(0,float(Vtrace.t_stop - Vtrace.t_start), int((Vtrace.t_stop - Vtrace.t_start)*Vtrace.sampling_rate))
#     t = np.array(t)
#     plt.plot(t,v, '+', label=filename)
#     plt.legend()
#     Sno = Sno+1
#     if Sno>30:
#         break
#     print(Sno)
######################################################

############### Plotting feature correlation ############
features_df = pd.read_csv('Deepanjali Data/features.csv', delimiter='\t', index_col=0)
corr_features = np.abs(features_df.ix[:,4:].corr())


def plot_corr(corr_features):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cmap = cm.get_cmap('jet', 30)
    cax = ax1.imshow(corr_features, interpolation="nearest", cmap=cmap)
    ax1.grid(True)
    plt.title('Absolute Feature correlations')
    labels=list(corr_features.columns)
    ax1.set_xticks(range(46))
    ax1.set_yticks(range(46))
    ax1.set_xticklabels(labels,fontsize=6, rotation=90)
    ax1.set_yticklabels(labels,fontsize=6)
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    fig.colorbar(cax, ticks=[0,0.5,1])

plot_corr(corr_features)
plt.show(block=False)
