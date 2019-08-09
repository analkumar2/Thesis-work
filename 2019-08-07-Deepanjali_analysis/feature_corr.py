# exec(open('Deepanjali data/feature_corr.py').read())

# Plots a correlation graph but with absolute values of the correlations

from matplotlib import pyplot as plt
from matplotlib import cm as cm
import pandas as pd

# exec(open('Deepanjali data/Analysis_final.py').read())
features_pd = pd.read_csv("Deepanjali data/WT step input cells/features_WT.csv", sep = '\t', index_col = 0)
features_pd = features_pd.dropna(1) # drops features which are not present in any of the cells.

charchar = features_pd.iloc[:,np.arange(4,len(features_pd.columns))].copy()
# charchar = charchar.drop(columns=['ADP_AP1_amp_150.0','ADP_APp_amp_150.0','mAHP_AP1_dur_150.0', 'mAHP_AP1_amp_150.0'])

def correlation_plot(df):
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    cmap = cm.get_cmap('jet', 30)
    cax = ax1.imshow(np.abs(charchar.corr()), interpolation="nearest", cmap=cmap)
    ax1.grid(True)
    plt.title('Feature correlations')
    labels=list(charchar.columns)
    ax1.set_xticks(range(len(charchar.columns)))
    ax1.set_yticks(range(len(charchar.columns)))
    ax1.set_xticklabels(labels, rotation=90)
    ax1.set_yticklabels(labels)
    # Add colorbar, make sure to specify tick locations to match desired ticklabels
    fig.colorbar(cax)
    plt.show()

correlation_plot(charchar)

charchar_corr = np.abs(charchar.corr())
