# exec(open('Deepanjali data/PCA_tSNE.py').read())

import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects

import seaborn as sns
sns.set_style('darkgrid')
sns.set_palette('muted')
sns.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 2.5})
RS = 123

features_WT = pd.read_csv('Deepanjali data/WT step input cells/features_WT.csv',sep='\t')
features_KO = pd.read_csv('Deepanjali data/KO step input cells/features_FXS.csv',sep='\t')
features = pd.concat([features_WT,features_KO], ignore_index=True)

y_train = np.concatenate([np.ones(27),np.zeros(17)])
X_train = features.drop(['Unnamed: 0', 'Cell name', 'Sampling rate', 'stim_start','stim_end', 'Rinput', 'Cm'], axis=1)
X_train['Type'] = y_train
X_train = X_train.sample(frac=1).reset_index(drop=True)
y_train = X_train['Type']
X_train = X_train.drop(['Type'], axis=1)
X_train = X_train.dropna(axis=1)


# Utility function to visualize the outputs of PCA and t-SNE

def fashion_scatter(x, colors):
    # choose a color palette with seaborn.
    num_classes = len(np.unique(colors))
    palette = np.array(sns.color_palette("hls", num_classes))

    # create a scatter plot.
    f = plt.figure(figsize=(8, 8))
    ax = plt.subplot(aspect='equal')
    sc = ax.scatter(x[:,0], x[:,1], lw=0, s=40, c=palette[colors.astype(np.int)])
    plt.xlim(-25, 25)
    plt.ylim(-25, 25)
    ax.axis('off')
    ax.axis('tight')

    # add the labels for each digit corresponding to the label
    txts = []

    for i in range(num_classes):

        # Position of each label at median of data points.

        xtext, ytext = np.median(x[colors == i, :], axis=0)
        txt = ax.text(xtext, ytext, str(i), fontsize=24)
        txt.set_path_effects([
            PathEffects.Stroke(linewidth=5, foreground="w"),
            PathEffects.Normal()])
        txts.append(txt)

    return f, ax, sc, txts


x_subset = X_train
y_subset = y_train

print(np.unique(y_subset))

from sklearn.decomposition import PCA

time_start = time.time()

pca = PCA(n_components=6)
pca_result = pca.fit_transform(x_subset)


pca_df = pd.DataFrame(columns = ['pca1','pca2','pca3','pca4'])

pca_df['pca1'] = pca_result[:,0]
pca_df['pca2'] = pca_result[:,1]
pca_df['pca3'] = pca_result[:,2]
pca_df['pca4'] = pca_result[:,3]
pca_df['pca5'] = pca_result[:,4]
pca_df['pca6'] = pca_result[:,5]

print('Variance explained per principal component: {}'.format(pca.explained_variance_ratio_))


top_two_comp = pca_df[['pca1','pca2']] # taking first and second principal component

fashion_scatter(top_two_comp.values,y_subset) # Visualizing the PCA output
plt.show()


#### t-SNE ####
from sklearn.manifold import TSNE
import time
time_start = time.time()

fashion_tsne = TSNE(random_state=RS).fit_transform(x_subset)

fashion_scatter(fashion_tsne, y_subset)
plt.show()

## 3d plotting##
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(pca_df['pca1'][y_train==1], pca_df['pca2'][y_train==1], pca_df['pca3'][y_train==1])
ax.scatter(pca_df['pca1'][y_train==0], pca_df['pca2'][y_train==0], pca_df['pca3'][y_train==0])
plt.show()
