# exec(open('Deepanjali data/visual_clustering.py').read())

# Compares and plots features between WT and FXS experimental recordings

import pandas as pd
import numpy as np


features_WT = pd.read_csv('Deepanjali data/WT step input cells/features_WT.csv',sep='\t')
features_KO = pd.read_csv('Deepanjali data/KO step input cells/features_FXS.csv',sep='\t')
features = pd.concat([features_WT,features_KO])

plt.close()
plt.clf()
for cc in range(5, len(features.columns)):
    print(features_KO.iloc[:,cc].name)
    plt.plot(np.zeros(len(features_KO)),features_KO.iloc[:,cc], '.')
    plt.plot(np.ones(len(features_WT)),features_WT.iloc[:,cc], '.')
    plt.xlim([-1,2])
    plt.title(features_KO.iloc[:,cc].name)
    plt.xticks([0,1], ['KO', 'WT'])
    plt.ylabel('Appropriate units')
    # plt.legend()
    plt.show()
    # plt.savefig(f'tempplots/{features_KO.iloc[:,cc].name}.png', bbox_inches='tight')
    plt.close()
    plt.clf()
