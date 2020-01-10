# exec(open('PCA_NaKDR_kinetics.py').read())

#First collect all the kinetics parameters
# Colloect all of their ranges in which the parameter was allowed to vary
# then for each parameter, calulate normalized variance = variance*12/(b-a)**2
# Do PCA and calculate normalized variance = variance*12/(b-a)**2 for each component
# IN both the cases the parameters with the least noramlized variance is our guy

import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import os
import MOOSEModel_4 as mm
from pprint import *



Models_temp = {}
for file in os.listdir('repModels'):
	exec(open('Modelfiles/outputModels_dict_'+file.split('_')[0]+'.py').read())
	Models_temp[file.split('.')[0]] = Models[file.split('_')[1].split('.')[0]]

rdes = mm.plotModel(Models_temp['2383402058_Model34'], 150e-12)