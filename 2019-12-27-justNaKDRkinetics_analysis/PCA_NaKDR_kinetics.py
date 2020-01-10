# exec(open('PCA_NaKDR_kinetics.py').read())

#First collect all the kinetics parameters
# Colloect all of their ranges in which the parameter was allowed to vary
# then for each parameter, calulate normalized variance = variance*12/(b-a)**2
# Do PCA and calculate normalized variance = variance*12/(b-a)**2 for each component
# IN both the cases the parameters with the least noramlized variance is our guy

import time
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import os
import seaborn as sbn
sbn.set()
import MOOSEModel_4 as mm
from pprint import *



Models_temp = {}
for file in os.listdir('repModels'):
	exec(open('Modelfiles/outputModels_dict_'+file.split('_')[0]+'.py').read())
	Models_temp[file.split('.')[0]] = Models[file.split('_')[1].split('.')[0]]

# rdes = mm.plotModel(Models_temp['2383402058_Model34'], 150e-12)

## The original bounds of the free parameters
sm_area = 14.2e-9
bounds_freeparams = {}
bounds_freeparams['Cm'] = (80e-12, 200e-12)
bounds_freeparams['Rm'] = (150e6, 1000e6)
bounds_freeparams['Em'] = (-0.1, -0.02)
bounds_freeparams['Na_Chan_Gbar'] = (7*sm_area, 1000*sm_area)
bounds_freeparams['K_DR_Chan_Gbar'] = (0.00001*sm_area, 1000*sm_area)
# bounds_freeparams['K_A_Chan_Gbar'] = (0.00001*sm_area, 1000*sm_area)
bounds_freeparams['Ek'] = (-0.100, -0.080)
bounds_freeparams['ENa'] = (0.050, 0.100)
bounds_freeparams['Na_Chan_X_vhalf'] = (-0.050,-0.024)
bounds_freeparams['Na_Chan_X_slope'] = (0.002,0.015)
bounds_freeparams['Na_Chan_X_A'] = (-0.050,-0.017)
bounds_freeparams['Na_Chan_X_B'] = (0.002,0.040)
bounds_freeparams['Na_Chan_X_C'] = (0,0.1)
bounds_freeparams['Na_Chan_X_D'] = (0,0.1)
bounds_freeparams['Na_Chan_X_E'] = (0.010,0.050)
bounds_freeparams['Na_Chan_X_F'] = (0.0004,0.002)
bounds_freeparams['Na_Chan_Y_vhalf'] = (-0.050,-0.030)
bounds_freeparams['Na_Chan_Y_slope'] = (-0.010,-0.002)
bounds_freeparams['Na_Chan_Y_A'] = (-0.060,-0.035)
bounds_freeparams['Na_Chan_Y_B'] = (0.002,0.020)
bounds_freeparams['Na_Chan_Y_C'] = (0,0.1)
bounds_freeparams['Na_Chan_Y_D'] = (0,0.1)
bounds_freeparams['Na_Chan_Y_E'] = (0.002,0.030)
bounds_freeparams['Na_Chan_Y_F'] = (0.002,0.1)
bounds_freeparams['K_DR_Chan_X_vhalf'] = (-0.012,0.030)
bounds_freeparams['K_DR_Chan_X_slope'] = (0.002,0.020)
bounds_freeparams['K_DR_Chan_X_A'] = (-0.030,-0.008)
bounds_freeparams['K_DR_Chan_X_B'] = (0.003,0.045)
bounds_freeparams['K_DR_Chan_X_C'] = (0,0.1)
bounds_freeparams['K_DR_Chan_X_D'] = (0,0.1)
bounds_freeparams['K_DR_Chan_X_E'] = (0.005,0.055)
bounds_freeparams['K_DR_Chan_X_F'] = (0.002,0.4)

Models_df = pd.DataFrame(columns=bounds_freeparams.keys())

for mm in Models_temp.keys():
	gg = Models_temp[mm]['parameters']['Passive']
	kkNa = Models_temp[mm]['parameters']['Channels']['Na_Chan']
	kkK_DR = Models_temp[mm]['parameters']['Channels']['K_DR_Chan']
	# kkK_A = Models_temp[mm]['parameters']['Channels']['K_A_Chan']
	curr = {'Modelname':mm, 'Cm':gg['Cm'], 'Rm':gg['Rm'], 'Em':gg['Em'],
		'Na_Chan_Gbar':kkNa['Gbar'], 'K_DR_Chan_Gbar':kkK_DR['Gbar'],
		'Ek':kkK_DR['Erev'], 'ENa':kkNa['Erev'],
		'Na_Chan_X_vhalf':kkNa['gateX'][0], 'Na_Chan_X_slope':kkNa['gateX'][1], 'Na_Chan_X_A':kkNa['gateX'][2], 'Na_Chan_X_B':kkNa['gateX'][3], 'Na_Chan_X_C':kkNa['gateX'][4], 'Na_Chan_X_D':kkNa['gateX'][5], 'Na_Chan_X_E':kkNa['gateX'][6], 'Na_Chan_X_F':kkNa['gateX'][7],
		'Na_Chan_Y_vhalf':kkNa['gateY'][0], 'Na_Chan_Y_slope':kkNa['gateY'][1], 'Na_Chan_Y_A':kkNa['gateY'][2], 'Na_Chan_Y_B':kkNa['gateY'][3], 'Na_Chan_Y_C':kkNa['gateY'][4], 'Na_Chan_Y_D':kkNa['gateY'][5], 'Na_Chan_Y_E':kkNa['gateY'][6], 'Na_Chan_Y_F':kkNa['gateY'][7],
		'K_DR_Chan_X_vhalf':kkK_DR['gateX'][0], 'K_DR_Chan_X_slope':kkK_DR['gateX'][1], 'K_DR_Chan_X_A':kkK_DR['gateX'][2], 'K_DR_Chan_X_B':kkK_DR['gateX'][3], 'K_DR_Chan_X_C':kkK_DR['gateX'][4], 'K_DR_Chan_X_D':kkK_DR['gateX'][5], 'K_DR_Chan_X_E':kkK_DR['gateX'][6], 'K_DR_Chan_X_F':kkK_DR['gateX'][7],}
	Models_df = Models_df.append(curr, ignore_index=True)
Models_df = Models_df.set_index('Modelname')

nrm_var = {}
for keyy in bounds_freeparams.keys():
	nrm_var[keyy] = Models_df.var()[keyy]/((bounds_freeparams[keyy][1] - bounds_freeparams[keyy][0])**2/12)

pprint(nrm_var)

for ppp in ['Na_Chan_X_vhalf', 'Na_Chan_Y_vhalf', 'Na_Chan_X_F', 'Na_Chan_Y_F', 'K_DR_Chan_X_vhalf', 'K_DR_Chan_X_F']:
	plt.hist(Models_df[ppp])
	plt.locator_params(axis='x', nbins=5)
	plt.xlabel('Value of the parameter')
	plt.ylabel('Number of models')
	plt.title(ppp)
	# plt.show()
	# plt.savefig('/home/analkumar2/Thesis_work/Codes/2019-12-27-justNaKDRkinetics_analysis/'+ppp)
	plt.savefig(ppp)
	plt.clf()

for ppp in ['Em', 'K_DR_Chan_X_F', 'K_DR_Chan_X_slope', 'Na_Chan_Gbar', 'Na_Chan_X_vhalf', 'Na_Chan_Y_F']:
	plt.hist(Models_df[ppp])
	plt.locator_params(axis='x', nbins=10)
	plt.xlabel('Value of the parameter')
	plt.ylabel('Number of models')
	plt.title(ppp)
	plt.show()

# pca = PCA()
# pca_result = pca.fit_transform(Models_df)
