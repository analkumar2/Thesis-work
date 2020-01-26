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

def absDict(Dict):
	for key in Dict:
		Dict[key] = np.abs(Dict[key])
	return Dict

def domDict(Dict1, Dict2):
	truthtable = []
	for key in Dict1['Score']:
		truthtable.append(Dict1['Score'][key]>=Dict2['Score'][key])
	if all(truthtable):
		return Dict2
	else:
		return Dict1

def betterDict(Dict1, Dict2):
	if np.sum(list(Dict1['Score'].values()))> np.sum(list(Dict2['Score'].values())):
		return Dict2
	else:
		return Dict1


nan = 1000
mnum = 0
Modelslist = []
errors = []
for file in os.listdir('../../Output/2020-01-08-Parametersearch_withNaKDRKinetics/Output'):
	exec(open('../../Output/2020-01-08-Parametersearch_withNaKDRKinetics/Output/'+file).read())
	for model in Models:
		mnum = mnum+1
		Models[model]['Score'] = absDict(Models[model]['Score'])
		errors.append(np.sum(list(Models[model]['Score'].values())))
		Modelslist.append(Models[model])
		# bestModel = betterDict(bestModel, Models[model])
		pprint(file+'_'+model)
		# pprint(Models[model]['Score'])
		# mm.plotModel(Models[model], 150e-12)

sortidx = np.argsort(errors)
Modelslist = np.array(Modelslist)
errors = np.array(errors)
Modelslist = Modelslist[sortidx]
errors = errors[sortidx]

# rdes = mm.plotModel(Modelslist[0], 150e-12)


## The original bounds of the free parameters
sm_area = 14.2e-9
bounds_freeparams = {}
bounds_freeparams['Cm'] = np.array([80e-12, 200e-12])
bounds_freeparams['Rm'] = np.array([150e6, 1000e6])
bounds_freeparams['Em'] = np.array([-0.1, -0.05])
bounds_freeparams['Ca_B'] = np.array([1.8e7, 1.8e11])
bounds_freeparams['Ca_tau'] = np.array([0.001, 0.200])
bounds_freeparams['Ca_base'] = np.array([0.01e-3, 1e-3])
bounds_freeparams['K_BK_Chan_Gbar'] = np.array([0.00001, 8])*sm_area
bounds_freeparams['Ca_T_Chan_Gbar'] = np.array([0.00001, 1])*sm_area
bounds_freeparams['Ca_L_Chan_Gbar'] = np.array([0.00001, 1])*sm_area
bounds_freeparams['Ca_N_Chan_Gbar'] = np.array([0.00001, 1])*sm_area
bounds_freeparams['Na_Chan_Gbar'] = np.array([7, 1000])*sm_area
bounds_freeparams['Na_P_Chan_Gbar'] = np.array([0.0001, 1])*sm_area
bounds_freeparams['K_DR_Chan_Gbar'] = np.array([0.00001, 1000])*sm_area
bounds_freeparams['K_D_Chan_Gbar'] = np.array([0.00001, 0.5])*sm_area
bounds_freeparams['K_A_Chan_Gbar'] = np.array([0.00001, 1000])*sm_area
bounds_freeparams['K_M_Chan_Gbar'] = np.array([0.00001, 11])*sm_area
bounds_freeparams['K_SK_Chan_Gbar'] = np.array([0.00001, 11])*sm_area
bounds_freeparams['h_Chan_Gbar'] = np.array([0.00001, 2.5])*sm_area
# bounds_freeparams
# bounds_freeparams
bounds_freeparams['Ek'] = np.array([-0.100, -0.080])
bounds_freeparams['Eh'] = np.array([-0.050, -0.030])
bounds_freeparams['ENa'] = np.array([0.050, 0.100])
bounds_freeparams['ECa'] = np.array([0.120, 0.140])
bounds_freeparams['Na_Chan_X_vhalf'] = np.array([-0.038,-0.024])
bounds_freeparams['Na_Chan_X_slope'] = np.array([0.002,0.015])
bounds_freeparams['Na_Chan_X_A'] = np.array([-0.050,-0.017])
bounds_freeparams['Na_Chan_X_B'] = np.array([0.002,0.040])
bounds_freeparams['Na_Chan_X_C'] = np.array([0,0.1])
bounds_freeparams['Na_Chan_X_D'] = np.array([0,0.1])
bounds_freeparams['Na_Chan_X_E'] = np.array([0.010,0.050])
bounds_freeparams['Na_Chan_X_F'] = np.array([0.0004,0.002])
bounds_freeparams['Na_Chan_Y_vhalf'] = np.array([-0.050,-0.030])
bounds_freeparams['Na_Chan_Y_slope'] = np.array([-0.010,-0.002])
bounds_freeparams['Na_Chan_Y_A'] = np.array([-0.060,-0.035])
bounds_freeparams['Na_Chan_Y_B'] = np.array([0.002,0.020])
bounds_freeparams['Na_Chan_Y_C'] = np.array([0,0.1])
bounds_freeparams['Na_Chan_Y_D'] = np.array([0,0.1])
bounds_freeparams['Na_Chan_Y_E'] = np.array([0.002,0.030])
bounds_freeparams['Na_Chan_Y_F'] = np.array([0.002,0.07])
bounds_freeparams['K_DR_Chan_X_vhalf'] = np.array([-0.012,0.030])
bounds_freeparams['K_DR_Chan_X_slope'] = np.array([0.002,0.012])
bounds_freeparams['K_DR_Chan_X_A'] = np.array([-0.030,-0.008])
bounds_freeparams['K_DR_Chan_X_B'] = np.array([0.003,0.045])
bounds_freeparams['K_DR_Chan_X_C'] = np.array([0,0.1])
bounds_freeparams['K_DR_Chan_X_D'] = np.array([0,0.1])
bounds_freeparams['K_DR_Chan_X_E'] = np.array([0.005,0.055])
bounds_freeparams['K_DR_Chan_X_F'] = np.array([0.002,0.11])

def Modelslist_to_dict(Modelslist):
	Models_df = pd.DataFrame(columns=bounds_freeparams.keys())
	for i in np.arange(len(Modelslist)):
		print(i, end='\r')
		mm = Modelslist[i]
		passivee = mm['parameters']['Passive']
		Caa = mm['parameters']['Ca_Conc']
		kkK_BK = mm['parameters']['Channels']['K_BK_Chan']
		kkCa_T = mm['parameters']['Channels']['Ca_T_Chan']
		kkCa_L = mm['parameters']['Channels']['Ca_L_Chan']
		kkCa_N = mm['parameters']['Channels']['Ca_N_Chan']
		kkNa = mm['parameters']['Channels']['Na_Chan']
		kkNa_P = mm['parameters']['Channels']['Na_P_Chan']
		kkK_DR = mm['parameters']['Channels']['K_DR_Chan']
		kkK_D = mm['parameters']['Channels']['K_D_Chan']
		kkK_A = mm['parameters']['Channels']['K_A_Chan']
		kkK_M = mm['parameters']['Channels']['K_M_Chan']
		kkK_SK = mm['parameters']['Channels']['K_SK_Chan']
		kkh = mm['parameters']['Channels']['h_Chan']

		curr = {'Modelname':i, 'Cm':passivee['Cm'], 'Rm':passivee['Rm'], 'Em':passivee['Em'],
			'Ca_B':Caa['Ca_B'], 'Ca_tau':Caa['Ca_tau'], 'Ca_base':Caa['Ca_base'], 
			'K_BK_Chan_Gbar':kkK_BK['Gbar'], 'Ca_T_Chan_Gbar':kkCa_T['Gbar'], 'Ca_L_Chan_Gbar':kkCa_L['Gbar'], 'Ca_N_Chan_Gbar':kkCa_N['Gbar'], 'Na_Chan_Gbar':kkNa['Gbar'], 'Na_P_Chan_Gbar':kkNa_P['Gbar'], 'K_DR_Chan_Gbar':kkK_DR['Gbar'], 'K_D_Chan_Gbar':kkK_D['Gbar'], 'K_A_Chan_Gbar':kkK_A['Gbar'], 'K_M_Chan_Gbar':kkK_M['Gbar'], 'K_SK_Chan_Gbar':kkK_SK['Gbar'], 'h_Chan_Gbar':kkh['Gbar'],
			'Ek':kkK_DR['Erev'], 'Eh':kkh['Erev'], 'ENa':kkNa['Erev'], 'ECa':kkCa_T['Erev'],
			'Na_Chan_X_vhalf':kkNa['gateX'][0], 'Na_Chan_X_slope':kkNa['gateX'][1], 'Na_Chan_X_A':kkNa['gateX'][2], 'Na_Chan_X_B':kkNa['gateX'][3], 'Na_Chan_X_C':kkNa['gateX'][4], 'Na_Chan_X_D':kkNa['gateX'][5], 'Na_Chan_X_E':kkNa['gateX'][6], 'Na_Chan_X_F':kkNa['gateX'][7],
			'Na_Chan_Y_vhalf':kkNa['gateY'][0], 'Na_Chan_Y_slope':kkNa['gateY'][1], 'Na_Chan_Y_A':kkNa['gateY'][2], 'Na_Chan_Y_B':kkNa['gateY'][3], 'Na_Chan_Y_C':kkNa['gateY'][4], 'Na_Chan_Y_D':kkNa['gateY'][5], 'Na_Chan_Y_E':kkNa['gateY'][6], 'Na_Chan_Y_F':kkNa['gateY'][7],
			'K_DR_Chan_X_vhalf':kkK_DR['gateX'][0], 'K_DR_Chan_X_slope':kkK_DR['gateX'][1], 'K_DR_Chan_X_A':kkK_DR['gateX'][2], 'K_DR_Chan_X_B':kkK_DR['gateX'][3], 'K_DR_Chan_X_C':kkK_DR['gateX'][4], 'K_DR_Chan_X_D':kkK_DR['gateX'][5], 'K_DR_Chan_X_E':kkK_DR['gateX'][6], 'K_DR_Chan_X_F':kkK_DR['gateX'][7],}
		Models_df = Models_df.append(curr, ignore_index=True)
	Models_df = Models_df.set_index('Modelname')
	return Models_df

Models_df = Modelslist_to_dict(Modelslist)


# ### uniformity score####
# nrm_var = {}
# for keyy in bounds_freeparams.keys():
# 	nrm_var[keyy] = Models_df.var()[keyy]/((bounds_freeparams[keyy][1] - bounds_freeparams[keyy][0])**2/12)

# pprint(nrm_var)

# for ppp in ['Na_Chan_X_vhalf', 'Na_Chan_Y_vhalf', 'Na_Chan_X_F', 'Na_Chan_Y_F', 'K_DR_Chan_X_vhalf', 'K_DR_Chan_X_F']:
# 	plt.hist(Models_df[ppp])
# 	plt.locator_params(axis='x', nbins=5)
# 	plt.xlabel('Value of the parameter')
# 	plt.ylabel('Number of models')
# 	plt.title(ppp)
# 	# plt.show()
# 	# plt.savefig('/home/analkumar2/Thesis_work/Codes/2019-12-27-justNaKDRkinetics_analysis/'+ppp)
# 	plt.savefig(ppp)
# 	plt.clf()

# ## Analysis all
# for ppp in bounds_freeparams.keys():
# 	plt.hist(Models_df[ppp], range=(bounds_freeparams[ppp][0], bounds_freeparams[ppp][1]))
# 	plt.locator_params(axis='x', nbins=5)
# 	plt.xlabel('Value of the parameter')
# 	plt.ylabel('Number of models')
# 	plt.title(ppp)
# 	# plt.show()
# 	# plt.savefig('/home/analkumar2/Thesis_work/Codes/2019-12-27-justNaKDRkinetics_analysis/'+ppp)
# 	plt.savefig(ppp)
# 	plt.clf()

## Analysis only for AHP score<3
ModelslistAHPl3 = []
for i in range(len(Modelslist)):
	if Modelslist[i]['Score']['AHP_APp_amp_1.5e-10']<3:
		ModelslistAHPl3.append(Modelslist[i])

Models_dfAHPl3 = Modelslist_to_dict(ModelslistAHPl3)
for ppp in bounds_freeparams.keys():
	plt.hist(Models_dfAHPl3[ppp], range=(bounds_freeparams[ppp][0], bounds_freeparams[ppp][1]))
	plt.locator_params(axis='x', nbins=5)
	plt.xlabel('Value of the parameter')
	plt.ylabel('Number of models')
	plt.title(ppp)
	# plt.show()
	# plt.savefig('/home/analkumar2/Thesis_work/Codes/2019-12-27-justNaKDRkinetics_analysis/'+ppp)
	plt.savefig('Output_AHPl3/'+ppp)
	plt.clf()

## Analysis only for AHP score<3
Modelslist_AP1_width_l3 = []
for i in range(len(Modelslist)):
	if Modelslist[i]['Score']['AP1_width_1.5e-10']<3:
		Modelslist_AP1_width_l3.append(Modelslist[i])

Models_df_AP1_width_l3 = Modelslist_to_dict(Modelslist_AP1_width_l3)
for ppp in bounds_freeparams.keys():
	plt.hist(Models_df_AP1_width_l3[ppp], range=(bounds_freeparams[ppp][0], bounds_freeparams[ppp][1]))
	plt.locator_params(axis='x', nbins=5)
	plt.xlabel('Value of the parameter')
	plt.ylabel('Number of models')
	plt.title(ppp)
	# plt.show()
	# plt.savefig('/home/analkumar2/Thesis_work/Codes/2019-12-27-justNaKDRkinetics_analysis/'+ppp)
	plt.savefig('Output_(AP1_width_1.5e-10)_l3/'+ppp)
	plt.clf()

# for ppp in ['Em', 'K_DR_Chan_X_F', 'K_DR_Chan_X_slope', 'Na_Chan_Gbar', 'Na_Chan_X_vhalf', 'Na_Chan_Y_F']:
# 	plt.hist(Models_df[ppp])
# 	plt.locator_params(axis='x', nbins=10)
# 	plt.xlabel('Value of the parameter')
# 	plt.ylabel('Number of models')
# 	plt.title(ppp)
# 	plt.show()

# pca = PCA()
# pca_result = pca.fit_transform(Models_df)



######## Correlation between Ca_base and Gbar of SK################
Ca_base_list = []
Gbar_SK_list = []
for i in range(len(Modelslist)):
	if i>10000:
		break
	print(i, end='\r')
	Ca_base_list.append(Modelslist[i]['parameters']['Ca_Conc']['Ca_base'])
	Gbar_SK_list.append(Modelslist[i]['parameters']['Channels']['K_SK_Chan']['Gbar'])

plt.scatter(Ca_base_list, Gbar_SK_list)

plt.xlabel('Ca_base (mM)')
plt.ylabel('Gbar SK channel (S)')
plt.show()
