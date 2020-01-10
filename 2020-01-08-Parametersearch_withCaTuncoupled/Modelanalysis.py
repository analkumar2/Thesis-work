## exec(open('Modelanalysis.py').read())

import os
import numpy as np
import MOOSEModel_4 as mm
from pprint import pprint
import pandas as pd
import matplotlib.pyplot as plt

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
bestModel = {'Score': {'E_rest_0': -3.656096578281632, 'Input resistance': 2.618358563498228, 'Cell capacitance': -1.7160117451340426, 'AP1_amp_1.5e-10': 3.0073474446191715, 'APp_amp_1.5e-10': 4.28262079037893, 'AP1_width_1.5e-10': 3.546910013201938, 'APp_width_1.5e-10': 1.4433108636905843, 'AP1_thresh_1.5e-10': 8.069615517720345, 'APp_thresh_1.5e-10': 4.185651154692111, 'AP1_lat_1.5e-10': -100.59041747250744, 'ISI1_1.5e-10': 7.053352164588987, 'ISIl_1.5e-10': 0.8692310824929669, 'ISIavg_1.5e-10': 3.486661451037165, 'freq_1.5e-10': 4.8424949723257225, 'Adptn_id_1.5e-10': -5.258617151824255, 'fAHP_AP1_amp_1.5e-10': -5.572340116707596, 'fAHP_APp_amp_1.5e-10': -6.796472779709998, 'mAHP_stimend_amp_1.5e-10': -6.433196215719613, 'sAHP_stimend_amp_1.5e-10': -9.035656538623803, 'AHP_AP1_amp_1.5e-10': -5.843274022334083, 'AHP_APp_amp_1.5e-10': -7.679046600479444, 'AHP_AP1_time_1.5e-10': 1.1216178165245565, 'AHP_APp_time_1.5e-10': -1.0980387798168443}, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.6976647596774438e-10, 'Rm': 445660177.7570957, 'Em': -0.06083873307789175}, 'Channels': {'K_BK_Chan': {'Gbar': 4.8547260081837e-09, 'Kinetics': '../../Compilations/Kinetics/K_BK_Chan_(Migliore2018)', 'Erev': -0.08784832728715974}, 'Ca_T_Chan': {'Gbar': 5.870570859348553e-10, 'Kinetics': '../../Compilations/Kinetics/Ca_T_Chan_(Migliore2018)', 'Erev': 0.12903466644471448}, 'Ca_L_Chan': {'Gbar': 9.197015186069505e-09, 'Kinetics': '../../Compilations/Kinetics/Ca_L_Chan_(Migliore2018)', 'Erev': 0.12903466644471448}, 'Ca_N_Chan': {'Gbar': 3.110487249083809e-09, 'Kinetics': '../../Compilations/Kinetics/Ca_N_Chan_(Migliore2018)', 'Erev': 0.12903466644471448}, 'Na_Chan': {'Gbar': 1.8417242546181758e-06, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.09815969870703357, 'gateX': [-0.024839180357274598, 0.006840872243896084, -0.02625973279336557, 0.005816100913287417, 0.061806304205302226, 0.06434498333539079, 0.04110188846240265, 0.0007202810242443193], 'gateY': [-0.04669752782154454, -0.0038254098311545325, -0.04464387570075286, 0.0173863901291127, 0.0709166306215429, 0.07689728407370325, 0.02574520703900539, 0.024125398501901638]}, 'Na_P_Chan': {'Gbar': 1.0148518112707323e-08, 'Kinetics': '../../Compilations/Kinetics/Na_P_Chan_(Migliore2018)', 'Erev': 0.09815969870703357}, 'K_DR_Chan': {'Gbar': 1.3318413631216184e-05, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.08784832728715974, 'gateX': [0.02031655320492873, 0.005113587950698997, -0.025998968840198908, 0.02886546749687531, 0.06733486234092359, 0.0032422979293821653, 0.04433434609999344, 0.10827707744692375]}, 'K_D_Chan': {'Gbar': 4.121786813786762e-10, 'Kinetics': '../../Compilations/Kinetics/K_D_Chan_(Migliore2018)', 'Erev': -0.08784832728715974}, 'K_A_Chan': {'Gbar': 2.657957471910208e-06, 'Kinetics': '../../Compilations/Kinetics/K_A_Chan_(Migliore2018)_ghk', 'Erev': -0.08784832728715974}, 'K_M_Chan': {'Gbar': 4.1658868047006385e-08, 'Kinetics': '../../Compilations/Kinetics/K_M_Chan_(Migliore2018)', 'Erev': -0.08784832728715974}, 'K_SK_Chan': {'Gbar': 1.281494894541936e-07, 'Kinetics': '../../Compilations/Kinetics/K_SK_Chan_(Migliore2018)', 'Erev': -0.08784832728715974}, 'h_Chan': {'Gbar': 1.3297186404470057e-08, 'Kinetics': '../../Compilations/Kinetics/h_Chan_(Migliore2018)', 'Erev': -0.032676131451251055}}, 'Ca_Conc': {'Ca_B': 146607090524.42236, 'Ca_tau': 0.0023588604208249846, 'Ca_base': 5.405348001826576e-05, 'Kinetics': '../../Compilations/Kinetics/Ca_Conc_(Common)'}}}
bestModel['Score'] = absDict(bestModel['Score'])

nan = 1000
mnum = 0
Modelslist = []
errors = []
for file in os.listdir('OutputModels1'):
	exec(open('OutputModels1/'+file).read())
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

mm.plotModel(Modelslist[0], 150e-12)