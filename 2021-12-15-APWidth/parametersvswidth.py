## exec(open('parametercrawl_v6.py').read())
## parallel -u 'sleep {1}; for i in $(seq 2); do echo $i; time python3 parametercrawl_v7.py; done' ::: $(seq 3)
## randomly selects a free parameter value and then stores the model with the scores
## THe features are aslo stored
## Samples log-space instead of uniform
## THis version removes the unnecessary channels altogether
## v2: Added K_A_Chan. THus the total channels now become - Na, K_DR, K_A, K_SK, Ca_T
## v3: SK kinetics is now a free parameter
## v4: SK kinetics is stored also now.
## v5: Uses custom versions of Na, KDR and KA. Dummy model specification system changed
## v6: Simplified
## v7: Optimized for speed - writes both WT and KO files
## v8: Optimized for speed - AUtomatically sets up RM and Em

from neo.io import AxonIO
import time
import pandas as pd
import plotexpv3 as pex
import MOOSEModel_18 as mm

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from copy import deepcopy
import numpy.random as nr
import featuresv29_nonallen as fts
from pprint import pprint
import scipy.signal as scs
import sys
import os
import io
import scipy.signal as scs
import scipy
import scipy.interpolate as sciip
from pprint import pprint
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

# Setting seed
seeed = nr.randint(2 ** 32 - 1)
# seeed = 110043517
nr.seed(seeed)
print(seeed)
f = open("OutputModels/outputModels_dict_"+str(seeed)+".py","a+")
f.write("# exec(open('OutputModels/outputModels_dict_"+str(seeed)+".py').read())\n\n")
f.write("Models = {} \n\n")


# Importing the means and std
features_ess_WT = pd.read_csv("features_ess_WT_v3.csv", sep="\t", index_col=0)
meanfeaturesWT = features_ess_WT.mean()
stdfeaturesWT = features_ess_WT.std()
minfeaturesWT = features_ess_WT.min()
maxfeaturesWT = features_ess_WT.max()

features_ess_KO = pd.read_csv("features_ess_KO_v3.csv", sep="\t", index_col=0)
meanfeaturesKO = features_ess_KO.mean()
stdfeaturesKO = features_ess_KO.std()
minfeaturesKO = features_ess_KO.min()
maxfeaturesKO = features_ess_KO.max()


#dummy Model
sm_area = 14.2e-9
dummyModel = {
    "Parameters": {
        "notes": "",
        "Morphology": {"sm_len": 6.73077545020806e-05, "sm_diam": 6.73077545020806e-05},
        "Passive": {
            "Cm": 1.19e-10,
            "Rm": 160496106.7347296,
            "Em": -0.05071961460259502,
        },
        "Channels": {
            "Na_Chan": {
                "Gbar": 0.00012480475550619003*1.5,
                "Erev": 0.06,
                "Kinetics": "../../Compilations/Kinetics/Na_Chan_Custom4",
                "KineticVars": {
                    "m_vhalf_inf":-31.6e-3, "m_slope_inf":6.8e-3, "m_A":-3.65E-02, "m_B":2.00E-02, "m_C":1.61E-02, "m_D":5.47E-02, "m_E":3.11E-02, "m_F":6.40E-04,
                    "h_vhalf_inf":-66e-3, "h_slope_inf":-5.3e-3, "h_A":-0.04560699, "h_B":0.00433522, "h_C":0.01197575, "h_D":0.02617791, "h_E":0.00853832, "h_F":0.03900321,
                    # "s_vhalf_inf":-45e-3, "s_slope_inf":-6e-3, "s_A":1, "s_B":0.001, "s_C":0.001, "s_D":0.5, "s_E":0.001, "s_F":1,
                    # "m_vhalf_inf": -0.0316,
                    # "h_vhalf_inf": -0.066,
                    "s_vhalf_inf": -0.033,
                },
            },
            "K_DR_Chan": {
                "Gbar": 1.0502259538910637e-7,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_DR_Chan_Custom3",
                "KineticVars": {
                    "n_vhalf_inf":0.013, "n_slope_inf":0.0087666, "n_A":1.26E-02, "n_B":1.73E-02, "n_C":0, "n_D":0, "n_E":3.43E-02, "n_F":1.02E-01,
                    "n_F": 0.00306,
                    # "n_vhalf_inf": 0.0112,
                    # "n_slope_inf": 0.017,
                    # "n_A": -8.78e-3,
                    # "n_B": 5.63e-2,
                    # "n_C": 0,
                    # "n_D": 0,
                    # "n_E": 2.65e-2,
                    # "n_F": 1.05e-2,
                },
            },
            "K_A_Chan": {
                "Gbar": 1.008422244061249e-06,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_A_Chan_Custom3",
                "KineticVars": {
                    "n_vhalf_inf":0.0112, "n_slope_inf":0.017, "n_A":-8.78e-3, "n_B":5.63e-2, "n_C":0, "n_D":0, "n_E":2.65e-2, "n_F":1.05e-2,
                    "l_vhalf_inf":-0.056, "l_slope_inf":-0.00877, "l_min":0.002, "l_m":0.26, "l_cm":0.050,
                    "n_vhalf_inf": 0.025,
                },
            },
            "K_M_Chan": {
                "Gbar": 4.016076778584836e-09*2,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_M_Chan_Custom1",
                "KineticVars": {"factor": 3.3e-05},
            },
            "h_Chan": {
                "Gbar": 5.3739087243907273e-11,
                "Erev": -0.04,
                "Kinetics": "../../Compilations/Kinetics/h_Chan_Custom1",
                "KineticVars": {},
            },
            # "Ca_L_Chan": {
            #     "Gbar": 1.0330355973445023e-10,
            #     "Erev": 0.14,
            #     "Kinetics": "../../Compilations/Kinetics/Ca_L_Chan_Custom1",
            # },
            # "K_SK_Chan": {
            #     "Gbar": 2.0605218247275846e-09,
            #     "Erev": -0.09,
            #     "Kinetics": "../../Compilations/Kinetics/K_SK_Chan_Custom4",
            # },
        },
        "Ca_Conc": {
            "Ca_B": 75427936887.46373,
            "Ca_tau": 0.038,
            "Ca_base": 8e-05,
            "Kinetics": "../../Compilations/Kinetics/Ca_Conc_(Common)",
        },
    },
}

# # Load the passive models. Em, Rm, Cm
# exec(open("properpasModels_l3.py").read())

# Define constants not to be changed
F = 96485.3329
elecPlotDt = 0.00005
elecDt = 0.00005

# Generating base model
# sys.stdout = open(os.devnull, "w")
rdes = mm.generateModel(dummyModel, 150e-12)  # Dummy model to load the channel kinetics
rdes.buildModel()
# sys.stdout = sys.__stdout__

# Actual for loop
num_of_iterations = 2500
mnum = 1
for channame in list(dummyModel["Parameters"]["Channels"].keys()):
    for param in list(dummyModel["Parameters"]["Channels"][channame]["KineticVars"].keys()):
        ttt = time.time()
        # sys.stdout = open(os.devnull, "w")
        Modeltemp = deepcopy(dummyModel)
        Modeltemp["Scores"] = {}
        # Modeltemp["Parameters"] = {}
        Modeltemp["Parameters"]["notes"] = ""
        # if any(xx in channame for xx in ["vhalf", "A", "C", 'D']):
        #     Modeltemp["Parameters"]["Channels"][channame]["KineticVars"][param] = Modeltemp["Parameters"]["Channels"][channame]["KineticVars"][param] - 0.001
        # else:
        #     Modeltemp["Parameters"]["Channels"][channame]["KineticVars"][param] = Modeltemp["Parameters"]["Channels"][channame]["KineticVars"][param] * 0.99

        W_Erest = nr.uniform(minfeaturesWT['E_rest_0'], maxfeaturesWT['E_rest_0'])
        W_Rin = nr.uniform(minfeaturesWT['Input resistance'], maxfeaturesWT['Input resistance']+20e6)
        W_Gtot = 1/W_Rin
        soma = moose.element('/model/elec/soma')
        sf_area = np.pi*soma.length*soma.diameter

        I_throughChan_sum = 0
        A_Gtot = 0
        for chan in moose.wildcardFind("/library/##[CLASS==HHChannel]"):
            if chan.useConcentration == 1:
                Gbar = Modeltemp["Parameters"]["Channels"][chan.name]["Gbar"]
                gateX = moose.element(f"{chan.path}/gateX")
                gateY = moose.element(f"{chan.path}/gateY")
                gateZ = moose.element(f"{chan.path}/gateZ")
                G = Gbar
                if chan.Xpower > 0:
                    DIVSX = np.linspace(gateX.min, gateX.max, len(gateX.tableA))
                    finterpX = scipy.interpolate.interp1d(DIVSX, gateX.tableA / gateX.tableB)
                    gateinfX = finterpX(W_Ca)
                    G = G * gateinfX ** chan.Xpower
                if chan.Ypower > 0:
                    DIVSY = np.linspace(gateY.min, gateY.max, len(gateY.tableA))
                    finterpY = scipy.interpolate.interp1d(DIVSY, gateY.tableA / gateY.tableB)
                    gateinfY = finterpY(W_Ca)
                    G = G * gateinfY ** chan.Ypower
                if chan.Zpower > 0:
                    DIVSZ = np.linspace(gateZ.min, gateZ.max, len(gateZ.tableA))
                    finterpZ = scipy.interpolate.interp1d(DIVSZ, gateZ.tableA / gateZ.tableB)
                    gateinfZ = finterpZ(W_Ca)
                    G = G * gateinfZ ** chan.Zpower
                I_throughChan = G*(W_Erest - Modeltemp["Parameters"]["Channels"][chan.name]["Erev"])
                I_throughChan_sum = I_throughChan_sum + I_throughChan
                A_Gtot = A_Gtot + G
            else:
                Gbar = Modeltemp["Parameters"]["Channels"][chan.name]["Gbar"]
                gateX = moose.element(f"{chan.path}/gateX")
                gateY = moose.element(f"{chan.path}/gateY")
                gateZ = moose.element(f"{chan.path}/gateZ")
                G = Gbar
                if chan.Xpower > 0:
                    DIVSX = np.linspace(gateX.min, gateX.max, len(gateX.tableA))
                    finterpX = scipy.interpolate.interp1d(DIVSX, gateX.tableA / gateX.tableB)
                    gateinfX = finterpX(W_Erest)
                    G = G * gateinfX ** chan.Xpower
                if chan.Ypower > 0:
                    DIVSY = np.linspace(gateY.min, gateY.max, len(gateY.tableA))
                    finterpY = scipy.interpolate.interp1d(DIVSY, gateY.tableA / gateY.tableB)
                    gateinfY = finterpY(W_Erest)
                    G = G * gateinfY ** chan.Ypower
                if chan.Zpower > 0:
                    DIVSZ = np.linspace(gateZ.min, gateZ.max, len(gateZ.tableA))
                    finterpZ = scipy.interpolate.interp1d(DIVSZ, gateZ.tableA / gateZ.tableB)
                    gateinfZ = finterpZ(W_Erest)
                    G = G * gateinfZ ** chan.Zpower
                I_throughChan = G*(W_Erest - Modeltemp["Parameters"]["Channels"][chan.name]["Erev"])
                I_throughChan_sum = I_throughChan_sum + I_throughChan
                A_Gtot = A_Gtot + G
            # print(chan, G)

        Gl = W_Gtot - A_Gtot
        if Gl <=0:
            # sys.stdout = sys.__stdout__
            print(seeed, i, Gl)
            continue

        Modeltemp["Parameters"]["Passive"]["Rm"] = 1/Gl
        Modeltemp["Parameters"]["Passive"]["Em"] = (I_throughChan_sum + Gl*W_Erest)/Gl

        ######################################################################
        tvec, Vmvec, Cavec = mm.runModel(Modeltemp, CurrInjection=150e-12)

        tt = tvec
        vv = Vmvec
        ii = np.zeros(len(tt))
        ii[(tt >= 1) & (tt <= 1.5)] = 150e-12
        sweep_ext = EphysSweepFeatureExtractor(
            t=tt, v=vv * 1e3, i=ii * 1e12, filter=len(tt) / tt[-1] / 2500
        )
        sweep_ext.process_spikes()
        p1 = sweep_ext.spike_feature('peak_index')[0]
        f2 = sciip.interp1d(vv[p1-10:p1+1], tt[p1-10:p1+1])
        f3 = sciip.interp1d(vv[p1:p1+30+1], tt[p1:p1+30+1])
        AP1_width = f3(0) - f2(0)

        print(channame,param, Modeltemp["Parameters"]["Channels"][channame]["KineticVars"][param], f'{AP1_width =}')
        plt.plot(tvec, Vmvec)
        plt.show()

        # peaks = scs.find_peaks(x=Vmvec, height=None, threshold=None, distance=None, prominence=0.050, width=None, wlen=None, rel_height=0.5, plateau_size=None)
        # if len(peaks[0])<4:
        #     # sys.stdout = sys.__stdout__
        #     print(seeed, channame, param, peaks[0])
        #     continue

        # Features = fts.modelfeatures(Modeltemp, stim_start=1, stim_end=1.5, apa=False)
        # ScoresWT = fts.modelscore(Modeltemp, meanfeaturesWT, stdfeaturesWT, minfeature=minfeaturesWT, maxfeature=maxfeaturesWT, modelfeature=Features, apa=False)
        # ScoresKO = fts.modelscore(Modeltemp, meanfeaturesKO, stdfeaturesKO, minfeature=minfeaturesKO, maxfeature=maxfeaturesKO, modelfeature=Features, apa=False)
        # # print(Scores)
        # # print('\n')
        # if ScoresWT == False:
        #     # sys.stdout = sys.__stdout__
        #     print(seeed, channame, param, ScoresWT)
        #     continue

        # # print(Modeltemp)

        # #####################################################################
        # # Model = {}
        # Modeltemp["ScoresWT"] = ScoresWT
        # Modeltemp["ScoresKO"] = ScoresKO
        # Modeltemp["Features"] = Features
        # Modeltemp["ErrorWT"] = np.sum(np.array(list(Modeltemp["ScoresWT"].values())) ** 2)/len(Modeltemp["ScoresWT"].values())
        # Modeltemp["ErrorKO"] = np.sum(np.array(list(Modeltemp["ScoresKO"].values())) ** 2)/len(Modeltemp["ScoresKO"].values())

        Modeltemp["AP1_width"] = AP1_width

        f.write("Models['Model" + str(channame) +"_" + str(param) + "'] = " + str(Modeltemp) + "\n\n")
        mnum = mnum + 1
        # sys.stdout = sys.__stdout__
        # print(seeed, i, time.time()-ttt)

f.close()
