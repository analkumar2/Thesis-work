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
            "Rm": 609705941.9371204,
            "Em": 0.06022479405782892,
        },
        "Channels": {
            "Na_Chan": {
                "Gbar": 1.87e-4,
                "Erev": 0.06,
                "Kinetics": "../../Compilations/Kinetics/Na_Chan_Custom4",
                "KineticVars": {
                    "m_vhalf_inf": -0.0316,
                    "h_vhalf_inf": -0.066,
                    "s_vhalf_inf": -0.033,
                },
            },
            "K_DR_Chan": {
                "Gbar": 1.0502259538910637e-7,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_DR_Chan_Custom3",
                "KineticVars": {
                    # "n_F": 0.00306,
                    # "n_vhalf_inf": 0.0112,
                    # "n_slope_inf": 0.017,
                    "n_A": -8.78e-3,
                    "n_B": 5.63e-2,
                    "n_C": 0,
                    "n_D": 0,
                    "n_E": 2.65e-2,
                    "n_F": 1.05e-2,
                },
            },
            "K_A_Chan": {
                "Gbar": 1.008422244061249e-06,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_A_Chan_Custom3",
                "KineticVars": {"n_vhalf_inf": 0.025},
            },
            "K_M_Chan": {
                "Gbar": 8.032e-09,
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
sys.stdout = open(os.devnull, "w")
rdes = mm.generateModel(dummyModel, 150e-12)  # Dummy model to load the channel kinetics
rdes.buildModel()
sys.stdout = sys.__stdout__

# Actual for loop
num_of_iterations = 2500
mnum = 1
for i in range(num_of_iterations):
    # print(f'iternum = {i}')
    ttt = time.time()
    sys.stdout = open(os.devnull, "w")
    Modeltemp = deepcopy(dummyModel)
    sm_area = 14.2e-9  # Does not matter as a random value of Gbars is choosen anyway

    ##########Initializing each parameter randomly##############
    Modeltemp["Scores"] = {}
    # Modeltemp["Parameters"] = {}
    Modeltemp["Parameters"]["notes"] = ""
    # Modeltemp["Parameters"]["Passive"] = {}
    # Modeltemp["Parameters"]["Passive"]["Cm"] = 119e-12  # nr.uniform(40e-12, 300e-12)
    # Modeltemp["Parameters"]["Passive"]["RM"] = nr.uniform(1, 3)
    # Modeltemp["Parameters"]["Passive"]["Em"] = nr.uniform(-0.065, 0)
    # Modeltemp["Parameters"]["Channels"] = {}
    # for chan in moose.wildcardFind(
    #     "/model/elec/soma/#[CLASS==HHChannel2D]"
    # ) + moose.wildcardFind("/model/elec/soma/#[CLASS==ZombieHHChannel]"):
    #     Modeltemp["Parameters"]["Channels"][chan.name] = {}
    #     Modeltemp["Parameters"]["Channels"][chan.name]["Gbar"] = moose.element(
    #         f"/model/elec/soma/{chan.name}"
    #     ).Gbar
    #     Modeltemp["Parameters"]["Channels"][chan.name]["Erev"] = moose.element(
    #         f"/model/elec/soma/{chan.name}"
    #     ).Ek
    #     Modeltemp["Parameters"]["Channels"][chan.name]["Kinetics"] = [
    #         i for i in np.ravel(rdes.chanProtoList) if chan.name + "()" in i
    #     ][0][::-1].partition(".")[2][::-1]
    # Modeltemp["Parameters"]["Ca_Conc"] = {}
    # Modeltemp["Parameters"]["Ca_Conc"]["Ca_B"] = 10 ** nr.uniform(
    #     7.26, 11.26
    # )  # nr.uniform(1.8e7, 1.8e11)
    # Modeltemp["Parameters"]["Ca_Conc"]["Ca_tau"] = nr.uniform(0.001, 0.200) #0.150  # nr.uniform(0.001, 0.200)
    # Modeltemp["Parameters"]["Ca_Conc"]["Ca_base"] = 50e-6  # nr.uniform(1e-5, 10e-5)
    # Modeltemp["Parameters"]["Ca_Conc"]["Kinetics"] = [
    #     i for i in np.ravel(rdes.chanProtoList) if "Ca_Conc()" in i
    # ][0][::-1].partition(".")[2][::-1]
    # Modeltemp["Parameters"]["Channels"]["Ca_T_Chan"]["Gbar"] = (
    #     10 ** nr.uniform(-4, 2.3) * sm_area
    # )  # nr.uniform(0.00001, 20)*sm_area
    Modeltemp["Parameters"]["Channels"]["Na_Chan"]["Gbar"] = 10 ** (
        nr.uniform(-5, -3)
    )
    Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["Gbar"] = (
        10 ** nr.uniform(-8, -6)
    )
    Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["Gbar"] = (
        10 ** nr.uniform(-7, -5)
    )
    Modeltemp["Parameters"]["Channels"]["K_M_Chan"]["Gbar"] = (
        10 ** nr.uniform(-9.5, -7)
    )
    Modeltemp["Parameters"]["Channels"]["h_Chan"]["Gbar"] = (
        10 ** nr.uniform(-12, -10)
    )
    # Modeltemp["Parameters"]["Channels"]["K_SK_Chan"]["Gbar"] = nr.uniform(0.0129, 1.29) * sm_area #0.129*sm_area #(
        # 10 ** nr.uniform(-3, 2.5) * sm_area)  # nr.uniform(0.00001, 1000)*sm_area
    EK = -0.090  # nr.uniform(-0.090, -0.080)
    Eh = -0.040  # nr.uniform(-0.050, -0.030)
    ENa = 0.060  # nr.uniform(0.050, 0.100)
    ECa = 0.120  # nr.uniform(0.120, 0.140)
    # Modeltemp["Parameters"]["Channels"]["Ca_T_Chan"]["Erev"] = ECa
    Modeltemp["Parameters"]["Channels"]["Na_Chan"]["Erev"] = ENa
    Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["Erev"] = EK
    Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["Erev"] = EK
    Modeltemp["Parameters"]["Channels"]["K_M_Chan"]["Erev"] = EK
    Modeltemp["Parameters"]["Channels"]["h_Chan"]["Erev"] = Eh
    # Modeltemp["Parameters"]["Channels"]["K_SK_Chan"]["Erev"] = EK

    # shiftKin = nr.uniform(0, 0.035)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"] = {}
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["m_vhalf_inf"] = -31.6e-3+nr.uniform(-0.010, 0.010)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["m_slope_inf"] = 0.008*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["m_A"] = -0.043+nr.uniform(-0.005, 0.005)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["m_F"] = 6.40e-04*10**nr.uniform(-1,1)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["h_vhalf_inf"] = -0.066+nr.uniform(-0.015, 0.015)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["h_slope_inf"] = -0.0053*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["h_A"] = -0.051+nr.uniform(-0.005, 0.005)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["h_F"] = 0.0965*10**nr.uniform(-1,1)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["s_vhalf_inf"] = -33e-3+nr.uniform(-0.010, 0.015)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["s_slope_inf"] = -0.006*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]["s_D"] = 0.5*nr.uniform(0.25,2)

    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"] = {}
    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]["n_vhalf_inf"] = 0.013+nr.uniform(-0.005, 0.005)
    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]["n_slope_inf"] = 0.0088*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]["n_A"] = 1.26E-02+nr.uniform(-0.005, 0.005)
    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]["n_B"] = 1.73E-02*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]["n_E"] = 3.43e-02 * 0.5*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]["n_F"] = 1.02e-1 * 10**nr.uniform(-3,-1)

    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"] = {}
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["n_vhalf_inf"] = 0.0228+nr.uniform(-0.005, 0.005)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["n_slope_inf"] = 0.017*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["n_A"] = -8.78e-3+nr.uniform(-0.025,0.025)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["n_B"] = 5.63e-2*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["n_E"] = 0.04725089723705451 * 0.1*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["n_F"] = 0.02526175002229309 * 0.05*10**nr.uniform(-1,1)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["l_vhalf_inf"] = -0.045+nr.uniform(-0.005,0.005)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["l_slope_inf"] = -0.0157*nr.uniform(0.5,2)
    # Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]["l_m"] = 0.26*10**nr.uniform(-1,1)

    # Modeltemp["Parameters"]["Channels"]["K_SK_Chan"]["KineticVars"] = {}
    # Modeltemp["Parameters"]["Channels"]["K_SK_Chan"]["KineticVars"]["cac"] = 500e-6
    #####################################################################
    W_Erest = nr.uniform(minfeaturesWT['E_rest_0'], maxfeaturesWT['E_rest_0'])
    W_Rin = nr.uniform(minfeaturesWT['Input resistance'], maxfeaturesWT['Input resistance'])
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
        print(chan, G)

    Gl = W_Gtot - A_Gtot
    if Gl <=0:
        sys.stdout = sys.__stdout__
        print(seeed, i, Gl)
        continue

    Modeltemp["Parameters"]["Passive"]["Rm"] = 1/Gl
    Modeltemp["Parameters"]["Passive"]["Em"] = (I_throughChan_sum + Gl*W_Erest)/Gl

    ######################################################################
    tvec, Vmvec, Cavec = mm.runModel(Modeltemp, CurrInjection=150e-12)
    peaks = scs.find_peaks(x=Vmvec, height=None, threshold=None, distance=None, prominence=0.050, width=None, wlen=None, rel_height=0.5, plateau_size=None)
    if len(peaks[0])<4:
        sys.stdout = sys.__stdout__
        print(seeed, i, peaks[0])
        continue

    Features = fts.modelfeatures(Modeltemp, stim_start=1, stim_end=1.5, apa=False)
    ScoresWT = fts.modelscore(Modeltemp, meanfeaturesWT, stdfeaturesWT, minfeature=minfeaturesWT, maxfeature=maxfeaturesWT, modelfeature=Features, apa=False)
    ScoresKO = fts.modelscore(Modeltemp, meanfeaturesKO, stdfeaturesKO, minfeature=minfeaturesKO, maxfeature=maxfeaturesKO, modelfeature=Features, apa=False)
    # print(Scores)
    # print('\n')
    if ScoresWT == False:
        sys.stdout = sys.__stdout__
        print(seeed, i, ScoresWT)
        continue

    # print(Modeltemp)

    #####################################################################
    # Model = {}
    Modeltemp["ScoresWT"] = ScoresWT
    Modeltemp["ScoresKO"] = ScoresKO
    Modeltemp["Features"] = Features
    Modeltemp["ErrorWT"] = np.sum(np.array(list(Modeltemp["ScoresWT"].values())) ** 2)/len(Modeltemp["ScoresWT"].values())
    Modeltemp["ErrorKO"] = np.sum(np.array(list(Modeltemp["ScoresKO"].values())) ** 2)/len(Modeltemp["ScoresKO"].values())
    # Model["Parameters"] = deepcopy(Modeltemp["Parameters"])
    # Model["Parameters"] = {}
    # Model["Parameters"]["notes"] = ""
    # Model["Parameters"]["Morphology"] = {}
    # Model["Parameters"]["Morphology"]["sm_len"] = moose.element(
    #     "/model/elec/soma"
    # ).length
    # Model["Parameters"]["Morphology"]["sm_diam"] = moose.element(
    #     "/model/elec/soma"
    # ).diameter
    # Model["Parameters"]["Passive"] = {}
    # Model["Parameters"]["Passive"]["Cm"] = moose.element("/model/elec/soma").Cm
    # Model["Parameters"]["Passive"]["Rm"] = moose.element("/model/elec/soma").Rm
    # Model["Parameters"]["Passive"]["Em"] = moose.element("/model/elec/soma").Em
    # Model["Parameters"]["Channels"] = {}
    # for chan in moose.wildcardFind(
    #     "/model/elec/soma/#[CLASS==HHChannel2D]"
    # ) + moose.wildcardFind("/model/elec/soma/#[CLASS==ZombieHHChannel]"):
    #     Model["Parameters"]["Channels"][chan.name] = {}
    #     Model["Parameters"]["Channels"][chan.name]["Gbar"] = moose.element(
    #         f"/model/elec/soma/{chan.name}"
    #     ).Gbar
    #     Model["Parameters"]["Channels"][chan.name]["Kinetics"] = [
    #         i for i in np.ravel(rdes.chanProtoList) if chan.name + "()" in i
    #     ][0][::-1].partition(".")[2][::-1]
    #     Model["Parameters"]["Channels"][chan.name]["Erev"] = moose.element(
    #         f"/model/elec/soma/{chan.name}"
    #     ).Ek
    # # Model["Parameters"]["Ca_Conc"] = {}
    # Model["Parameters"]["Ca_Conc"]["Ca_B"] = moose.element("/model/elec/soma/Ca_conc").B
    # Model["Parameters"]["Ca_Conc"]["Ca_tau"] = moose.element(
    #     "/model/elec/soma/Ca_conc"
    # ).tau
    # Model["Parameters"]["Ca_Conc"]["Ca_base"] = moose.element(
    #     "/model/elec/soma/Ca_conc"
    # ).Ca_base
    # Model["Parameters"]["Ca_Conc"]["Kinetics"] = [
    #     i for i in np.ravel(rdes.chanProtoList) if "Ca_Conc()" in i
    # ][0][::-1].partition(".")[2][::-1]
    # Model["Parameters"]["Channels"]["Na_Chan"]["KineticVars"] = Modeltemp["Parameters"]["Channels"]["Na_Chan"]["KineticVars"]
    # Model["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"] = Modeltemp["Parameters"]["Channels"]["K_DR_Chan"]["KineticVars"]
    # Model["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"] = Modeltemp["Parameters"]["Channels"]["K_A_Chan"]["KineticVars"]
    # Model["Parameters"]["Channels"]["K_SK_Chan"]["KineticVars"] = Modeltemp["Parameters"]["Channels"]["K_SK_Chan"]["KineticVars"]

    f.write("Models['Model" + str(mnum) + "'] = " + str(Modeltemp) + "\n\n")
    mnum = mnum + 1
    sys.stdout = sys.__stdout__
    print(seeed,mnum)
    # print(seeed, i, time.time()-ttt)

f.close()
