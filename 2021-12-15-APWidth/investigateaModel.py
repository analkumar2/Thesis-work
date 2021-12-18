import numpy as np
import matplotlib.pyplot as plt
import MOOSEModel_18 as mm
import moose
import featuresv29_nonallen as fts
import scipy
import pandas as pd
import scipy.interpolate as sciip
from pprint import pprint
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

Modeltemp = {
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
                "Gbar": 0.00012480475550619003*1.5*0.6,
                "Erev": 0.06,
                "Kinetics": "../../Compilations/Kinetics/Na_Chan_Custom4",
                "KineticVars": {
                    "m_vhalf_inf": -0.0316,
                    "h_vhalf_inf": -0.066,
                    "h_C": 0.01197575-0.01,
                    "s_vhalf_inf": -0.045,
                },
            },
            "K_DR_Chan": {
                "Gbar": 1.0502259538910637e-7*10,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_DR_Chan_Custom3",
                "KineticVars": {
                    "n_F": 0.00306,
                    # "n_vhalf_inf": 0.013,
                    # "n_A": 1.26E-02,
                },
            },
            "K_A_Chan": {
                "Gbar": 1.008422244061249e-06,
                "Erev": -0.09,
                "Kinetics": "../../Compilations/Kinetics/K_A_Chan_Custom3",
                "KineticVars": {"n_vhalf_inf": 0.025, "n_slope_inf": 0.017},
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

# t, V, Ca = mm.runModel(Modeltemp, 150e-12, refreshKin=True)
# ARin, ACin = fts.calcRinCin(Modeltemp, refreshKin=False)
rdes = mm.generateModel(Modeltemp, 150e-12)

# soma = moose.element("/model/elec/soma")
# sf_area = np.pi * soma.length * soma.diameter
# sf_area = np.pi*Modeltemp["Parameters"]["Morphology"]["sm_len"]*Modeltemp["Parameters"]["Morphology"]["sm_diam"]
W_Em = -0.065
W_Ca = 80e-6

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
        I_throughChan = G * (
            W_Em - Modeltemp["Parameters"]["Channels"][chan.name]["Erev"]
        )
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
            gateinfX = finterpX(W_Em)
            G = G * gateinfX ** chan.Xpower
        if chan.Ypower > 0:
            DIVSY = np.linspace(gateY.min, gateY.max, len(gateY.tableA))
            finterpY = scipy.interpolate.interp1d(DIVSY, gateY.tableA / gateY.tableB)
            gateinfY = finterpY(W_Em)
            G = G * gateinfY ** chan.Ypower
        if chan.Zpower > 0:
            DIVSZ = np.linspace(gateZ.min, gateZ.max, len(gateZ.tableA))
            finterpZ = scipy.interpolate.interp1d(DIVSZ, gateZ.tableA / gateZ.tableB)
            gateinfZ = finterpZ(W_Em)
            G = G * gateinfZ ** chan.Zpower
        I_throughChan = G * (
            W_Em - Modeltemp["Parameters"]["Channels"][chan.name]["Erev"]
        )
        I_throughChan_sum = I_throughChan_sum + I_throughChan
        A_Gtot = A_Gtot + G
    print(chan, G)

Gin_ = 1e-8 - A_Gtot
if Gin_ <= 0:
    raise SomeError("G leak cannot be -ve")

Modeltemp["Parameters"]["Passive"]["Rm"] = 1 / Gin_
Modeltemp["Parameters"]["Passive"]["Em"] = (I_throughChan_sum + Gin_ * W_Em) / Gin_

print(f'{Modeltemp["Parameters"]["Passive"] = }')

ARin, ACin = fts.calcRinCin(Modeltemp, refreshKin=False)
# Rin_calculate = 1 / A_Gtot
# print(f"{Rin_calculate = }")
print(ARin, ACin)

t150, V150, Ca150 = mm.runModel(Modeltemp, 150e-12, refreshKin=False)
plt.plot(t150, V150)
# t, V, Ca = mm.runModel(Modeltemp, 300e-12, refreshKin=False)
# plt.plot(t, V)
# t, V, Ca = mm.runModel(Modeltemp, 0, refreshKin=False)
# plt.plot(t, V)
plt.show()

# # Importing the means and std
# features_ess_WT = pd.read_csv("features_ess_WT_v3.csv", sep="\t", index_col=0)
# meanfeaturesWT = features_ess_WT.mean()
# stdfeaturesWT = features_ess_WT.std()
# minfeaturesWT = features_ess_WT.min()
# maxfeaturesWT = features_ess_WT.max()

# features_ess_KO = pd.read_csv("features_ess_KO_v3.csv", sep="\t", index_col=0)
# meanfeaturesKO = features_ess_KO.mean()
# stdfeaturesKO = features_ess_KO.std()
# minfeaturesKO = features_ess_KO.min()
# maxfeaturesKO = features_ess_KO.max()

# Features = fts.modelfeatures(Modeltemp, refreshKin=False)
# ScoreKO = fts.modelscore(
#     Modeltemp,
#     meanfeaturesKO,
#     stdfeaturesKO,
#     minfeature=minfeaturesKO,
#     maxfeature=maxfeaturesKO,
#     modelfeature=Features,
#     apa=False,
#     refreshKin=False,
# )

# pprint(Features)
# pprint(ScoreKO)

tt = t150
vv = V150
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
print(f'{AP1_width = }')
