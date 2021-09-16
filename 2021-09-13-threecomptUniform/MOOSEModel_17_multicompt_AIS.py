# exec(open('MOOSEModel_15.py').read())
# Has been modified heavily to be used in Parametersearch_again folder. Do not use as an api
# v3 takes care of Ca_B while giving an option to modify kinetics
# v5 added option to free SK channel kinetics. HHChannel2D cannot be parameterized right now.
# v6 Added a vclamp option
# v7 Instead of manually using ChanGate, now kinetic variabless from supported channel rdesigneur files can be directly changed. No backward compatibility
# v8 Fixed the vclamp where the output ws V instead of I
# v9 No longer prints the annoying rdesigneur statements
# v10 deleting a moose element which does not exists now gives a seg fault (current BhallaLab master branch). Modifies the scrit to avoid getting seg faults
# v11 kineticVars issue resolved
# v12 added refreshkinetics or not option. If using the kinetics as used in the kinetics file directly, refreshkin=False. Not properly implemented as of now. Use v11 for now
# v13 Has a synaptic input function
# v14 fixed an issue where ca dynamics B term was not properly changed.
# v14_custom2 for CV_EPSPfreq_amp
# v15 There has been several v14s. v15 compiles all these versions finally
# v16 Now also works if Ca_conc mechanisms not present in the model
# v17 Now you can specify the file and modelnumber to plot directly
# v17_multicompt_uniform For model with ion channels only in AIS. In the parameterdict, specify morphology swc file. Pass gbar, cm, rm, ra instead of their absolute versions. vclamp does not work right now.

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import xmltodict
import sys
import os
import io
import importlib
import numpy.random as nr
import argparse
from pprint import pprint

sys.path.insert(1, "../../Compilations/Kinetics")


elecid_ori = None
elecPlotDt = 0.0001
elecDt = 0.00001


def Parameterdict_parser(Parameterdict):
    """
    Parses Parameterdict and returns rdesigneur function Parameters

    Arguements:
    Parameterdict -- A valid Parameterdict
    """
    depth = 0.1
    F = 96485.3329
    Model = Parameterdict["Parameters"]

    Parameters = {}
    # sm_area = (
    #     np.pi
    #     * float(Model["Morphology"]["sm_diam"])
    #     * float(Model["Morphology"]["sm_len"])
    # )
    # cellProto = [
    #     [
    #         "somaProto",
    #         "soma",
    #         float(Model["Morphology"]["sm_diam"]),
    #         float(Model["Morphology"]["sm_len"]),
    #     ]
    # ]
    cellProto = [Model["Morphology"], 'elec'],

    chanProto = []
    chanDistrib = []
    chd = Model["Channels"]
    for channel in chd.keys():
        chanProto.append([chd[channel]["Kinetics"] + "." + channel + "()", channel])
        chanDistrib.append(
            [channel, "axon_1_2", "Gbar", str(chd[channel]["gbar"])]
        )
        Parameters[channel + "_Erev"] = chd[channel]["Erev"]

    if "Ca_Conc" in Model.keys():
        chanProto.append([Model["Ca_Conc"]["Kinetics"] + ".Ca_Conc()", "Ca_conc"])
        chanDistrib.append(
            [
                "Ca_conc",
                "axon_1_2",
                "CaBasal",
                str(Model["Ca_Conc"]["Ca_base"]),
                "tau",
                str(Model["Ca_Conc"]["Ca_tau"]),
            ]
        )

    passiveDistrib = [
        [
            "#",
            "RM",
            str(Model["Passive"]["RM"]),
            "CM",
            str(Model["Passive"]["CM"]),
            "initVm",
            str(-0.065),
            "Em",
            str(Model["Passive"]["Em"]),
            "RA",
            str(Model["Passive"]["RA"]),
        ]
    ]

    Parameters["cellProto"] = cellProto
    Parameters["chanProto"] = chanProto
    Parameters["chanDistrib"] = chanDistrib
    Parameters["passiveDistrib"] = passiveDistrib
    if "Ca_Conc" in Model.keys():
        Parameters["Ca_B"] = float(Model["Ca_Conc"]["Ca_B"])

    return Parameters


def generateModel(
    Parameterdict, CurrInjection=150e-12, vClamp=None, refreshKin=True, syn=False, synwg=0.05, synfq=5
):
    """
    Returns in-silico model current clamp. Except Ca_B everything is set up

    Arguements:
    Parameterdict -- A valid Parameterdict
    CurrInjection -- Current clamp level, float
    """
    if syn:
        moose.seed()
    global elecid_ori
    Parameters = Parameterdict_parser(Parameterdict)
    depth = 0.1
    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    mooselelist = moose.le()
    if "/model" in mooselelist:
        moose.delete("/model")
    if "/Graphs" in mooselelist:
        moose.delete("/Graphs")

    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        # moose.delete("/model")
        # moose.delete("/Graphs")
        rdes.elecid = moose.element(elecid_ori)
    except:
        pass

    if syn:
        synGbar = 1
    else:
        synGbar = 1e-8

    print(Parameters["cellProto"])
    if vClamp:
        rdes = rd.rdesigneur(
            elecPlotDt=elecPlotDt,
            elecDt=elecDt,
            cellProto=Parameters["cellProto"],
            chanProto=Parameters["chanProto"],
            passiveDistrib=Parameters["passiveDistrib"],
            chanDistrib=Parameters["chanDistrib"],
            stimList=[["soma", "1", ".", "vclamp", vClamp]],
            plotList=[
                ["soma", "1", "vclamp", "current", "Soma holding current"],
                ["soma", "1", "Ca_conc", "Ca", "Soma Calcium concentration"],
            ],
        )
    else:
        rdes = rd.rdesigneur(
            elecPlotDt=elecPlotDt,
            elecDt=elecDt,
            cellProto=Parameters["cellProto"],
            chanProto=Parameters["chanProto"] + [["make_glu()", "glu"]],
            passiveDistrib=Parameters["passiveDistrib"],
            chanDistrib=Parameters["chanDistrib"]
            + [["glu", "soma", "Gbar", f"{synGbar}"]],
            stimList=[
                [
                    "soma",
                    "1",
                    ".",
                    "inject",
                    f"(t>={preStimTime} && t<={preStimTime+injectTime}) ? {CurrInjection} : 0",
                ],
                ["soma", f"{synwg}", "glu", "randsyn", f"{synfq}"],
            ],
            plotList=[
                # ["axon_1_2", "1", ".", "Vm", "AIS Membrane potential MOOSE"],
                ["soma", "1", ".", "Vm", "soma Membrane potential MOOSE"],
                ["soma", "1", "Ca_conc", "Ca", "Soma Calcium concentration"],
                # ["soma", "1", "K_A_Chan", "Gk", "Soma K_A_Chan Gbar"],
            ],
        )

    if refreshKin:  #True if in the new run, the kinetics was changed.
        for chan in Parameterdict["Parameters"]["Channels"].keys():
            imm = Parameterdict["Parameters"]["Channels"][chan]["Kinetics"].split("/")[
                -1
            ]
            exec(f"import {imm}")
            exec(f"importlib.reload({imm})")

            if "KineticVars" in Parameterdict["Parameters"]["Channels"][chan].keys():
                for var in Parameterdict["Parameters"]["Channels"][chan][
                    "KineticVars"
                ].keys():
                    valuee = Parameterdict["Parameters"]["Channels"][chan][
                        "KineticVars"
                    ][var]
                    exec(f"{imm}.{var} = {valuee}")
            exec(f"{imm}.{chan}('{chan}')")

    for chan in moose.wildcardFind("/library/#[CLASS==HHChannel]"):
        moose.element(f"/library/{chan.name}").Ek = Parameters[chan.name + "_Erev"]

    for chan in moose.wildcardFind("/library/#[CLASS==HHChannel2D]"):
        moose.element(f"/library/{chan.name}").Ek = Parameters[chan.name + "_Erev"]

    # Setup clock table to record time
    clk = moose.element("/clock")
    moose.Neutral("Graphs")
    plott = moose.Table("/Graphs/plott")
    moose.connect(plott, "requestOut", clk, "getCurrentTime")

    print("MOOSE Model generated")
    elecid_ori = rdes.elecid.path
    return rdes


def runModel(
    Parameterdict,
    CurrInjection=150e-12,
    vClamp=None,
    refreshKin=True,
    Truntime=None,
    syn=False,
    synwg=0.05, synfq=5
):
    """
    CurrInjection: in A. Put None if using vClamp
    """
    if syn:
        moose.seed()
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")

    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    rdes = generateModel(
        Parameterdict,
        CurrInjection=CurrInjection,
        vClamp=vClamp,
        refreshKin=refreshKin,
        syn=syn, synwg=synwg, synfq=synfq
    )
    rdes.buildModel()

    somaelements = moose.le("/model/elec/soma")
    if "/model[0]/elec[0]/soma[0]/vclamp" in somaelements:
        moose.element("/model/elec/soma/vclamp").gain = (
            moose.element("/model/elec/soma").Cm / elecDt
        )
        moose.element("/model/elec/soma/vclamp").tau = 5 * elecDt
        moose.element("/model/elec/soma/vclamp").ti = elecDt * 2
        moose.element("/model/elec/soma/vclamp").td = 0

    # Setting Ca_conc B value
    AISelements = moose.le("/model/elec/axon_1_2")
    # Parameters = Parameterdict_parser(Parameterdict)
    # if "/model[0]/elec[0]/axon_1_2[0]/Ca_conc" in AISelements:
    #     moose.element("/model/elec/axon_1_2/Ca_conc").B = Parameters["Ca_B"]
    #     # moose.element('/model/elec/soma/Ca_conc').B *= 2
    #     # moose.element('/model/elec/soma/Ca_conc').B = 0

    Parameters = Parameterdict_parser(Parameterdict)
    for Ca_concelement in moose.wildcardFind("/library/#[CLASS==ZombieCaConc]"):
        Ca_concelement.B = Parameters["Ca_B"]

    moose.reinit()
    if Truntime is None:
        moose.start(preStimTime + injectTime + postStimTime)
    else:
        moose.start(Truntime)

    # rdes.display()
    Vmvec = moose.element("/model/graphs/plot0").vector
    tvec = moose.element("/Graphs/plott").vector
    if "/model[0]/elec[0]/soma[0]/Ca_conc" in AISelements:
        Cavec = moose.element("/model/graphs/plot1").vector
    else:
        Cavec = None

    sys.stdout = old_stdout
    return [tvec, Vmvec, Cavec]


def plotModel(
    Parameterdict,
    CurrInjection=150e-12,
    vClamp=None,
    refreshKin=True,
    Truntime=None,
    syn=False,
    synwg=0.05, synfq=5
):
    """
    Returns in-silico model current clamp

    Arguements:
    Parameterdict -- A valid Parameterdict address, string
    CurrInjection -- Current clamp level, float
    """
    moose.seed()
    # old_stdout = sys.stdout
    # sys.stdout = open(os.devnull, "w")

    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    rdes = generateModel(
        Parameterdict,
        CurrInjection=CurrInjection,
        vClamp=vClamp,
        refreshKin=refreshKin,
        syn=syn, synwg=synwg, synfq=synfq
    )
    rdes.buildModel()

    somaelements = moose.le("/model/elec/soma")
    if "/model[0]/elec[0]/soma[0]/vclamp" in somaelements:
        moose.element("/model/elec/soma/vclamp").gain = (
            moose.element("/model/elec/soma").Cm / elecDt
        )
        moose.element("/model/elec/soma/vclamp").tau = 5 * elecDt
        moose.element("/model/elec/soma/vclamp").ti = elecDt * 2
        moose.element("/model/elec/soma/vclamp").td = 0

    # # Setting Ca_conc B value
    AISelements = moose.le("/model/elec/axon_1_2")
    # Parameters = Parameterdict_parser(Parameterdict)
    # if "/model[0]/elec[0]/axon_1_2[0]/Ca_conc" in AISelements:
    #     moose.element("/model/elec/axon_1_2/Ca_conc").B = Parameters["Ca_B"]
    #     # moose.element('/model/elec/soma/Ca_conc').B *= 2
    #     # moose.element('/model/elec/soma/Ca_conc').B = 0

    Parameters = Parameterdict_parser(Parameterdict)
    for Ca_concelement in moose.wildcardFind("/library/#[CLASS==ZombieCaConc]"):
        Ca_concelement.B = Parameters["Ca_B"]

    moose.reinit()
    if Truntime is None:
        moose.start(preStimTime + injectTime + postStimTime)
    else:
        moose.start(Truntime)

    # sys.stdout = old_stdout

    rdes.display()
    return rdes


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str)
    parser.add_argument('-M', type=str, default='Model1', nargs='?')
    parser.add_argument('-I', type=float, default=150e-12, nargs='?')
    args = parser.parse_args()
    exec(open(args.file).read())
    plotModel(Models[args.M], args.I)
    import featuresv26_nonallen_AIS as fts
    from pprint import pprint
    A = fts.modelfeatures(Models[args.M], stim_start=1, stim_end=1.5, apa=False)
    pprint(A)