#exec(open('MOOSEModel.py').read())
#Has been modified heavily to be used in Parametersearch_again folder. Do not use as an api

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import xmltodict

def parameterfile_parser(Parameterdict):
    """
    Parses Parameterfile and returns rdesigneur function parameters

    Arguements:
    Parameterfile -- A valid parameterfile address, string
    """
    depth = 0.1
    F = 96485.3329
    Model = Parameterdict['parameters']

    Parameters = {}
    sm_area = np.pi*float(Model['Morphology']['sm_diam'])*float(Model['Morphology']['sm_len'])
    cellProto = [['somaProto', 'soma', float(Model['Morphology']['sm_diam']), float(Model['Morphology']['sm_len'])]]

    chanProto = []
    chanDistrib = []
    chd = Model['Channels']
    for channel in chd.keys():
        chanProto.append([ chd[channel]['kinetics'][:-2]+channel+'()' , channel ])
        chanDistrib.append([ channel, 'soma', 'Gbar', str(chd[channel]['Gbar']/sm_area) ])

    chanProto.append([ Model['Ca_Conc']['kinetics'][:-2]+'Ca_Conc()' , 'Ca_conc' ])
    chanDistrib.append([ 'Ca_conc', 'soma', 'CaBasal', str(Model['Ca_Conc']['Ca_base']), 'tau', str(Model['Ca_Conc']['Ca_tau']) ])

    passiveDistrib = [['soma', 'RM', Model['Model']['segment']['Passive']['@RM'], 'CM', Model['Model']['segment']['Passive']['@CM'], 'initVm', str(-0.065), 'Em', Model['Model']['segment']['Passive']['@Em']]]

    Parameters['cellProto'] = cellProto
    Parameters['chanProto'] = chanProto
    Parameters['chanDistrib'] = chanDistrib
    Parameters['passiveDistrib'] = passiveDistrib
    Parameters['Ca_B'] = float(Model['Model']['segment']['Ca_Conc']['@Ca_B'])

    return Parameters


def generateModel(Parameterfile, CurrInjection):
    """
    Returns in-silico model current clamp

    Arguements:
    Parameterfile -- A valid parameterfile address, string
    CurrInjection -- Current clamp level, float
    """
    Parameters = parameterfile_parser(Parameterfile)
    elecPlotDt = 0.0001
    elecDt = 0.00005
    depth = 0.1
    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        moose.delete('/model')
    except:
        pass

    rdes = rd.rdesigneur(
        elecPlotDt = elecPlotDt,
        elecDt = elecDt,
        cellProto = Parameters['cellProto'],
        chanProto = Parameters['chanProto'],
        passiveDistrib = Parameters['passiveDistrib'],
        chanDistrib = Parameters['chanDistrib'],
        stimList = [['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {CurrInjection} : 0']],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
        ],
    )

    rdes.buildModel()

    #Setup clock table to record time
    clk = moose.element('/clock')
    plott = moose.Table('/model/graphs/plott')
    moose.connect(plott, 'requestOut', clk, 'getCurrentTime')

    #Setting Ca_conc B value
    try:
        moose.element('/model/elec/soma/Ca_conc').B = Parameters['Ca_B']
        # moose.element('/model/elec/soma/Ca_conc').B *= 2
        # moose.element('/model/elec/soma/Ca_conc').B = 0
    except:
        pass

    print('MOOSE Model generated')
    return rdes

def runModel(Parameterfile, CurrInjection):
    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    generateModel(Parameterfile, CurrInjection)
    moose.reinit()
    moose.start(preStimTime+injectTime+postStimTime)

    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector

    return [tvec, Vmvec]

def plotModel(Parameterfile, CurrInjection):
    """
    Returns in-silico model current clamp

    Arguements:
    Parameterfile -- A valid parameterfile address, string
    CurrInjection -- Current clamp level, float
    """
    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    rdes = generateModel(Parameterfile, CurrInjection)
    moose.reinit()
    moose.start(preStimTime+injectTime+postStimTime)

    rdes.display()
    return rdes


if __name__ == '__main__':
    rdes = plotModel('Modelparameters/dummyModel.xml', 150e-12)
    tvec, Vmvec = runModel('Modelparameters/dummyModel.xml', 150e-12)
