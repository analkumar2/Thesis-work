#exec(open('MOOSEModel_2.py').read())
#Has been modified heavily to be used in Parametersearch_again folder. Do not use as an api

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import xmltodict

def Parameterdict_parser(Parameterdict):
    """
    Parses Parameterdict and returns rdesigneur function parameters

    Arguements:
    Parameterdict -- A valid Parameterdict
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
        chanProto.append([ chd[channel]['Kinetics']+'.'+channel+'()' , channel ])
        chanDistrib.append([ channel, 'soma', 'Gbar', str(chd[channel]['Gbar']/sm_area) ])
        Parameters[channel+'_Erev'] = chd[channel]['Erev']

    chanProto.append([ Model['Ca_Conc']['Kinetics']+'.Ca_Conc()' , 'Ca_conc' ])
    chanDistrib.append([ 'Ca_conc', 'soma', 'CaBasal', str(Model['Ca_Conc']['Ca_base']), 'tau', str(Model['Ca_Conc']['Ca_tau']) ])

    passiveDistrib = [['soma', 'Rm', str(Model['Passive']['Rm']), 'Cm', str(Model['Passive']['Cm']), 'initVm', str(-0.065), 'Em', str(Model['Passive']['Em'])]]

    Parameters['cellProto'] = cellProto
    Parameters['chanProto'] = chanProto
    Parameters['chanDistrib'] = chanDistrib
    Parameters['passiveDistrib'] = passiveDistrib
    Parameters['Ca_B'] = float(Model['Ca_Conc']['Ca_B'])


    return Parameters


def generateModel(Parameterdict, CurrInjection):
    """
    Returns in-silico model current clamp

    Arguements:
    Parameterdict -- A valid Parameterdict
    CurrInjection -- Current clamp level, float
    """
    Parameters = Parameterdict_parser(Parameterdict)
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
    for chan in moose.wildcardFind('/model/elec/soma/#[CLASS==HHChannel2D]') + moose.wildcardFind('/model/elec/soma/#[CLASS==ZombieHHChannel]'):
        moose.element(f"/model/elec/soma/{chan.name}").Ek = Parameters[chan.name+'_Erev']
        moose.element(f"/library/{chan.name}").Ek = Parameters[chan.name+'_Erev']

    #Setup clock table to record time
    clk = moose.element('/clock')
    plott = moose.Table('/model/graphs/plott')
    moose.connect(plott, 'requestOut', clk, 'getCurrentTime')

    #Setting Ca_conc B value
    try:
        moose.element('/model/elec/soma/Ca_conc').B = Parameters['Ca_B']
        # moose.element('/model/elec/soma/Ca_conc').B *= 2
        # moose.element('/model/elec/soma/Ca_conc').B = 0
        moose.element('/library/Ca_conc').B = Parameters['Ca_B']
    except:
        pass

    print('MOOSE Model generated')
    return rdes

def runModel(Parameterdict, CurrInjection):
    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    generateModel(Parameterdict, CurrInjection)
    moose.reinit()
    moose.start(preStimTime+injectTime+postStimTime)

    Vmvec=moose.element('/model/graphs/plot0').vector
    tvec=moose.element('/model/graphs/plott').vector

    return [tvec, Vmvec]

def plotModel(Parameterdict, CurrInjection):
    """
    Returns in-silico model current clamp

    Arguements:
    Parameterdict -- A valid Parameterdict address, string
    CurrInjection -- Current clamp level, float
    """
    preStimTime = 1
    injectTime = 0.5
    postStimTime = 1

    rdes = generateModel(Parameterdict, CurrInjection)
    moose.reinit()
    moose.start(preStimTime+injectTime+postStimTime)

    rdes.display()
    return rdes


if __name__ == '__main__':
    rdes = plotModel(Models['Model1'], 150e-12)
    tvec, Vmvec = runModel(Models['Model1'], 150e-12)
