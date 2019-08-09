
import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

[moose.delete(x) for x in ['/model']]

diameter = sm_diam
elecPlotDt = 0.00005
RA = 1
RM = 1
CM = 0.03
Em = -0.07

rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    cellProto = [
        ['6compartment/Compartments.swc','elec']
    ],
    passiveDistrib = [
        ['soma', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
        ['apical#', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
        ['dend#', 'RM', str(RM), 'RA', str(RA), 'CM', str(CM), 'initVm', str(Em), 'Em', str(Em)],
    ],
    stimList = [
        ['soma', '1', '.', 'vclamp', f'-0.065 + (t>{1} && t<{1.5}) * 0.050' ],
    ],

    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
    ],
)

rdes.buildModel()
sm_len = moose.element('/model/elec/soma').length
sm_diam = moose.element('/model/elec/soma').diameter
sm_area = np.pi*sm_diam*sm_len

try:
    moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass

moose.reinit()
moose.start( 2 )
rdes.display()
