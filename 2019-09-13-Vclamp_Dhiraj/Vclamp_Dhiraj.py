#exec(open('Vclamp_Dhiraj.py').read())
import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

sm_diam = 20e-6
sm_len = 20e-6
dend_diam = 4e-6
dend_len = 500e-6
numDend = 10
sm_area = np.pi*sm_len*sm_diam
elecPlotDt = 50e-6
elecDt = 10e-6
RM = 1
CM = 0.01
RA = 1.5
Vrest = -0.065

try:
    moose.delete('/model')
    moose.delete('/library')
except:
    pass
rdes = rd.rdesigneur(
    elecPlotDt = elecPlotDt,
    elecDt = elecDt,
    # cellProto = [['ballAndStick', 'soma', sm_diam, sm_len, dend_diam, dend_len, numDend]],
    cellProto = [['somaProto', 'soma', sm_diam, sm_len]],
    passiveDistrib = [
        ['#', 'RM', str(RM), 'CM', str(CM), 'initVm', str(Vrest), 'Em', str(Vrest)],
    ],
    stimList = [['soma', '1', '.', 'vclamp', '-0.065 + (t>0.1 && t<0.2) * 0.02' ]],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
    ]
)
rdes.buildModel()
try:
    moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecDt
    moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecDt
    moose.element( '/model/elec/soma/vclamp' ).ti = elecDt*2
    moose.element( '/model/elec/soma/vclamp' ).td = 0
except:
    pass
moose.reinit()
moose.start( 0.3 )
# V2 = moose.element('/model/graphs/plot1').vector
rdes.display()
