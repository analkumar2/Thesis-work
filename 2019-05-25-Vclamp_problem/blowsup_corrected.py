import moose
import rdesigneur as rd
moose.delete('/model')
rdes = rd.rdesigneur(
    cellProto = [['somaProto', 'soma', 20e-6, 200e-6]],
    stimList = [['soma', '1', '.', 'vclamp', '-0.065 + (t>0.1 && t<0.2) * 0.02' ]],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
    ]
)
rdes.buildModel()
moose.element( '/model/elec/soma/vclamp' ).gain *= sm_area/7.85e-7
moose.reinit()
moose.start( 0.3 )
rdes.display()
