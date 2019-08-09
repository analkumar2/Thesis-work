#exec(open('/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work/Bianchi2002_WT_CA1/test_vclamp.py').read())
import moose
import rdesigneur as rd

try:
    [moose.delete(x) for x in ['/model', '/library']]
except:
    pass
    
    
rdes = rd.rdesigneur(
    cellProto = [['somaProto', 'soma', 20e-6, 200e-6]],
    stimList = [
        ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<8) * 0.02' ],
        # ['soma', '1', '.', 'inject', '(t>=3.0 && t<=8.003) ? 150e-12 : 0'],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma membrane potential'],
        ['soma', '1', 'vclamp', 'current', 'Soma holding current'],
    ]
)

rdes.buildModel()
# moose.element( '/model/elec/soma/vclamp' ).gain *= 0.02

# addmsg1 = moose.Mstring( KsAHP.path + '/addmsg3' )
# addmsg1.value = '../Ca_conc    concOut    . concen'

# soma = moose.element('/model/elec/soma')
# clamp = moose.element('/model/elec/soma/vclamp')
# moose.connect(clamp, 'currentOut', soma, 'injectMsg')
# moose.connect(soma, 'VmOut', clamp, 'sensedIn')

moose.reinit()
moose.start( 10 )
rdes.display()