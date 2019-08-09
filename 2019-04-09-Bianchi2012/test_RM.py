#exec(open('/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work/Bianchi2002_WT_CA1/test_RM.py').read())
import moose
import rdesigneur as rd

try:
    [moose.delete(x) for x in ['/model', '/library']]
except:
    pass
    
sm_diam = 20e-6
sm_len = 200e-6
sm_sarea = np.pi*sm_diam*sm_len
    
rdes = rd.rdesigneur(
    cellProto = [['somaProto', 'soma', 20e-6, 200e-6]],
    passiveDistrib = [
        # ['soma', 'RM', '20000', 'Cm', '120e-12', 'initVm', str(Em), 'Em', str(Em)],
        ['soma', 'Rm', str(20000/sm_sarea), 'Cm', '120e-12', 'initVm', str(Em), 'Em', str(Em)]
    ],
    stimList = [
        # ['soma', '1', '.', 'vclamp', '-0.065 + (t>3 && t<8) * 0.02' ],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Soma membrane potential'],
    ]
)

rdes.buildModel()
moose.reinit()
moose.start( 10 )
rdes.display()