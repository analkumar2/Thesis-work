# exec(open('changingmInf.py').read())
# 3.2.0-6addb0964  --Vinu's version
# 3.2.0.dev20191210 --My version

import moose
import numpy as np
import pylab
import rdesigneur as rd


# Deleting any previous run of the model
try:
    [moose.delete(x) for x in ['/model', '/library']]
    # moose.delete('/model')
except:
    pass
######################################

rdes = rd.rdesigneur(
    chanProto = [['Na_Chan_(Migliore2018).Na_Chan()', 'Na'], ['K_A_Chan_(Migliore2018)_ghk.K_A_Chan()', 'K']],
    # chanProto = [['make_HH_Na()', 'Na'], ['make_HH_K()', 'K']],
    chanDistrib = [
        ['K', 'soma', 'Gbar', '2000' ],
        ['Na', 'soma', 'Gbar', '100' ],],
    stimList = [['soma', '1', '.', 'inject', '(t>0.1 && t<0.2) * 1e-8' ]],
    plotList = [['soma', '1', '.', 'Vm', 'Membrane potential']]
)

####Initial run
print('Initial run')
elecid_ori = rdes.elecid.path
rdes.buildModel()
moose.reinit()
moose.start( 0.3 )
rdes.display()

####Changing Na mINf kinetics########
Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
v = np.linspace(Vmin,Vmax, Vdivs)

def Na_Chan_minf(v,vshiftm, slopem):
    #vshiftm and slopem should be dimensionless but think of them to be in mV terms
    mAlpha = (0.182 * (v*1e3- (-38+vshiftm)))/(1-(np.exp(-(v*1e3- (-38+vshiftm))/slopem)))
    mBeta  = (0.124 * (-v*1e3 + (-38+vshiftm)))/(1-(np.exp(-(-v*1e3 + (-38+vshiftm))/slopem)))
    mInf = mAlpha/(mAlpha + mBeta)
    return mInf

tmA = moose.element('/library/Na/gateX').tableA
tmB = moose.element('/library/Na/gateX').tableB
moose.element('/library/Na/gateX').tableA = Na_Chan_minf(v,-1,6)*tmB
############################################

print('Second run with Na minf changed')
moose.delete('/model')
rdes.elecid = moose.element(elecid_ori)
rdes.buildModel()
moose.reinit()
moose.start( 0.3 )
rdes.display()
