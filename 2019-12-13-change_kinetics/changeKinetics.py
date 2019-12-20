# exec(open('changeKinetics.py').read())
# 3.2.0-6addb0964  --Vinu's version
# 3.2.0.dev20191210 --My version

import moose
import numpy as np
import pylab
import rdesigneur as rd
import MOOSEModel_3

exec(open('../2019-12-05-Parametersearch_again/Modelparameters_temp/outputModels_dict_4179468691.py').read())
rdes = MOOSEModel_3.generateModel(Models['Model38'], 150e-12)

####Initial run
print('Initial run')
elecid_ori = rdes.elecid.path
rdes.buildModel()

#Setting Ca_conc B value
Parameters = MOOSEModel_3.Parameterdict_parser(Models['Model38'])
try:
    moose.element('/model/elec/soma/Ca_conc').B = Parameters['Ca_B']
    # moose.element('/library/Ca_conc').B *= 2
    # moose.element('/library/Ca_conc').B = 0
except:
    pass

moose.reinit()
moose.start( 2.5 )
rdes.display()

def wwrapper(Na_vshiftm, Na_slopem):
    ####Changing Na mINf kinetics########
    Vmin = -0.100
    Vmax = 0.100
    Vdivs = 3000
    v = np.linspace(Vmin,Vmax, Vdivs)

    def Na_Chan_m(v,vshiftm, slopem):
        qt = 2.3**((34-21)/10)
        #vshiftm and slopem should be dimensionless but think of them to be in mV terms
        mAlpha = (0.182 * (v*1e3- (-38+vshiftm)))/(1-(np.exp(-(v*1e3- (-38+vshiftm))/slopem)))
        mBeta  = (0.124 * (-v*1e3 + (-38+vshiftm)))/(1-(np.exp(-(-v*1e3 + (-38+vshiftm))/slopem)))
        mInf = mAlpha/(mAlpha + mBeta)
        mTau = (1/(mAlpha + mBeta))/qt
        return [mInf,mTau*1e-3]

    def Na_Chan_h(v,vshifth, slopeh):
        qt = 2.3**((34-21)/10)
        #vshiftm and slopem should be dimensionless but think of them to be in mV terms
        hAlpha = (-0.015 * (v*1e3- (-66+vshifth)))/(1-(np.exp((v*1e3- (-66+vshifth))/slopeh)))
        hBeta  = (-0.015 * (-v*1e3 +(-66+vshifth)))/(1-(np.exp((-v*1e3 +(-66+vshifth))/slopeh)))
        hInf = hAlpha/(hAlpha + hBeta)
        hTau = (1/(hAlpha + hBeta))/qt
        return [hInf,hTau*1e-3]

    tmA = moose.element('/library/Na_Chan/gateX').tableA
    tmB = moose.element('/library/Na_Chan/gateX').tableB
    mInf, mTau = Na_Chan_minf(v,Na_vshiftm, Na_slopem)
    moose.element('/library/Na_Chan/gateX').tableA = mInf/mTau
    moose.element('/library/Na_Chan/gateX').tableB = 1/mTau
    # moose.element('/library/Na_Chan/gateX').tableA = tmA
    ############################################

    print('Second run with Na minf changed')
    moose.delete('/model')
    rdes.elecid = moose.element(elecid_ori)
    rdes.buildModel()

    #Setting Ca_conc B value
    Parameters = MOOSEModel_3.Parameterdict_parser(Models['Model38'])
    try:
        moose.element('/model/elec/soma/Ca_conc').B = Parameters['Ca_B']
        # moose.element('/library/Ca_conc').B *= 2
        # moose.element('/library/Ca_conc').B = 0
    except:
        pass

    moose.reinit()
    moose.start( 2.5 )
    rdes.display()
