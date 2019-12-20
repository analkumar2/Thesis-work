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

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
    # alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

Parameters = MOOSEModel_3.Parameterdict_parser(Models['Model38'])
Na_Chan_mInf, Na_Chan_mTau = ChanGate(v,-0.038,0.0072, -0.0365,0.020,0.0161,0.0547,0.0311,0.00064)
Na_Chan_hInf, Na_Chan_hTau = ChanGate(v,-0.050,-0.004, -0.0456,0.00433,0.01198,0.0262,0.00854,0.039)
K_DR_Chan_mInf, K_DR_Chan_mTau = ChanGate(v,0.013,0.0088, 0.0125,0.0173,0,0,0.0341,0.1022)
def wwrapper():
    moose.element('/library/Na_Chan/gateX').tableA = Na_Chan_mInf/Na_Chan_mTau
    moose.element('/library/Na_Chan/gateX').tableB = 1/Na_Chan_mTau

    moose.element('/library/Na_Chan/gateY').tableA = Na_Chan_hInf/Na_Chan_hTau
    moose.element('/library/Na_Chan/gateY').tableB = 1/Na_Chan_hTau

    moose.element('/library/K_DR_Chan/gateX').tableA = K_DR_Chan_mInf/K_DR_Chan_mTau
    moose.element('/library/K_DR_Chan/gateX').tableB = 1/K_DR_Chan_mTau
    ############################################

    print('Second run with Na minf changed')
    moose.delete('/model')
    rdes.elecid = moose.element(elecid_ori)
    rdes.buildModel()

    #Setting Ca_conc B value
    try:
        moose.element('/model/elec/soma/Ca_conc').B = Parameters['Ca_B']
        # moose.element('/library/Ca_conc').B *= 2
        # moose.element('/library/Ca_conc').B = 0
    except:
        pass

    moose.reinit()
    moose.start( 2.5 )
    rdes.display()
