#Author: Anal Kumar
#Uses the parameters in Traub et al, 1991 to model a reduced CA3 neuron
#Morphology was already set up by Harshavardan
#Free paprmeters are B( or thickness) and Ca_scale. which are to be varied

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

EREST = -0.060
# EREST = 0
F = 96485.3329

#R is the radius of the compartment. B is the scaling term.
#length is the length of the compartment.
#The function finds the thickness of the Ca shell when the scaling
#term is given.
def findThickness(R, length, B):
    thick = R - np.sqrt(R**2 - 1.0/(2*B*F*np.pi*length))
    return(str(thick))

try:
    [moose.delete(x) for x in ['/model', '/library']]
except:
    pass
    
rdes = rd.rdesigneur(
    cellProto = [
        ['Compartments.swc','elec']
    ],
        
    chanProto = [
        ['Channelprotos.NaChan()', 'Na_Chan'],
        ['Channelprotos.KdrChan()', 'Kdr_Chan'],
        ['Channelprotos.KaChan()', 'Ka_Chan'],
        ['Channelprotos.CaConc()', 'Ca_conc'],
        ['Channelprotos.CaChan()', 'Ca_Chan'],        
        ['Channelprotos.KahpChan()', 'Kahp_Chan'],
        ['Channelprotos.KcChan()', 'Kc_Chan'],
    ],
        
    passiveDistrib = [
        ['soma', 'Rm', '301002256.4', 'Ra', '2223716.4', 'Cm', '9.96E-11', 'initVm', str(EREST), 'Em', str(EREST)],
        ['apical#', 'RM', '1', 'RA', '1', 'CM', '0.03', 'initVm', str(EREST), 'Em', str(EREST)],
        ['dend#', 'RM', '1', 'RA', '1', 'CM', '0.03', 'initVm', str(EREST), 'Em', str(EREST)],
    ],
        
    chanDistrib = [
        # soma
        ['Na_Chan', 'soma', 'Gbar', '300'],
        ['Kdr_Chan', 'soma', 'Gbar', '250'],
        ['Ka_Chan', 'soma', 'Gbar', '50'],
        ['Ca_Chan', 'soma', 'Gbar', '40'],
        ['Ca_conc', 'soma', 'thick', findThickness(16.25e-6, 3.25e-05, 17.402e12)],
        ['Kahp_Chan', 'soma', 'Gbar', '8'],
        ['Kc_Chan', 'soma', 'Gbar', '100'],
        
        # apical dendrtites
        ['Na_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 150 : ((p>=240e-6 && p<=360e-6) ? 200 : 0)'],
        ['Kdr_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 100 : ((p>=120e-6 && p<=240e-6) ? 50 : ((p>=240e-6 && p<=360e-6) ? 200 : 0))'],
        ['Ca_conc', 'apical#', 'thick', '(p<=120e-6) ? %s : %s' %(findThickness(2.89e-6, 0.00012, 26.404e12), findThickness(2.89e-6, 0.00012, 5.941e12))],
        ['Ca_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 80 : ((p>=120e-6 && p<=240e-6) ? 50 : ((p>=240e-6 && p<=360e-6) ? 170 : ((p>=360e-6 && p<=720e-6) ? 70 : ((p>=720e-6 && p<=1080e-6) ? 50 : 0))))'],
        ['Kahp_Chan', 'apical#', 'Gbar', '(p<=1080e-6) ? 8 : 0'],
        ['Kc_Chan', 'apical#', 'Gbar', '(p<=120e-6) ? 200 : ((p>=120e-6 && p<=240e-6) ? 50 : ((p>=240e-6 && p<=360e-6) ? 150 : ((p>=360e-6 && p<=1080e-6) ? 50 : 0)))'],
        
        # basal dendrites
        ['Na_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 150 : ((p>=220e-6 && p<=330e-6) ? 200 : 0)'],
        ['Kdr_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 100 : ((p>=110e-6 && p<=220e-6) ? 50 : ((p>=220e-6 && p<=330e-6) ? 200 : 0))'],
        ['Ca_conc', 'dend#', 'thick', '(p<=110e-6) ? %s : %s' %(findThickness(2.42e-6, 0.00011, 34.530e12), findThickness(2.42e-6, 0.00011, 7.769e12))],
        ['Ca_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 80 : ((p>=110e-6 && p<=220e-6) ? 50 : ((p>=220e-6 && p<=330e-6) ? 120 : ((p>=330e-6 && p<=550e-6) ? 70 : ((p>=550e-6 && p<=770e-6) ? 50 : 0))))'],
        ['Kahp_Chan', 'dend#', 'Gbar', '(p<=770e-6) ? 8 : 0'],
        ['Kc_Chan', 'dend#', 'Gbar', '(p<=110e-6) ? 200 : ((p>=110e-6 && p<=220e-6) ? 50 : ((p>=220e-6 && p<=330e-6) ? 100 : ((p>=330e-6 && p<=770e-6) ? 50 : 0)))'],
    
    ],
        
    stimList = [
        ['soma', '1', '.', 'inject', '(t>=3.000 && t<4.000) ? 150e-12 : 0' ],
    ],
    
    # plotList = [
        #soma
        # ['soma', '1', '.', 'Vm', 'Soma Membrane potential'],        
        # ['soma', '1', '.', 'Im', 'Soma Membrane Current'],
        # ['soma', '1', 'Ca_Chan', 'Ik', 'Soma Calcium current'],
        # ['soma', '1', 'Ca_Chan', 'Gk', 'Soma Calcium conductance'],
        # ['soma', '1', 'Ca_conc', 'Ca', 'Soma Calcium concentration'],
        # ['soma', '1', 'Na_Chan', 'Ik', 'Soma Sodium current'],
        # ['soma', '1', 'Na_Chan', 'Gk', 'Soma Sodium conductance'],
        # ['soma', '1', 'Kdr_Chan', 'Ik', 'Soma Kdr current'],
        # ['soma', '1', 'Kdr_Chan', 'Gk', 'Soma Kdr conductance'],
        # ['soma', '1', 'Ka_Chan', 'Ik', 'Soma Ka current'],
        # ['soma', '1', 'Ka_Chan', 'Gk', 'Soma Ka conductance'],
        # ['soma', '1', 'Kc_Chan', 'Ik', 'Soma Kc current'],
        # ['soma', '1', 'Kc_Chan', 'Gk', 'Soma Kc conductance'],
        # ['soma', '1', 'Kahp_Chan', 'Ik', 'Soma Kahp current'],
        # ['soma', '1', 'Kahp_Chan', 'Gk', 'Soma Kahp conductance'],
        # ['soma', '1', 'Kahp_Chan', 'Z', 'Soma Kahp channel Z gate'],
        
        #apical 3
        # ['apical_1_2', '1', '.', 'Vm', 'Apical 3 Membrane potential'],
        # ['apical_1_2', '1', 'Ca_conc', 'Ca', 'Apical 3 Calcium concentration'],
        
        #apical 6
        # ['apical_1_5', '1', '.', 'Vm', 'Apical 6 Membrane potential'],
    # ]
        
    # moogList = [
        # ['#', '1', '.', 'Vm', 'Soma potential'],
    # ]
)

rdes.buildModel()


moose.reinit()

data = moose.Neutral('/data')
# somaKdrcurr = moose.Table('/data/somaKdrcurr')
# somaKacurr = moose.Table('/data/somaKacurr')
# somaNacurr = moose.Table('/data/somaNacurr')
# somaCacurr = moose.Table('/data/somaCacurr')
somaKahpZgate = moose.Table('/data/somaKahpZgate')
somaVmdata = moose.Table('/data/somaVmdata')

# somaKdr = moose.element('/model/elec/soma/Kdr_Chan')
# somaKa = moose.element('/model/elec/soma/Ka_Chan')
# somaNa = moose.element('/model/elec/soma/Na_Chan')
# somaCa = moose.element('/model/elec/soma/Ca_Chan')
somaKahp = moose.element('/model/elec/soma/Kahp_Chan')
somaVm = moose.element('/model/elec/soma')

# moose.connect(somaKdrcurr, 'requestOut', somaKdr, 'getIk')
# moose.connect(somaKacurr, 'requestOut', somaKa, 'getIk')
# moose.connect(somaNacurr, 'requestOut', somaNa, 'getIk')
# moose.connect(somaCacurr, 'requestOut', somaCa, 'getIk')
moose.connect(somaKahpZgate, 'requestOut', somaKahp, 'getZ')
moose.connect(somaVmdata, 'requestOut', somaVm, 'getVm')

moose.start( 10 )

# plt.plot(somaKdrcurr.vector+somaKacurr.vector+somaNacurr.vector+somaCacurr.vector)
# plt.plot(somaKahpZgate.vector)
rdes.display()

somaVm = somaVmdata.vector
somaVmdata1 = np.gradient(somaVm,0.0001)
plt.figure(10)
plt.plot(np.linspace(0,10,100000), somaVm)
plt.xlabel('time (s)')
plt.ylabel('Membrane potential (V)')
plt.figure(11)
plt.plot(np.linspace(0,10,100000), somaVmdata1)
plt.show()