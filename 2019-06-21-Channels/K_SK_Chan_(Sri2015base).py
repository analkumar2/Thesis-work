# K_SK channel taken from mod files of Srikanth2015 (check once again)
# Problems: qt not proper

import numpy as np
import pickle
import pandas as pd
import moose

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.055
EK = -0.090
Eh = -0.030
ECa = 0.140
Em = -0.065

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
dV = (Vmax-Vmin)/Vdivs
v = np.arange(Vmin,Vmax+dV, dV)
Camin = 1e-12
Camax = 3e-3
Cadivs = 3000
dCa = (Camax-Camin)/Cadivs
ca = np.arange(Camin,Camax+dCa, dCa)

def K_SK_Chan(name):
    K_SK = moose.HHChannel( '/library/' + name )
    K_SK.Ek = EK
    K_SK.Gbar = 300.0*SOMA_A
    K_SK.Gk = 0.0
    K_SK.Xpower = 0.0
    K_SK.Ypower = 0.0
    K_SK.Zpower = 3.0
    K_SK.useConcentration = 1

    n=4
    cai=50.e-6
    a0=1.3e13
    b0=.5e-2
    tfactor = 1
    cahalf=140.04e-6
    k=0.10857

    alp = a0*ca**n
    a = alp
    tau = tfactor*1/(a + b0)
    inf=1/(1+np.exp((np.log(cahalf/ca))/k))

    zgate = moose.element( K_SK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = inf/tau*1e3
    zgate.tableB = 1.0/tau*1e3

    addmsg3 = moose.Mstring( K_SK.path + '/addmsg3' )
    addmsg3.value = '../Ca_conc    concOut    . concen'
    return K_SK
