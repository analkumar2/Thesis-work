# Ca_T channel taken from mod files of Srikanth2015 (check once again)
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

def Ca_T_Chan(name):
    Ca_T = moose.HHChannel( '/library/' + name )
    Ca_T.Ek = ECa
    Ca_T.Gbar = 300.0*SOMA_A
    Ca_T.Gk = 0.0
    Ca_T.Xpower = 2.0
    Ca_T.Ypower = 1.0
    Ca_T.Zpower = 0.0

    gCa_Tbar=.003
    cao = 2
    q10 = 5
    mmin=0.2
    hmin=10
    a0h =0.015
    zetah = 3.5
    vhalfh = -75
    gmh=0.6
    a0m =0.04
    zetam = 2
    vhalfm = -28
    gmm=0.61
    vhm=-60
    vhh=-85
    mtfactor=1
    htfactor=1

    alph = np.exp(0.0378*zetah*(v*1e3-vhalfh))
    beth = np.exp(0.0378*zetah*gmh*(v*1e3-vhalfh))
    alpmt = np.exp(0.0378*zetam*(v*1e3-vhalfm))
    betmt = np.exp(0.0378*zetam*gmm*(v*1e3-vhalfm))
    qt=q10**((celsius-25)/10)

    minf = (1/(1 + np.exp(-(v*1e3-vhm)/10)))
    mtau = mtfactor*betmt/(qt*a0m*(1+alpmt))
    mtau[mtau<mmin] = mmin

    hinf = (1/(1 + np.exp((v*1e3-vhh)/10)))
    htau =  htfactor*beth/(a0h*(1+alph))
    htau[htau<hmin] = hmin

    xgate = moose.element( Ca_T.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau*1e3
    xgate.tableB = 1.0/mtau*1e3

    ygate = moose.element( Ca_T.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/htau*1e3
    ygate.tableB = 1.0/htau*1e3

    addmsg1 = moose.Mstring( Ca_T.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return Ca_T
