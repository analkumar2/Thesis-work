# Ca_N channel taken from mod files of Srikanth2015 (check once again)
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

def Ca_N_Chan(name):
    Ca_N = moose.HHChannel( '/library/' + name )
    Ca_N.Ek = ECa
    Ca_N.Gbar = 300.0*SOMA_A
    Ca_N.Gk = 0.0
    Ca_N.Xpower = 2.0
    Ca_N.Ypower = 1.0
    Ca_N.Zpower = 1.0
    Ca_N.useConcentration = 1
    Ca_N.instant = 4

    ki=.001
    cai = 50e-6
    cao = 10
    mtfactor = 1
    htfactor = 1
    vhalfm = 19.88
    vhalfh = 39.0

    alph = 1.6e-4*np.exp(-v*1e3/48.4)
    beth = 1/(np.exp((-v*1e3+vhalfh)/10.)+1.)
    alpm = 0.1967*(-1.0*v*1e3+vhalfm)/(np.exp((-1.0*v*1e3+vhalfm)/10.0)-1.0)
    betm = 0.046*np.exp(-v*1e3/20.73)
    a = alpm
    taum = mtfactor*1/(a + betm)
    minf = a*mtfactor*1/(a + betm)
    a = alph
    tauh = htfactor*1/(a + beth)
    hinf = a*htfactor*1/(a + beth)

    xgate = moose.element( Ca_N.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum*1e3
    xgate.tableB = 1.0/taum*1e3

    ygate = moose.element( Ca_N.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh*1e3
    ygate.tableB = 1.0/tauh*1e3

    zgate = moose.element( Ca_N.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( Ca_N.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( Ca_N.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return Ca_N
