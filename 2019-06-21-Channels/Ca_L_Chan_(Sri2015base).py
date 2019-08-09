# Ca_L channel taken from mod files of Srikanth2015 (check once again)
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

def Ca_L_Chan(name):
    Ca_L = moose.HHChannel( '/library/' + name )
    Ca_L.Ek = ECa
    Ca_L.Gbar = 300.0*SOMA_A
    Ca_L.Gk = 0.0
    Ca_L.Xpower = 1.0
    Ca_L.Ypower = 0.0
    Ca_L.Zpower = 1.0
    Ca_L.useConcentration = 1
    Ca_L.instant = 4

    ki=.001
    cai = 50e-6
    cao = 2
    vhalfa = -27.01
    mtfactor = 1

    alpm = 0.055*(vhalfa - v*1e3)/(np.exp((vhalfa-v*1e3)/3.8) - 1)
    betm = 0.94*np.exp((-63.01-v*1e3)/17)
    a = alpm
    taum = mtfactor*1/(5*(a+betm))
    minf = a/(a+betm)

    xgate = moose.element( Ca_L.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum*1e3
    xgate.tableB = 1.0/taum*1e3

    zgate = moose.element( Ca_L.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( Ca_L.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( Ca_L.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return Ca_L
