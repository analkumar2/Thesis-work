# Ca_R channel taken from mod files of Srikanth2015 (check once again)
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

def Ca_R_Chan(name):
    Ca_R = moose.HHChannel( '/library/' + name )
    Ca_R.Ek = ECa
    Ca_R.Gbar = 300.0*SOMA_A
    Ca_R.Gk = 0.0
    Ca_R.Xpower = 3.0
    Ca_R.Ypower = 1.0
    Ca_R.Zpower = 0.0

    q10  = 3
    taum_exp = 0.92
    z = 2
    vhalfm = 3
    vhalfh = -39
    mtfactor = 1
    htfactor = 1

    qt=q10**(-(celsius-22)/10)
    minf = 1/(1+np.exp(-(v*1e3- vhalfm)/8.3))
    taum = (mtfactor*qt*taum_exp + 0*v)
    hinf = 1/(1+np.exp( (v*1e3-vhalfh)/9.2))
    tauh = (htfactor*qt*53 + 0*v)


    xgate = moose.element( Ca_R.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum*1e3
    xgate.tableB = 1.0/taum*1e3

    ygate = moose.element( Ca_R.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh*1e3
    ygate.tableB = 1.0/tauh*1e3

    addmsg1 = moose.Mstring( Ca_R.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return Ca_R
