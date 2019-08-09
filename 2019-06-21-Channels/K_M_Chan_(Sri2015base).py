# K_M channel taken from mod files of Srikanth2015 (check once again)

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

def K_M_Chan(name):
    K_M = moose.HHChannel( '/library/' + name )
    K_M.Ek = EK
    K_M.Gbar = 300.0*SOMA_A
    K_M.Gk = 0.0
    K_M.Xpower = 1.0
    K_M.Ypower = 0.0
    K_M.Zpower = 0.0

    vhalfl=-40
    kl=-10
    vhalft=-42
    a0t=0.009
    zetat=7
    gmt=.4
    q10=5
    b0=60
    st=1
    tfactor = 1

    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))
    qt=q10**((celsius-35)/10)
    inf = (1/(1 + np.exp((v*1e3-vhalfl)/kl)))
    a = alpt
    tau = tfactor*(b0 + bett)/(a0t*(1+a))

    xgate = moose.element( K_M.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = inf/tau*1e3
    xgate.tableB = 1.0/tau*1e3
    return K_M
