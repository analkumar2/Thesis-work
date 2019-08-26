# Na channel taken from mod files of Migliore2018: na3n.mod
# Problems: q10 same for both X and Y gates.

import numpy as np
import pickle
import pandas as pd
import moose

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.092
EK = -0.100
Eh = -0.030
ECa = 0.140
Em = -0.065

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 1e-12
Camax = 1
Cadivs = 400
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 3.0
    Na.Ypower = 1.0
    Na.Zpower = 1

    sh2   = 0
    gbar = 0.010
    tha  =  -30
    qa   = 7.2
    Ra   = 0.4
    Rb   = 0.124
    thi1  = -45
    thi2  = -45
    qd   = 1.5
    qg   = 1.5
    mmin=0.02
    hmin=0.5
    q10=2
    Rg   = 0.01
    Rd   = .03
    qq   = 10
    tq   = -55
    thinf  = -50
    qinf  = 4
    vhalfs=-60
    a0s=0.0003
    zetas=12
    gms=0.2
    smax=10
    vvh=-58
    vvs=2
    a2=1

    def trap0(v,th,a,q):
        if np.abs(v*1e3-th) > 1e-6:
            return a * (v*1e3 - th) / (1 - np.exp(-(v*1e3 - th)/q))
        else:
            return a * q

    qt=q10**((celsius-24)/10)
    a = np.array([trap0(vm,tha+sh2,Ra,qa) for vm in v])
    b = np.array([trap0(-vm,-tha-sh2,Rb,qa) for vm in v])
    mtau = 1/(a+b)/qt
    mtau[mtau<mmin] = mmin
    minf = a/(a+b)

    a = np.array([trap0(vm,thi1+sh2,Rd,qd) for vm in v])
    b = np.array([trap0(-vm,-thi2-sh2,Rg,qg) for vm in v])
    htau =  1/(a+b)/qt
    htau[htau<hmin] = hmin
    hinf = 1/(1+np.exp((v*1e3-thinf-sh2)/qinf))

    c = 1/(1+np.exp((v*1e3-vvh-sh2)/vvs))
    sinf = c+a2*(1-c)
    alps = np.exp(1.e-3*zetas*(v*1e3-vhalfs-sh2)*9.648e4/(8.315*(273.16+celsius)))
    bets = np.exp(1.e-3*zetas*gms*(v*1e3-vhalfs-sh2)*9.648e4/(8.315*(273.16+celsius)))
    taus = bets/(a0s*(1+alps))
    taus[taus<smax] = smax

    xgate = moose.element( Na.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau*1e3
    xgate.tableB = 1.0/mtau*1e3

    ygate = moose.element( Na.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/htau*1e3
    ygate.tableB = 1.0/htau*1e3

    zgate = moose.element( Na.path + '/gateZ' )
    zgate.min = Vmin
    zgate.max = Vmax
    zgate.divs = Vdivs
    zgate.tableA = sinf/taus*1e3
    zgate.tableB = 1.0/taus*1e3
    return Na
