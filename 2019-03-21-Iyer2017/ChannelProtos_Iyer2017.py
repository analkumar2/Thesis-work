#Author: Anal Kumar
#Defines the channel kinetics

import moose
import rdesigneur as rd
import numpy as np
import csv

SOMA_A = 4*np.pi*0.5e-6**2
F = 96485.3329
R = 8.314
Temp = 307.15
dt = 0.05e-3
ENa = 0.055
EK = -0.075 #EK for KDR is set to -0.077
ECa = 0.50
Em = -0.045
Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
dV = (Vmax-Vmin)/Vdivs
Camin = 0
Camax = 3e-3
Cadivs = 3000
dCa = (Camax-Camin)/Cadivs


def Na_SChan(name):
    Na_S = moose.HHChannel( '/library/' + name )
    Na_S.Ek = ENa
    Na_S.Gbar = 300.0*SOMA_A
    Na_S.Gk = 0.0
    Na_S.Xpower = 3.0
    Na_S.Ypower = 1.0
    Na_S.Zpower = 1.0

    tha  =  -25e-3
    qa   = 7.2e-3
    Ra   = 0.4e3
    Rb   = 0.124e3
    thi1  = -45e-3
    thi2  = -45e-3
    qd   = 1.5e-3
    qg   = 1.5e-3
    mmin=0.02
    hmin=0.5
    q10=3
    Rg   = 0.01e3
    Rd   = 0.03e3
    qq   = 10e-3
    tq   = -55e-3
    thinf  = -50e-3
    qinf  = 2e-3
    vhalfs=-60e-3
    a0s=0.0003e3
    zetas=12
    gms = 0.2
    smax=10e-3
    vvh=-58e-3
    vvs=2e-3
    ar2=1
    celsius = 34

    v = np.arange(Vmin,Vmax+dV, dV)
    alpv = 1/(1+np.exp((v-vvh)/vvs))
    alps = np.exp(zetas*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    bets = np.exp(zetas*gms*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    qt=q10**((celsius-24)/10)

    a = Ra * (v - tha) / (1 - np.exp(-(v - tha)/qa))
    a[abs(v-tha)<=1e-9] = Ra*qa
    b = Rb * (-v + tha) / (1 - np.exp(-(-v + tha)/qa))
    b[abs(v-tha)<=1e-9] = Rb*qa
    mtau = 1e-3/(a+b)/qt
    mtau[mtau<mmin*1e-3]=mmin*1e-3
    minf = a/(a+b)

    a = Rd * (v - thi1) / (1 - np.exp(-(v - thi1)/qd))
    a[abs(v-tha)<=1e-6] = Rd*qd
    b = Rg * (-v + thi2) / (1 - np.exp(-(-v + thi2)/qg))
    b[abs(v-tha)<=1e-6] = Rb*qa
    htau = 1e-3/(a+b)/qt
    htau[htau<hmin*1e-3]=hmin*1e-3
    hinf = 1/(1+np.exp((v-thinf)/qinf))

    c=alpv
    sinf = c+ar2*(1-c)
    taus = bets/(a0s*(1+alps))
    taus[taus<smax*1e-3]=smax*1e-3

    xgate = moose.element( Na_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau
    xgate.tableB = 1.0/mtau

    ygate = moose.element( Na_S.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/htau
    ygate.tableB = 1.0/htau

    zgate = moose.element( Na_S.path + '/gateZ' )
    zgate.min = Vmin
    zgate.max = Vmax
    zgate.divs = Vdivs
    zgate.tableA = sinf/taus
    zgate.tableB = 1.0/taus
    return Na_S

def KDR_SChan(name):
    KDR_S = moose.HHChannel( '/library/' + name )
    KDR_S.Ek = -0.072
    KDR_S.Gbar = 300.0*SOMA_A
    KDR_S.Gk = 0.0
    KDR_S.Xpower = 1.0
    KDR_S.Ypower = 0.0
    KDR_S.Zpower = 0.0

    vhalfn=13e-3
    a0n=0.02e3
    zetan=-3
    gmn=0.7
    nmax=2
    q10=1
    celsius = 34

    v = np.arange(Vmin,Vmax+dV, dV)
    alpn = np.exp(zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    qt=q10**((celsius-24)/10)
    a = alpn
    ninf = 1/(1+a)
    taun = betn/(qt*a0n*(1+a))
    taun[taun<nmax*1e-3]=nmax*1e-3

    xgate = moose.element( KDR_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun
    xgate.tableB = 1/taun
    return KDR_S

def KA_SChan(name):
    KA_S = moose.HHChannel( '/library/' + name )
    KA_S.Ek = EK
    KA_S.Gbar = 300.0*SOMA_A
    KA_S.Gk = 0.0
    KA_S.Xpower = 4.0
    KA_S.Ypower = 1.0
    KA_S.Zpower = 0.0

    v = np.arange(Vmin,Vmax+dV, dV)
    ninf = 1/(1+np.exp((60e-3-v-42e-3)/15e-3))
    taun = 10e-3 + 0*v
    linf = 1/(1+np.exp(-(v+43e-3)/20e-3))
    taul = 1.1e-3 + 2e-3*np.exp(-(v*1e3+50)**2/50)

    xgate = moose.element( KA_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun
    xgate.tableB = 1.0/taun

    ygate = moose.element( KA_S.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = linf/taul
    ygate.tableB = 1.0/taul
    return KA_S

def CaL_SChan(name):
    CaL_S = moose.HHChannel( '/library/' + name )
    CaL_S.Ek = ECa
    CaL_S.Gbar = 300.0*SOMA_A
    CaL_S.Gk = 0.0
    CaL_S.Xpower = 1.0
    CaL_S.Ypower = 0.0
    CaL_S.Zpower = 0.0

    v = np.arange(Vmin,Vmax+dV, dV)
    minf = 1/(1+np.exp((-v+55e-3)/5e-3))
    taum = 1.5e-3 + 18e-3*np.exp(-(v*1e3+45)**2/625)

    xgate = moose.element( CaL_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum
    xgate.tableB = 1.0/taum

    return CaL_S

def KSK_Chan(name):
    KSK = moose.HHChannel( '/library/' + name )
    KSK.Ek = EK
    KSK.Gbar = 300.0*SOMA_A
    KSK.Gk = 0.0
    KSK.Xpower = 0.0
    KSK.Ypower = 0.0
    KSK.Zpower = 1.0
    KSK.useConcentration = 1
    KSK.instant = 4

    ca = np.arange(Camin,Camax+dCa, dCa)
    KSK_km = 0.0002
    zgate = moose.element( KSK.path + '/gateZ' )
    zgate.tableA = ca**4/(ca**4+KSK_km**4)
    zgate.tableB = 1 + 0*ca

    addmsg4 = moose.Mstring( KSK.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return KSK

def Ca_Conc ( name ):
    Ca = moose.CaConc( '/library/' + name )
    taur	= 200e-3

    Ca.tau = 7*taur
    return Ca
