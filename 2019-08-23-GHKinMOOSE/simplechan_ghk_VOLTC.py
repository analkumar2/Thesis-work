# Simple channel simplechan_ghk.mod
# exec(open('simplechan_ghk.py').read())

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
# dV = (Vmax-Vmin)/Vdivs
# v = np.arange(Vmin,Vmax, dV)
v = np.linspace(Vmin,Vmax, Vdivs)
Camin = 1e-12
Camax = 3
Cadivs = 4000
# dCa = (Camax-Camin)/Cadivs
# ca = np.arange(Camin,Camax, dCa)
ca = np.linspace(Camin,Camax, Cadivs)

def simp_Chan(name):
    simp = moose.HHChannel2D( '/library/' + name )
    simp.Ek = ECa
    simp.Gbar = 300.0*SOMA_A
    simp.Gk = 0.0
    simp.Xpower = 1.0
    simp.X = 0
    simp.Ypower = 0.0
    simp.Zpower = 1.0
    simp.Xindex = 'VOLT_C1_INDEX'
    # simp.Yindex = 'VOLT_INDEX'
    simp.Zindex = 'VOLT_C1_INDEX'
    simp.instant = 4

    xgate = moose.element( simp.path + '/gateX' )
    xgate.xminA = Vmin
    xgate.xmaxA = Vmax
    xgate.xdivsA = Vdivs
    xgate.yminA = Camin
    xgate.ymaxA = Camax
    xgate.ydivsA = Cadivs
    xgate.xminB = Vmin
    xgate.xmaxB = Vmax
    xgate.xdivsB = Vdivs
    xgate.yminB = Camin
    xgate.ymaxB = Camax
    xgate.ydivsB = Cadivs

    zgate = moose.element( simp.path + '/gateZ' )
    zgate.xminA = Vmin
    zgate.xmaxA = Vmax
    zgate.xdivsA = Vdivs
    zgate.yminA = Camin
    zgate.ymaxA = Camax
    zgate.ydivsA = Cadivs
    zgate.xminB = Vmin
    zgate.xmaxB = Vmax
    zgate.xdivsB = Vdivs
    zgate.yminB = Camin
    zgate.ymaxB = Camax
    zgate.ydivsB = Cadivs

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
    gmm=0.1
    qt=q10**((celsius-25)/10)
    a = 0.2*(-1.0*v*1e3+19.26)/(np.exp((-1.0*v*1e3+19.26)/10.0)-1.0)
    b = 0.009*np.exp(-v*1e3/22.03)
    minf = a/(a+b)
    alpmt = np.exp(0.0378*zetam*(v*1e3-vhalfm))
    betmt = np.exp(0.0378*zetam*gmm*(v*1e3-vhalfm))
    mtau = betmt/(qt*a0m*(1+alpmt))
    mtau[mtau<mmin]=mmin
    print(min(1/mtau))
    tblA = np.zeros([Vdivs,Cadivs])
    tblB = np.zeros([Vdivs,Cadivs])
    for i in np.arange(Cadivs):
        tblA[:,i] = minf/mtau
        tblB[:,i] = 1/mtau
        print(i, end='\r')
    xgate.tableA = tblA*1e3
    xgate.tableB = tblB*1e3


    z = 2
    T = celsius+273.15
    cao = 2
    ezfrt  = np.exp(z*v*F/R/T)
    # ezfrt = np.exp(z*v*1e3*293.15/25/T) #Used in NEURON. gives very slightly different values
    tblA = np.zeros([Vdivs,Cadivs])
    tblB = np.zeros([Vdivs,Cadivs])
    for i in np.arange(Cadivs):
        tblA[:,i] = v*(ca[i]/cao*ezfrt-1)/(ezfrt-1)/(v-ECa)
        print(i, end='\r')
    tblB = tblA*0 +1
    zgate.tableA = tblA
    zgate.tableB = tblB

    addmsg4 = moose.Mstring( simp.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc concOut . concen'

    addmsg2 = moose.Mstring( simp.path + '/addmsg2' )
    addmsg2.value = '. IkOut ../Ca_conc current'
    return simp
