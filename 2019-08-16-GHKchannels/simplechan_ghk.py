# Simple channel simplechan_ghk.mod

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
Vdivs = 10
dV = (Vmax-Vmin)/Vdivs
v = np.arange(Vmin,Vmax+dV, dV)
Camin = 1e-12
Camax = 3e-3
Cadivs = 20
dCa = (Camax-Camin)/Cadivs
ca = np.arange(Camin,Camax+dCa, dCa)

def simp_Chan(name):
    simp = moose.HHChannel2D( '/library/' + name )
    simp.Ek = EK
    simp.Gbar = 300.0*SOMA_A
    simp.Gk = 0.0
    simp.Xpower = 1.0
    simp.Ypower = 0.0
    simp.Zpower = 0.0
    simp.Xindex = 'VOLT_INDEX'
    # simp.Yindex = 'VOLT_INDEX'
    # simp.Zindex = 'VOLT_C1_INDEX'

    xgate = moose.element( simp.path + '/gateX' )
    xgate.xminA = Vmin
    xgate.xmaxA = Vmax
    xgate.xdivsA = Vdivs
    xgate.yminA = Vmin
    xgate.ymaxA = Vmax
    xgate.ydivsA = Vdivs
    xgate.xminB = Vmin
    xgate.xmaxB = Vmax
    xgate.xdivsB = Vdivs
    xgate.yminB = Vmin
    xgate.ymaxB = Vmax
    xgate.ydivsB = Vdivs

    # ygate = moose.element( simp.path + '/gateY' )
    # ygate.xminA = 0
    # ygate.xmaxA = 0
    # ygate.xdivsA = 0
    # ygate.yminA = Vmin
    # ygate.ymaxA = Vmax
    # ygate.ydivsA = Vdivs
    # ygate.xminB = 0
    # ygate.xmaxB = 0
    # ygate.xdivsB = 0
    # ygate.yminB = Vmin
    # ygate.ymaxB = Vmax
    # ygate.ydivsB = Vdivs

    # zgate = moose.element( simp.path + '/gateZ' )
    # zgate.xminA = Vmin
    # zgate.xmaxA = Vmax
    # zgate.xdivsA = Vdivs
    # zgate.yminA = Camin
    # zgate.ymaxA = Camax
    # zgate.ydivsA = Cadivs
    # zgate.xminB = Vmin
    # zgate.xmaxB = Vmax
    # zgate.xdivsB = Vdivs
    # zgate.yminB = Camin
    # zgate.ymaxB = Camax
    # zgate.ydivsB = Cadivs

    cai = 50.e-6
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
    x = Vmin
    tblA = np.zeros([2,xgate.ydivsA])
    tblB = np.zeros([2,len(v)])
    for i in np.arange(2):
        tblA[i] = minf/mtau*np.ones(len(v))
        tblB[i] = 1/mtau*np.ones(len(v))
        x = x + dV
        print(i, end='\r')
    xgate.tableA = tblA*1e3
    xgate.tableB = tblB*1e3

    # x = Vmin
    # tblA = np.zeros([len(v),len(ca)])
    # tblB = np.zeros([len(v),len(ca)])
    # for i in np.arange(len(v)):
    #     tblA[i] = np.array([alp(x,y) for y in ca])
    #     tblB[i] = np.add(tblA[i],[bet(x, y) for y in ca])/tfactor
    #     x = x + dV
    #     print(i, end='\r')
    #
    # zgate.tableA = tblA*1e3
    # zgate.tableB = tblB*1e3

    addmsg4 = moose.Mstring( simp.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return simp
