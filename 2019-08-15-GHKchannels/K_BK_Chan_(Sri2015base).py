# K_BK channel taken from mod files of Srikanth2015 (check once again)
# Problems: qt not proper

import numpy as np
import pickle
import pandas as pd
import moose
import os

print(os.path.realpath(__file__))

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
Cadivs = 4000
dCa = (Camax-Camin)/Cadivs
ca = np.arange(Camin,Camax+dCa, dCa)

with open('Channels/K_BK_Chan_tableA.txt', 'rb') as f:
    tblA = pickle.load(f)
with open('Channels/K_BK_Chan_tableB.txt', 'rb') as f:
    tblB = pickle.load(f)

def K_BK_Chan(name):
    K_BK = moose.HHChannel2D( '/library/' + name )
    K_BK.Ek = EK
    K_BK.Gbar = 300.0*SOMA_A
    K_BK.Gk = 0.0
    K_BK.Xpower = 0.0
    K_BK.Ypower = 0.0
    K_BK.Zpower = 1.0
    K_BK.Zindex = 'VOLT_C1_INDEX'

    # cai = 5.e-5
    # d1 = .84
    # d2 = 1.
    # k1 = .48e-3
    # k2 = .13e-6
    # abar = .28
    # bbar = .48
    # st=1
    # tfactor = 1
    # F_KC = F/1000.0
    #
    # def exp1(k,d,v):
    #     return k*np.exp(-2*d*F_KC*v*1e3/R/(273.15 + celsius))
    # def alp(v,c):
    #     return c*abar/(c + exp1(k1,d1,v))
    # def bet(v,c):
    #     return bbar/(1 + c/exp1(k2,d2,v))

    zgate = moose.element( K_BK.path + '/gateZ' )
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

    zgate.tableA = tblA
    zgate.tableB = tblB

    addmsg4 = moose.Mstring( K_BK.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return K_BK
