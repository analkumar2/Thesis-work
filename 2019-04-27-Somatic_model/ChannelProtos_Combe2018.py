# exec(open('Somatic model/ChannelProtos_Combe2018.py').read())
#Author: Anal Kumar
#Defines the channel kinetics

import moose
import rdesigneur as rd
import numpy as np
import csv

SOMA_A = 3.32e-9
F = 96485.3329
R = 8.314
Temp = 307.15
dt = 0.05e-3
ENa = 0.050
EK = -0.080 #EK for K_DR is set to -0.077
Eh = -0.010
ECa = 0.140
Em = -0.070
celsius = 34

Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
dV = (Vmax-Vmin)/Vdivs
v = np.arange(Vmin,Vmax+dV, dV)
Camin = 1e-12
Camax = 10e-3
Cadivs = 4000
dCa = (Camax-Camin)/Cadivs
ca = np.arange(Camin,Camax+dCa, dCa)

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 3.0
    Na.Ypower = 1.0
    Na.Zpower = 0.0

    tha  =  -25
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
    Rd   = 0.03
    qq   = 10
    tq   = -55
    thinf  = -50
    qinf  = 1
    vhalfs=-60
    a0s=0.0003
    zetas=12
    gms = 0.2
    smax=10
    vvh=-58
    vvs=2
    ar2=1
    celsius = 34

    alpv = 1/(1+np.exp((v*1e3-vvh)/vvs))
    alps = np.exp(1.e-3*zetas*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    bets = np.exp(1.e-3*zetas*gms*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    qt=q10**((celsius-24)/10)

    a = Ra * (v*1e3 - tha) / (1 - np.exp(-(v*1e3 - tha)/qa))
    a[abs(v*1e3-tha)<=1e-6] = Ra*qa
    b = Rb * (-v*1e3 + tha) / (1 - np.exp(-(-v*1e3 + tha)/qa))
    b[abs(-v*1e3+tha)<=1e-6] = Rb*qa
    mtau = 1/(a+b)/qt
    mtau[mtau<mmin]=mmin
    minf = a/(a+b)

    a = Rd * (v*1e3 - thi1) / (1 - np.exp(-(v*1e3 - thi1)/qd))
    a[abs(v*1e3-tha)<=1e-6] = Rd*qd
    b = Rg * (-v*1e3 + thi2) / (1 - np.exp(-(-v*1e3 + thi2)/qg))
    b[abs(-v*1e3+tha)<=1e-6] = Rg*qg
    htau = 1/(a+b)/qt
    htau[htau<hmin]=hmin
    hinf = 1/(1+np.exp((v*1e3-thinf)/qinf))

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
    return Na

def K_DR_Chan(name):
    K_DR = moose.HHChannel( '/library/' + name )
    K_DR.Ek = -0.077
    K_DR.Gbar = 300.0*SOMA_A
    K_DR.Gk = 0.0
    K_DR.Xpower = 1.0
    K_DR.Ypower = 0.0
    K_DR.Zpower = 0.0

    vhalfn=13
    a0n=0.02
    zetan=-3
    gmn=0.7
    nmax=0.1
    q10=1
    celsius = 34

    alpn = np.exp(1.e-3*zetan*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(1.e-3*zetan*gmn*(v*1e3-13)*9.648e4/(8.315*(273.16+celsius)))
    qt=q10**((celsius-24)/10)
    a = alpn
    ninf = 1/(1+a)
    c = np.exp(-339.44*(v*1e3-13)/(8.315*(273.16+celsius)))
    taun = 0.7*betn/(qt*a0n*(1+c))
    taun[taun<nmax]=nmax

    xgate = moose.element( K_DR.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun*1e3
    xgate.tableB = 1/taun*1e3
    return K_DR

def K_A_Chan(name):
    K_A = moose.HHChannel( '/library/' + name )
    K_A.Ek = EK
    K_A.Gbar = 300.0*SOMA_A
    K_A.Gk = 0.0
    K_A.Xpower = 1.0
    K_A.Ypower = 1.0
    K_A.Zpower = 0.0

    vhalfn=11
    vhalfl=-56
    a0l=0.05
    a0n=0.05
    zetan=-1.5
    zetal=3
    gmn=0.55
    gml=1
    lmin=2
    nmin=0.1
    pw=-1
    tq=-40
    qq=5
    q10=5
    qtl=1

    qt = q10**((celsius-24)/10)
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    alpn = np.exp(1.e-3*zeta*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(1.e-3*zeta*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    alpl = np.exp(1.e-3*zetal*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    betl = np.exp(1.e-3*zetal*gml*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))

    a = alpn
    ninf = 1.0/(1 + a)
    taun = betn/(qt*a0n*(1+a))
    taun[taun<nmin] = nmin

    a = alpl
    linf = 1/(1+ a)
    taul = 0.26*(v*1e3+50)/qtl
    taul[taul<lmin/qtl] = lmin/qtl

    xgate = moose.element( K_A.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun*1e3
    xgate.tableB = 1.0/taun*1e3

    ygate = moose.element( K_A.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = linf/taul*1e3
    ygate.tableB = 1.0/taul*1e3
    return K_A

def K_M_Chan(name):
    K_M = moose.HHChannel( '/library/' + name )
    K_M.Ek = EK
    K_M.Gbar = 300.0*SOMA_A
    K_M.Gk = 0.0
    K_M.Xpower = 1.0
    K_M.Ypower = 0.0
    K_M.Zpower = 0.0

    vhalfl=-42
    kl=-4
    vhalft=-42
    a0t=0.04
    zetat=4
    gmt=0.7
    q10=5
    b0=60
    celsius = 34

    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))

    qt=q10**((celsius-35)/10)
    inf = (1.0/(1 + np.exp((v*1e3-vhalfl)/kl)))
    a = alpt
    tau = b0 + bett/(a0t*(1+a))

    xgate = moose.element( K_M.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = inf/tau*1e3
    xgate.tableB = 1.0/tau*1e3
    return K_M

def h_Chan(name):
    h = moose.HHChannel( '/library/' + name )
    h.Ek = Eh
    h.Gbar = 300.0*SOMA_A
    h.Gk = 0.0
    h.Xpower = 1.0
    h.Ypower = 0.0
    h.Zpower = 0.0

    ena    = 50
    K      = 8.5
    vhalf  = -81

    taun = 5*(1/(np.exp((v*1e3+145)/-17.5)+np.exp((v*1e3+16.8)/16.5)) + 5)
    taun[v*1e3>-30] = 1
    ninf = 1 - (1 / (1 + np.exp((vhalf - v*1e3)/K)))

    xgate = moose.element( h.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun*1e3
    xgate.tableB = 1.0/taun*1e3
    return h

def Ca_T_Chan(name):
    Ca_T = moose.HHChannel( '/library/' + name )
    Ca_T.Ek = ECa
    Ca_T.Gbar = 300.0*SOMA_A
    Ca_T.Gk = 0.0
    Ca_T.Xpower = 2.0
    Ca_T.Ypower = 1.0
    Ca_T.Zpower = 1.0
    Ca_T.useConcentration = 1
    Ca_T.instant = 4

    tBase = 23.5
    celsius = 22
    gcatbar = 0
    ki = 0.001
    cai = 5e-5
    cao = 2
    tfa = 1
    tfi = 0.68

    alph = 1.6e-4*np.exp(-(v*1e3+57)/19)
    beth = 1/(np.exp((-v*1e3+15)/10)+1.0)
    alpm = 0.1967*(-1.0*v*1e3+19.88)/(np.exp((-1.0*v*1e3+19.88)/10.0)-1.0)
    betm = 0.046*np.exp(-v*1e3/22.73)

    a = alpm
    taum = 1.0/(tfa*(a + betm))
    minf =  a/(a+betm)
    a = alph
    tauh = 1.0/(tfi*(a + beth))
    hinf = a/(a+beth)

    xgate = moose.element( Ca_T.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum*1e3
    xgate.tableB = 1.0/taum*1e3

    ygate = moose.element( Ca_T.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh*1e3
    ygate.tableB = 1.0/tauh*1e3

    zgate = moose.element( Ca_T.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    ca = np.arange(Camin,Camax+dCa, dCa)
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    #Does not activate calcium-dependent currents
    addmsg4 = moose.Mstring( Ca_T.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return Ca_T

def Ca_R_Chan(name):
    Ca_R = moose.HHChannel( '/library/' + name )
    Ca_R.Ek = ECa
    Ca_R.Gbar = 300.0*SOMA_A
    Ca_R.Gk = 0.0
    Ca_R.Xpower = 3.0
    Ca_R.Ypower = 1.0
    Ca_R.Zpower = 0.0

    celsius =34

    infm = 1 / (1 + np.exp((v*1e3+60)/(-3)))
    taum = 100+0*v
    infh = 1/ (1 + np.exp((v*1e3+62)/(1)))
    tauh = 5+0*v

    xgate = moose.element( Ca_R.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = infm/taum*1e3
    xgate.tableB = 1.0/taum*1e3

    ygate = moose.element( Ca_R.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = infh/tauh*1e3
    ygate.tableB = 1.0/tauh*1e3

    addmsg1 = moose.Mstring( Ca_R.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return Ca_R

def CaL_Chan(name):
    Ca_L = moose.HHChannel( '/library/' + name )
    Ca_L.Ek = ECa
    Ca_L.Gbar = 300.0*SOMA_A
    Ca_L.Gk = 0.0
    Ca_L.Xpower = 1.0
    Ca_L.Ypower = 0.0
    Ca_L.Zpower = 1.0
    Ca_L.useConcentration = 1
    Ca_L.instant = 4

    celsius= 34
    ki=.001
    cai = 50e-6
    cao = 2
    tfa = 5

    alpm = 0.055*(-27.01 - v*1e3)/(np.exp((-27.01-v*1e3)/3.8) - 1)
    betm =0.94*np.exp((-63.01-v*1e3)/17)
    a = alpm
    taum = 1/(tfa*(a+betm))
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
    ca = np.arange(Camin,Camax+dCa, dCa)
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( Ca_L.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( Ca_L.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return Ca_L

# #Modified to fit data from Hirschberg et. al., 1999.
# def K_SK_Chan(name):
#     K_SK = moose.HHChannel( '/library/' + name )
#     K_SK.Ek = EK
#     K_SK.Gbar = 300.0*SOMA_A
#     K_SK.Gk = 0.0
#     K_SK.Xpower = 0.0
#     K_SK.Ypower = 0.0
#     K_SK.Zpower = 1.0
#     K_SK.useConcentration = 1
#
#     celsius = 36
#     cac     = 0.56e-3
#     taumin  = 6.3e-3
#
#     ca = np.arange(Camin,Camax+dCa, dCa)
#     car = (ca/cac)**4.6
#     m_inf = car / ( 1 + car )
#     tau_m =  taumin + 0*ca
#
#     zgate = moose.element( K_SK.path + '/gateZ' )
#     zgate.min = Camin
#     zgate.max = Camax
#     zgate.divs = Cadivs
#     zgate.tableA = m_inf/tau_m
#     zgate.tableB = 1.0/tau_m
#
#     addmsg3 = moose.Mstring( K_SK.path + '/addmsg3' )
#     addmsg3.value = '../Ca_conc    concOut    . concen'
#     return K_SK

# Used in Combe2018
def K_SK_Chan(name):
    K_SK = moose.HHChannel( '/library/' + name )
    K_SK.Ek = EK
    K_SK.Gbar = 300.0*SOMA_A
    K_SK.Gk = 0.0
    K_SK.Xpower = 0.0
    K_SK.Ypower = 0.0
    K_SK.Zpower = 3.0
    K_SK.useConcentration = 1

    celsius = 36
    beta    = 0.03
    cac     = 0.025
    taumin  = 5

    tadj = 3**((celsius-22.0)/10)
    car = (ca/cac)**2
    m_inf = car / ( 1 + car )
    tau_m =  1 / beta / (1 + car) / tadj
    tau_m[tau_m < taumin] = taumin

    zgate = moose.element( K_SK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = m_inf/tau_m*1e3
    zgate.tableB = 1.0/tau_m*1e3

    addmsg3 = moose.Mstring( K_SK.path + '/addmsg3' )
    addmsg3.value = '../Ca_conc    concOut    . concen'
    return K_SK

def K_BK_Chan(name):
    K_BK = moose.HHChannel2D( '/library/' + name )
    K_BK.Ek = EK
    K_BK.Gbar = 300.0*SOMA_A
    K_BK.Gk = 0.0
    K_BK.Xpower = 0.0
    K_BK.Ypower = 0.0
    K_BK.Zpower = 1.0
    K_BK.Zindex = 'VOLT_C1_INDEX'

    celsius = 20
    gkbar = 0.01e4
    d1 =1
    d2 = 1.5
    k1 = 0.18
    k2 = 0.011
    bbar = 0.28
    abar = 0.48
    F_KC = F/1000.0

    def alp(v,c):
        return c*abar/(c + k1*np.exp(-2*d1*F_KC*v*1e3/R/(273.15 + celsius)))
    def bet(v,c):
        return bbar/(1 + c/(k2*np.exp(-2*d2*F*v/R/(273.15 + celsius))))

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

    x = Vmin
    tblA = np.zeros([len(v),len(ca)])
    tblB = np.zeros([len(v),len(ca)])
    for i in np.arange(len(v)):
        tblA[i] = [alp(x,y) for y in ca]
        tblB[i] = np.add(tblA[i],[bet(x, y) for y in ca])
        x = x + dV
        print(i, end='\r')

    zgate.tableA = tblA*1e3
    zgate.tableB = tblB*1e3

    addmsg4 = moose.Mstring( K_BK.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return K_BK

def Ca_Conc ( name ):
    Ca = moose.CaConc( '/library/' + name )
    depth    = 0.5e-6
    taur    = 1400e-3
    cainf    = 100e-6
    # B = 28789637.7
    Ca.tau = taur
    Ca.Ca_base = cainf
    Ca.diameter = 100e-6
    Ca.length = 100e-6
    Ca.thick = 100e-6
    return Ca
