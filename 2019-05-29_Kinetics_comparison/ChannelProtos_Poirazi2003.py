#exec(open('Kinetics_comparisons/ChannelProtos_Poirazi2003.py').read())
#This is the original taken from mod files. No extra temp compensation etc.
#In K_M_Chan, temperature adjustment is done in a wierd way.
#Author: Anal Kumar
#Defines the channel kinetics

import moose
import rdesigneur as rd
import numpy as np
import csv

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
celsius = 32
dt = 0.05e-3
ENa = 0.050
EK = -0.080 #-0.077 for K_DR
Eh = -0.010
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

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 2.0
    Na.Ypower = 1.0
    Na.Zpower = 1

    a0r = 0.0003
    b0r = 0.0003
    zetar = 12
	zetas = 12
    gmr = 0.2
    ar2 = 1.0
    taumin = 3
    vvs  = 2
    vhalfr = -60
    W = 0.016

    alpv = (1+ar2*np.exp((v*1e3-vhalfr)/vvs))/(1+np.exp((v*1e3-vhalfr)/vvs))
    alpr = np.exp(1.e-3*zetar*(v*1e3-vhalfr)*9.648e4/(8.315*(273.16+celsius)))
    betr = np.exp(1.e-3*zetar*gmr*(v*1e3-vhalfr)*9.648e4/(8.315*(273.16+celsius)))

    minf = 1 / (1 + np.exp((v*1e3 + 44)/(-3)))
    mtau = 0.05 + 0*v*1e3

    hinf = 1 / (1 + exp((v*1e3 + 49)/(3.5)))
    htau = 1 + 0*v*1e3

    sinf = alpv
    taus = betr/(a0r+b0r*alpr)
    taus[taus<taumin] = taumin

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

def K_DR_Chan(name):
    K_DR = moose.HHChannel( '/library/' + name )
    K_DR.Ek = -0.077
    K_DR.Gbar = 300.0*SOMA_A
    K_DR.Gk = 0.0
    K_DR.Xpower = 2.0
    K_DR.Ypower = 0.0
    K_DR.Zpower = 0.0

    ninf = 1 / (1 + exp((v*1e3 + 46.3)/(-3)))
    taun = 3.5 + 0*v*1e3

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

    qt = q10**((celsius-24)/10)
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    alpn = np.exp(1.e-3*zeta*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    betn = np.exp(1.e-3*zeta*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1 + alpn)
    taun = betn/(qt*a0n*(1+alpn))
    taun[taun<nmin] = nmin

    alpl = np.exp(1.e-3*zetal*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    betl = np.exp(1.e-3*zetal*gml*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    linf = 1/(1+ alpl)
    taul = 0.26*(v*1e3+50)
    taul[v<lmin] =lmin

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

	tha  = -30
    qa   = 9
    Ra   = 0.001
    Rb   = 0.001
    q10  = 2.3

    tadj = q10**((celsius - 23)/10)
    a = Ra * (v*1e3 - tha) / (1 - np.exp(-(v*1e3 - tha)/qa))
    b = -Rb * (v*1e3 - tha) / (1 - np.exp((v*1e3 - tha)/qa))
    taun = 1/(a+b)
    inf = a*taun
    tau = 1/(a+b)/tadj # As inf should not be affected by temp

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

    K      = 8.5
    vhalf  = -90

    taul = 2*(1/(np.exp((v*1e3+145)/-17.5)+np.exp((v*1e3+16.8)/16.5)) + 5)
    taul = np.array(taul)
    v = np.array(v)
    taul[v*1e3>-30] = 1
    linf = 1 - (1 / (1 + np.exp((vhalf - v*1e3)/K)))

    xgate = moose.element( h.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = linf/taul*1e3
    xgate.tableB = 1.0/taul*1e3
    return h

def Ca_T_Chan(name):
    Ca_T = moose.HHChannel( '/library/' + name )
    Ca_T.Ek = ECa
    Ca_T.Gbar = 300.0*SOMA_A
    Ca_T.Gk = 0.0
    Ca_T.Xpower = 2.0
    Ca_T.Ypower = 1.0
    Ca_T.Zpower = 1.0
    Ca_L.useConcentration = 1
    Ca_L.instant = 4

    ki = 0.001
    cai = 5.e-5
    cao = 2
    tfa = 1
    tfi = 0.68


    alph = 1.6e-4*np.exp(-(v*1e3+57)/19)
    beth = 1/(np.exp((-v*1e3+15)/10)+1.0)
    alpm = 0.1967*(-1.0*v*1e3+19.88)/(np.exp((-1.0*v*1e3+19.88)/10.0)-1.0)
    betm = 0.046*np.exp(-v*1e3/22.73)
    a = alpm
    taum = 1/(tfa*(a + betm))
    minf =  a/(a+betm)
    a = alph
    tauh = 1/(tfi*(a + beth))
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
    zgate.tableA = ki/(ki+cai)
    zgate.tableB = ki/(ki+cai)

    addmsg4 = moose.Mstring( Ca_L.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'

    #Does not update calcium current in Poirazi's model
    return Ca_T

def Ca_R_Chan(name):
    Ca_R = moose.HHChannel( '/library/' + name )
    Ca_R.Ek = ECa
    Ca_R.Gbar = 300.0*SOMA_A
    Ca_R.Gk = 0.0
    Ca_R.Xpower = 3.0
    Ca_R.Ypower = 1.0
    Ca_R.Zpower = 0.0

    minf = 1 / (1 + np.exp((v*1e3+60)/(-3)))
    taum = 100 + 0*v*1e3

    hinf = 1/ (1 + np.exp((v*1e3+62)/(1)))
    tauh = 5 + 0*v*1e3

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
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( Ca_L.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( Ca_L.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return Ca_L

def K_SK_Chan(name):
    K_SK = moose.HHChannel( '/library/' + name )
    K_SK.Ek = EK
    K_SK.Gbar = 300.0*SOMA_A
    K_SK.Gk = 0.0
    K_SK.Xpower = 0.0
    K_SK.Ypower = 0.0
    K_SK.Zpower = 3.0
    K_SK.useConcentration = 1

    cai     = 2.4e-5
    beta    = 0.03
    cac     = 0.025
    taumin  = 0.5

    tadj = 3**((celsius-22.0)/10)
    car = (cai/cac)**2
    m_inf = car / ( 1 + car )
    tau_m =  1 / beta / (1 + car) / tadj
    tau_m[tau_m<taumin] = taumin

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

    cai = 1e-3
    d1 = 0.84
    d2 = 1.0
    k1 = 0.18
    k2 = 0.011
    bbar = 0.28
    abar = 0.48
    F_KC = F/1000.0

    def exp1(k,d,v):
        return k*np.exp(-2*d*F_KC*v*1e3/R/(273.15 + celsius))
    def alp(v,c):
        return c*abar/(c + exp1(k1,d1,v))
    def bet(v,c):
        return bbar/(1 + c/exp1(k2,d2,v))

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
        tblA[i] = np.array([alp(x,y) for y in ca])
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
    depth    = 0.1e-6
    taur    = 200e-3
    cainf    = 100e-6
    # B = 28789637.7

    Ca.tau = taur/7
    Ca.Ca_base = cainf
    Ca.diameter = 500e-6
    Ca.length = 1000e-6
    Ca.thick = 177.9e-6
    return Ca
