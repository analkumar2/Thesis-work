#Author: Anal Kumar
#Defines the channel kinetics
# Delete from kinetics_comparison folder

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
EK = -0.080 #EK for KDR is set to -0.077
Eh = -0.010
ECa = 0.140
Em = -0.075
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
    q10=2
    Rg   = 0.01e3
    Rd   = 0.03e3
    qq   = 10e-3
    tq   = -55e-3
    thinf  = -50e-3
    qinf  = 1e-3
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
    KDR_S.Ek = -0.077
    KDR_S.Gbar = 300.0*SOMA_A
    KDR_S.Gk = 0.0
    KDR_S.Xpower = 1.0
    KDR_S.Ypower = 0.0
    KDR_S.Zpower = 0.0

    vhalfn=13e-3
    a0n=0.02e3
    zetan=-3
    gmn=0.7
    nmax=0.1
    q10=1
    celsius = 34

    v = np.arange(Vmin,Vmax+dV, dV)
    alpn = np.exp(zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    qt=q10**((celsius-24)/10)
    a = alpn
    ninf = 1/(1+a)
    c = np.exp(-339.44*(v*1e3-13)/(8.315*(273.16+celsius)))
    taun = 0.7*betn/(qt*a0n*(1+c))
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
    KA_S.Xpower = 1.0
    KA_S.Ypower = 1.0
    KA_S.Zpower = 0.0

    gkabar = 0
    q10 = 5
    qtl = 1
    a0l = 0.05e3
    a0n = 0.05e3
    nmin = 0.1e-3
    lmin = 2e-3
    celsius = 34
    gmn = 0.55
    gml = 1
    zetan = -1.5
    zetal = 3
    pw = -1
    tq = -40e-3
    qq = 5e-3
    vhalfl = -56e-3
    vhalfn = 11e-3

    v = np.arange(Vmin,Vmax+dV, dV)
    qt = q10**((celsius-24)/10)
    zeta = zetan+pw/(1+np.exp((v-tq)/qq))
    alpn = np.exp(zeta*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    alpl = np.exp(zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    betl = np.exp(zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius)))

    a = alpn
    ninf = 1.0/(1 + a)
    taun = betn/(qt*a0n*(1+a))
    taun[taun<nmin] = nmin

    a = alpl
    linf = 1/(1+ a)
    taul = 0.26*(v+50e-3)/qtl
    taul[taul<lmin/qtl] = lmin/qtl

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

def KM_Chan(name):
    KM = moose.HHChannel( '/library/' + name )
    KM.Ek = EK
    KM.Gbar = 300.0*SOMA_A
    KM.Gk = 0.0
    KM.Xpower = 1.0
    KM.Ypower = 0.0
    KM.Zpower = 0.0

    vhalfl=-42e-3
    kl=-4e-3
    vhalft=-42e-3
    a0t=0.04e3
    zetat=4
    gmt=0.7
    q10=5
    b0=60e-3
    celsius = 34

    v = np.arange(Vmin,Vmax+dV, dV)
    alpt = np.exp(0.0378e3*zetat*(v-vhalft))
    bett = np.exp(0.0378e3*zetat*gmt*(v-vhalft))

    qt=q10**((celsius-35)/10)
    inf = (1.0/(1 + np.exp((v-vhalfl)/kl)))
    a = alpt
    tau = b0 + bett/(a0t*(1+a))

    xgate = moose.element( KM.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = inf/tau
    xgate.tableB = 1.0/tau
    return KM

def h_Chan(name):
    h = moose.HHChannel( '/library/' + name )
    h.Ek = Eh
    h.Gbar = 300.0*SOMA_A
    h.Gk = 0.0
    h.Xpower = 1.0
    h.Ypower = 0.0
    h.Zpower = 0.0

    ena    = 50e-3
    K      = 8.5e-3
    vhalf  = -81e-3

    v = np.arange(Vmin,Vmax+dV, dV)
    taun = 1e-3+v*0
    taun[v<=-30e-3] = 5e-3*(1/(np.exp((v[v<=-30e-3]+145e-3)/-17.5e-3)+np.exp((v[v<=-30e-3]+16.8e-3)/16.5e-3)) + 5)
    ninf = 1 - (1 / (1 + np.exp((vhalf - v)/K)))

    xgate = moose.element( h.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun
    xgate.tableB = 1.0/taun
    return h

def CaT_Chan(name):
    CaT = moose.HHChannel( '/library/' + name )
    CaT.Ek = ECa
    CaT.Gbar = 300.0*SOMA_A
    CaT.Gk = 0.0
    CaT.Xpower = 2.0
    CaT.Ypower = 1.0
    CaT.Zpower = 1.0
    CaT.useConcentration = 1
    CaT.instant = 4

    tBase = 23.5
    celsius = 22
    gcatbar = 0
    ki = 0.001
    cai = 5e-5
    cao = 2
    tfa = 1
    tfi = 0.68

    v = np.arange(Vmin,Vmax+dV, dV)
    alph = 1.6e-1*np.exp(-(v+57e-3)/19e-3)
    beth = 1e3/(np.exp((-v+15e-3)/10e-3)+1.0)
    alpm = 0.1967e6*(-1.0*v+19.88e-3)/(np.exp((-1.0*v+19.88e-3)/10.0e-3)-1.0)
    betm = 0.046e3*np.exp(-v/22.73e-3)

    a = alpm
    taum = 1.0/(tfa*(a + betm))
    minf =  a/(a+betm)
    a = alph
    tauh = 1.0/(tfi*(a + beth))
    hinf = a/(a+beth)

    xgate = moose.element( CaT.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum
    xgate.tableB = 1.0/taum

    ygate = moose.element( CaT.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh
    ygate.tableB = 1.0/tauh

    zgate = moose.element( CaT.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    ca = np.arange(Camin,Camax+dCa, dCa)
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    #Does not activate calcium-dependent currents
    addmsg4 = moose.Mstring( CaT.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return CaT

def CaR_SChan(name):
    CaR_S = moose.HHChannel( '/library/' + name )
    CaR_S.Ek = ECa
    CaR_S.Gbar = 300.0*SOMA_A
    CaR_S.Gk = 0.0
    CaR_S.Xpower = 3.0
    CaR_S.Ypower = 1.0
    CaR_S.Zpower = 0.0

    celsius =34

    v = np.arange(Vmin,Vmax+dV, dV)
    infm = 1.0 / (1 + np.exp((v+60e-3)/(-3e-3)))
    taum = 100e-3+0*v
    infh = 1.0 / (1 + np.exp((v+62e-3)/(1e-3)))
    tauh = 5e-3+0*v

    xgate = moose.element( CaR_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = infm/taum
    xgate.tableB = 1.0/taum

    ygate = moose.element( CaR_S.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = infh/tauh
    ygate.tableB = 1.0/tauh

    addmsg1 = moose.Mstring( CaR_S.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return CaR_S

def CaR_DChan(name):
    pass
    return CaR_D

def CaL_SChan(name):
    CaL_S = moose.HHChannel( '/library/' + name )
    CaL_S.Ek = ECa
    CaL_S.Gbar = 300.0*SOMA_A
    CaL_S.Gk = 0.0
    CaL_S.Xpower = 1.0
    CaL_S.Ypower = 0.0
    CaL_S.Zpower = 1.0
    CaL_S.useConcentration = 1
    CaL_S.instant = 4

    celsius= 34
    ki=.001
    cai = 50e-6
    cao = 2
    tfa = 5

    v = np.arange(Vmin,Vmax+dV, dV)
    alpm = 0.055e6*(-27.01e-3 - v)/(np.exp((-27.01e-3-v)/3.8e-3) - 1)
    betm =0.94e3*np.exp((-63.01e-3-v)/17e-3)
    a = alpm
    taum = 1/(tfa*(a+betm))
    minf = a/(a+betm)

    xgate = moose.element( CaL_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum
    xgate.tableB = 1.0/taum

    zgate = moose.element( CaL_S.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    ca = np.arange(Camin,Camax+dCa, dCa)
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( CaL_S.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( CaL_S.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return CaL_S

#Modified to fit data from Hirschberg et. al., 1999.
def KSK_Chan(name):
    KSK = moose.HHChannel( '/library/' + name )
    KSK.Ek = EK
    KSK.Gbar = 300.0*SOMA_A
    KSK.Gk = 0.0
    KSK.Xpower = 0.0
    KSK.Ypower = 0.0
    KSK.Zpower = 1.0
    KSK.useConcentration = 1

    celsius = 36
    cac     = 0.56e-3
    taumin  = 6.3e-3

    ca = np.arange(Camin,Camax+dCa, dCa)
    car = (ca/cac)**4.6
    m_inf = car / ( 1 + car )
    tau_m =  taumin + 0*ca

    zgate = moose.element( KSK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = m_inf/tau_m
    zgate.tableB = 1.0/tau_m

    addmsg3 = moose.Mstring( KSK.path + '/addmsg3' )
    addmsg3.value = '../Ca_conc    concOut    . concen'
    return KSK

def KBK_Chan(name):
    KBK = moose.HHChannel2D( '/library/' + name )
    KBK.Ek = EK
    KBK.Gbar = 300.0*SOMA_A
    KBK.Gk = 0.0
    KBK.Xpower = 0.0
    KBK.Ypower = 0.0
    KBK.Zpower = 1.0
    KBK.Zindex = 'VOLT_C1_INDEX'

    celsius = 20
    gkbar = 0.01e4
    d1 =1
    d2 = 1.5
    k1 = 0.18
    k2 = 0.011
    bbar = 0.28e3
    abar = 0.48e3
    F_KC = F/1000.0

    def alp(v,c):
        return c*abar/(c + k1*np.exp(-2*d1*F*v/R/(273.15 + celsius)))
    def bet(v,c):
        return bbar/(1 + c/(k2*np.exp(-2*d2*F*v/R/(273.15 + celsius))))

    zgate = moose.element( KBK.path + '/gateZ' )
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
    tblA = np.zeros([Vdivs+1,Cadivs+2])
    tblB = np.zeros([Vdivs+1,Cadivs+2])
    for i in np.arange(Vdivs+1):
        tblA[i] = [alp(x,y) for y in np.arange(Camin,Camax+dCa,dCa)]
        tblB[i] = np.add(tblA[i],[bet(x, y) for y in np.arange(Camin,Camax+dCa,dCa)])
        x = x + dV
        print(i, end='\r')

    zgate.tableA = tblA*1e3
    zgate.tableB = tblB*1e3

    addmsg4 = moose.Mstring( KBK.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return KBK

def Ca_Conc ( name ):
    Ca = moose.CaConc( '/library/' + name )
    depth	= 0.5e-6
    taur	= 1400e-3
    cainf	= 100e-6
    # B = 28789637.7

    Ca.tau = taur
    Ca.Ca_base = cainf
    Ca.diameter = 500e-6
    Ca.length = 1000e-6
    Ca.thick = 177.9e-6
    return Ca
