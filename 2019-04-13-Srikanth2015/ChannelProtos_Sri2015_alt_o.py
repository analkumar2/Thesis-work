#exec(open('Srikanth2015/ChannelProtos_Sri2015_alt_o.py').read())
#This is the original without any tuning by Srikanth.
#Incomplete. Only Na, KDR, KA, KM, h, CaT, CaR, CaL
#After everyhting is done, check calcium connections
#Author: Anal Kumar
#Defines the channel kinetics

import moose
import rdesigneur as rd
import numpy as np
import csv

SOMA_A = 3.14e-8
F = 96485.3329
R = 8.314
Temp = 307.15
celsius = 34
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
Camin = 0
Camax = 3e-3
Cadivs = 3000
dCa = (Camax-Camin)/Cadivs
ca = np.arange(Camin,Camax+dCa, dCa)

def Na_SChan(name): # Hoffman et al.(1997) & Migliore et al. (1999) # Currently from Migliore1999 ModelDB.
    Na_S = moose.HHChannel( '/library/' + name )
    Na_S.Ek = ENa
    Na_S.Gbar = 300.0*SOMA_A
    Na_S.Gk = 0.0
    Na_S.Xpower = 3.0
    Na_S.Ypower = 1.0
    Na_S.Zpower = 1

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
    ar2=1

    def trap0(v,th,a,q):
        if np.abs(v-th) > 1e-6:
            return a * (v*1e3 - th) / (1 - np.exp(-(v*1e3 - th)/q))
        else:
            return a * q

    qt=q10**((celsius-24)/10)
    a = np.array([trap0(vm,tha,Ra,qa) for vm in v])
    b = np.array([trap0(-vm,-tha,Rb,qa) for vm in v])
    mtau = 1/(a+b)/qt*1e-3
    mtau[mtau<mmin*1e-3] = mmin*1e-3
    minf = a/(a+b)

    a = np.array([trap0(vm,thi1,Rd,qd) for vm in v])
    b = np.array([trap0(-vm,-thi2,Rg,qg) for vm in v])
    htau =  1/(a+b)/qt *1e-3
    htau[htau<hmin*1e-3] = hmin*1e-3
    hinf = 1/(1+np.exp((v*1e3-thinf)/qinf))

    c = 1/(1+np.exp((v*1e3-vvh)/vvs))
    sinf = c+ar2*(1-c)
    alps = np.exp(1.e-3*zetas*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    bets = np.exp(1.e-3*zetas*gms*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    taus = bets/(a0s*(1+alps))*1e-3
    taus[taus<smax*1e-3] = smax*1e-3

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

def KDR_SChan(name): # Hoffman et al.(1997) & Migliore et al. (1999) # Currently from Migliore1999.
    KDR_S = moose.HHChannel( '/library/' + name )
    KDR_S.Ek = EK
    KDR_S.Gbar = 300.0*SOMA_A
    KDR_S.Gk = 0.0
    KDR_S.Xpower = 1.0
    KDR_S.Ypower = 0.0
    KDR_S.Zpower = 0.0

    vhalfn=13
    a0n=0.02
    zetan=-3
    gmn=0.7
    nmax=2
    q10=1

    qt=q10**((celsius-24)/10)
    alpn = np.exp(1.e-3*zetan*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(1.e-3*zetan*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1+alpn)
    taun = betn/(qt*a0n*(1+alpn))*1e-3
    taun[taun<nmax*1e-3]=nmax*1e-3

    xgate = moose.element( KDR_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun
    xgate.tableB = 1/taun
    return KDR_S

def KA_SChan(name): # Hoffman et al.(1997) & Migliore et al. (1999) # Currently from Migliore1999.
    KA_S = moose.HHChannel( '/library/' + name )
    KA_S.Ek = EK
    KA_S.Gbar = 300.0*SOMA_A
    KA_S.Gk = 0.0
    KA_S.Xpower = 1.0
    KA_S.Ypower = 1.0
    KA_S.Zpower = 0.0

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
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    betn = np.exp(1.e-3*zeta*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1 + alpn)
    taun = betn/(qt*a0n*(1+alpn))*1e-3
    taun[taun<nmin*1e-3] = nmin*1e-3

    alpl = np.exp(1.e-3*zetal*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    betl = np.exp(1.e-3*zetal*gml*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    linf = 1/(1+ alpl)
    taul = 0.26*(v*1e3+50)/qtl*1e-3
    taul[v<lmin/qtl*1e-3] =lmin/qtl*1e-3

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

def KM_Chan(name): # Migliore et al. (2006)
    KM = moose.HHChannel( '/library/' + name )
    KM.Ek = EK
    KM.Gbar = 300.0*SOMA_A
    KM.Gk = 0.0
    KM.Xpower = 1.0
    KM.Ypower = 0.0
    KM.Zpower = 0.0

    xinf = 1/(1 + np.exp(-(v*1e3 + 35)/5))
    xtau = 1000/(3.3*np.exp((v*1e3 + 35)/40) + np.exp(-(v*1e3 + 35)/20))*1e-3

    xgate = moose.element( KM.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = xinf/xtau
    xgate.tableB = 1.0/xtau
    return KM

def h_Chan(name): # Magee (1998) and Poolos et al. (2002) # Currently from Poolos2002 paper.
    h = moose.HHChannel( '/library/' + name )
    h.Ek = Eh
    h.Gbar = 300.0*SOMA_A
    h.Gk = 0.0
    h.Xpower = 1.0
    h.Ypower = 0.0
    h.Zpower = 0.0

    vhalfl=-82
    vhalft=-75
    a0t=0.011
    zetal=4
    zetat=2.2
    gmt=.4
	q10=4.5
	qtl=1

    alpl = np.exp(0.0378*zetal*(v*1e3-vhalfl))
    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))
    qt=q10**((celsius-33)/10)
    linf = 1/(1+ alpl)
    taul = bett/(qtl*qt*a0t*(1+alpt))*1e-3

    xgate = moose.element( h.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = linf/taul
    xgate.tableB = 1.0/taul
    return h

def CaT_Chan(name): # Shah et al.(2011) #ModelDB
    CaT = moose.HHChannel( '/library/' + name )
    CaT.Ek = ECa
    CaT.Gbar = 300.0*SOMA_A
    CaT.Gk = 0.0
    CaT.Xpower = 2.0
    CaT.Ypower = 1.0
    CaT.Zpower = 0.0

	celsius = 25
	gcatbar=.003
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
	gmm=0.61
	vhm=-60
	vhh=-85

    alph = np.exp(0.0378*zetah*(v*1e3-vhalfh))
    beth = np.exp(0.0378*zetah*gmh*(v*1e3-vhalfh))
    alpmt = np.exp(0.0378*zetam*(v*1e3-vhalfm))
    betmt = np.exp(0.0378*zetam*gmm*(v*1e3-vhalfm))
    qt=q10**((celsius-25)/10)

    minf = (1/(1 + np.exp(-(v*1e3-vhm)/10)))
    mtau = betmt/(qt*a0m*(1+alpmt))*1e-3
    mtau[mtau<mmin*1e-3] = mtau=mmin*1e-3

    hinf = (1/(1 + np.exp((v*1e3-vhh)/10)))
    htau = beth/(a0h*(1+alph))*1e-3
	htau[htau<hmin*1e-3] = htau=hmin*1e-3

    xgate = moose.element( CaT.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau
    xgate.tableB = 1.0/mtau

    ygate = moose.element( CaT.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/htau
    ygate.tableB = 1.0/htau
    return CaT

def CaR_SChan(name): # Magee and Johnston (1995) and Poirazi et al. (2003) #Currently from Poirazi2003
    CaR_S = moose.HHChannel( '/library/' + name )
    CaR_S.Ek = ECa
    CaR_S.Gbar = 300.0*SOMA_A
    CaR_S.Gk = 0.0
    CaR_S.Xpower = 3.0
    CaR_S.Ypower = 1.0
    CaR_S.Zpower = 0.0

    celsius = 34

    infm = 1 / (1 + np.exp((v*1e3+60)/(-3)))
    taum = (100 + 0*v)*1e-3
    infh = 1/ (1 + np.exp((v*1e3+62)/(1)))
    tauh = (5 + 0*v)*1e-3

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

def CaL_SChan(name): # Magee and Johnston (1995) and Poirazi et al. (2003)
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

    alpm = 0.055*(-27.01 - v*1e3)/(np.exp((-27.01-v*1e3)/3.8) - 1)
    betm =0.94*np.exp((-63.01-v*1e3)/17)
    a = alpm
    taum = 1/(tfa*(a+betm))*1e-3
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
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( CaL_S.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( CaL_S.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return CaL_S

def CaN_SChan(name): # Migliore et al. (1995)
    pass
    return CaL_D

def KBK_Chan(name): # Moczydlowski and Latorre (1983)
    KsAHP = moose.HHChannel( '/library/' + name )
    KsAHP.Ek = EK
    KsAHP.Gbar = 300.0*SOMA_A
    KsAHP.Gk = 0.0
    KsAHP.Xpower = 0.0
    KsAHP.Ypower = 0.0
    KsAHP.Zpower = 3.0
    KsAHP.useConcentration = 1

    celsius = 36
    gbar    = 0.01
    beta    = 0.03e3
    cac     = 0.025
    taumin  = 0.5e-3

    tadj = 3**((celsius-22.0)/10)
    ca = np.arange(Camin,Camax+dCa, dCa)
    car = (ca/cac)**2
    m_inf = car / ( 1 + car )
    tau_m =  1 / beta / (1 + car) / tadj
    tau_m[tau_m<taumin] = taumin

    zgate = moose.element( KsAHP.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = m_inf/tau_m
    zgate.tableB = 1.0/tau_m

    addmsg3 = moose.Mstring( KsAHP.path + '/addmsg3' )
    addmsg3.value = '../Ca_conc    concOut    . concen'
    return KsAHP

def KSK_Chan(name): # Migliore et al. (1995)
    KmAHP = moose.HHChannel2D( '/library/' + name )
    KmAHP.Ek = EK
    KmAHP.Gbar = 300.0*SOMA_A
    KmAHP.Gk = 0.0
    KmAHP.Xpower = 0.0
    KmAHP.Ypower = 0.0
    KmAHP.Zpower = 1.0
    KmAHP.Zindex = 'VOLT_C1_INDEX'

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

    zgate = moose.element( KmAHP.path + '/gateZ' )
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

    zgate.tableA = tblA
    zgate.tableB = tblB

    addmsg4 = moose.Mstring( KmAHP.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return KmAHP

def Ca_Conc ( name ):
    Ca = moose.CaConc( '/library/' + name )
    depth	= 0.1e-6
    taur	= 200e-3
    cainf	= 100e-6
    # B = 28789637.7

    Ca.tau = 7*taur
    Ca.Ca_base = cainf
    Ca.diameter = 500e-6
    Ca.length = 1000e-6
    Ca.thick = 177.9e-6
    return Ca
