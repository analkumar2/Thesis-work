#exec(open('Kinetics_comparisons/ChannelProtos_Sri2015_base.py').read())
#This is the original without any tuning by Srikanth's original base code. No extra temp compensation etc.
#ECa and calcium decay dynamics not same as that in the paper.
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

def Na_Chan(name):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = ENa
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 3.0
    Na.Ypower = 1.0
    Na.Zpower = 1

    vhalfm  =  -30
    qa   = 7.2
    Ra   = 0.4
    Rb   = 0.124
    vhalfh  = -45
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
    stfactor = 1
    htfactor = 1
    mtfactor = 1

    def trap0(v,th,a,q):
        if np.abs(v*1e3-th) > 1e-6:
            return a * (v*1e3 - th) / (1 - np.exp(-(v*1e3 - th)/q))
        else:
            return a * q

    qt=q10**((celsius-24)/10)
    a = np.array([trap0(vm,vhalfm,Ra,qa) for vm in v])
    b = np.array([trap0(-vm,-vhalfm,Rb,qa) for vm in v])
    mtau = mtfactor*1/(a+b)/qt
    mtau[mtau<mmin] = mmin
    minf = a/(a+b)

    a = np.array([trap0(vm,vhalfh,Rd,qd) for vm in v])
    b = np.array([trap0(-vm,-vhalfh,Rg,qg) for vm in v])
    htau =  htfactor*1/(a+b)/qt
    htau[htau<hmin] = hmin
    hinf = 1/(1+np.exp((v*1e3-thinf)/qinf))

    c = 1/(1+np.exp((v*1e3-vvh)/vvs))
    sinf = c+ar2*(1-c)
    alps = np.exp(1.e-3*zetas*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    bets = np.exp(1.e-3*zetas*gms*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    taus = stfactor*bets/(a0s*(1+alps))
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

def K_DR_Chan(name):
    K_DR = moose.HHChannel( '/library/' + name )
    K_DR.Ek = EK
    K_DR.Gbar = 300.0*SOMA_A
    K_DR.Gk = 0.0
    K_DR.Xpower = 1.0
    K_DR.Ypower = 0.0
    K_DR.Zpower = 0.0

    vhalfn=13
    a0n=0.02
    zetan=-3
    gmn=0.7
    nmax=2
    q10=1
    ntfactor = 1

    qt=q10**((celsius-24)/10)
    alpn = np.exp(1.e-3*zetan*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(1.e-3*zetan*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1+alpn)
    taun = ntfactor*betn/(qt*a0n*(1+alpn))
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
    ntfactor = 1
    ltfactor = 1

    qt = q10**((celsius-24)/10)
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    alpn = np.exp(1.e-3*zeta*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    betn = np.exp(1.e-3*zeta*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1 + alpn)
    taun = ntfactor*betn/(qt*a0n*(1+alpn))
    taun[taun<nmin] = nmin

    alpl = np.exp(1.e-3*zetal*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    betl = np.exp(1.e-3*zetal*gml*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    linf = 1/(1+ alpl)
    taul = ltfactor*0.26*(v*1e3+50)/qtl
    taul[v<lmin/qtl] =lmin/qtl

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

def h_Chan(name):
    h = moose.HHChannel( '/library/' + name )
    h.Ek = Eh
    h.Gbar = 300.0*SOMA_A
    h.Gk = 0.0
    h.Xpower = 1.0
    h.Ypower = 0.0
    h.Zpower = 0.0

    vhalfl=-81
    kl=-8
    vhalft=-75
    a0t=0.011
    zetal=4
    zetat=2.2
    gmt=.4
    q10=4.5
    qtl=1
    ltfactor=1

    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))
    qt=q10**((celsius-33)/10)
    a = alpt
    linf = 1/(1 + np.exp(-(v*1e3-vhalfl)/kl))
    taul = ltfactor*bett/(qtl*qt*a0t*(1+a))

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
    Ca_T.Zpower = 0.0

    gCa_Tbar=.003
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
    mtfactor=1
    htfactor=1

    alph = np.exp(0.0378*zetah*(v*1e3-vhalfh))
    beth = np.exp(0.0378*zetah*gmh*(v*1e3-vhalfh))
    alpmt = np.exp(0.0378*zetam*(v*1e3-vhalfm))
    betmt = np.exp(0.0378*zetam*gmm*(v*1e3-vhalfm))
    qt=q10**((celsius-25)/10)

    minf = (1/(1 + np.exp(-(v*1e3-vhm)/10)))
    mtau = mtfactor*betmt/(qt*a0m*(1+alpmt))
    mtau[mtau<mmin] = mmin

    hinf = (1/(1 + np.exp((v*1e3-vhh)/10)))
    htau =  htfactor*beth/(a0h*(1+alph))
    htau[htau<hmin] = hmin

    xgate = moose.element( Ca_T.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/mtau*1e3
    xgate.tableB = 1.0/mtau*1e3

    ygate = moose.element( Ca_T.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/htau*1e3
    ygate.tableB = 1.0/htau*1e3

    addmsg1 = moose.Mstring( Ca_T.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return Ca_T

def Ca_R_Chan(name):
    Ca_R = moose.HHChannel( '/library/' + name )
    Ca_R.Ek = ECa
    Ca_R.Gbar = 300.0*SOMA_A
    Ca_R.Gk = 0.0
    Ca_R.Xpower = 3.0
    Ca_R.Ypower = 1.0
    Ca_R.Zpower = 0.0

    q10  = 3
    taum_exp = 0.92
    z = 2
    vhalfm = 3
    vhalfh = -39
    mtfactor = 1
    htfactor = 1

    qt=q10**(-(celsius-22)/10)
    minf = 1/(1+np.exp(-(v*1e3- vhalfm)/8.3))
    taum = (mtfactor*qt*taum_exp + 0*v)
    hinf = 1/(1+np.exp( (v*1e3-vhalfh)/9.2))
    tauh = (htfactor*qt*53 + 0*v)


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
    vhalfa = -27.01
    mtfactor = 1

    alpm = 0.055*(vhalfa - v*1e3)/(np.exp((vhalfa-v*1e3)/3.8) - 1)
    betm = 0.94*np.exp((-63.01-v*1e3)/17)
    a = alpm
    taum = mtfactor*1/(5*(a+betm))
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

def Ca_N_Chan(name):
    Ca_N = moose.HHChannel( '/library/' + name )
    Ca_N.Ek = ECa
    Ca_N.Gbar = 300.0*SOMA_A
    Ca_N.Gk = 0.0
    Ca_N.Xpower = 2.0
    Ca_N.Ypower = 1.0
    Ca_N.Zpower = 1.0
    Ca_N.useConcentration = 1
    Ca_N.instant = 4

    ki=.001
    cai = 50e-6
    cao = 10
    mtfactor = 1
    htfactor = 1
    vhalfm = 19.88
    vhalfh = 39.0

    alph = 1.6e-4*np.exp(-v*1e3/48.4)
    beth = 1/(np.exp((-v*1e3+vhalfh)/10.)+1.)
    alpm = 0.1967*(-1.0*v*1e3+vhalfm)/(np.exp((-1.0*v*1e3+vhalfm)/10.0)-1.0)
    betm = 0.046*np.exp(-v*1e3/20.73)
    a = alpm
    taum = mtfactor*1/(a + betm)
    minf = a*mtfactor*1/(a + betm)
    a = alph
    tauh = htfactor*1/(a + beth)
    hinf = a*htfactor*1/(a + beth)

    xgate = moose.element( Ca_N.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum*1e3
    xgate.tableB = 1.0/taum*1e3

    ygate = moose.element( Ca_N.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh*1e3
    ygate.tableB = 1.0/tauh*1e3

    zgate = moose.element( Ca_N.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( Ca_N.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( Ca_N.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return Ca_N

def K_SK_Chan(name):
    K_SK = moose.HHChannel( '/library/' + name )
    K_SK.Ek = EK
    K_SK.Gbar = 300.0*SOMA_A
    K_SK.Gk = 0.0
    K_SK.Xpower = 0.0
    K_SK.Ypower = 0.0
    K_SK.Zpower = 3.0
    K_SK.useConcentration = 1

    n=4
    cai=50.e-6
    a0=1.3e13
    b0=.5e-2
    tfactor = 1
    cahalf=140.04e-6
    k=0.10857

    alp = a0*ca**n
    a = alp
    tau = tfactor*1/(a + b0)
    inf=1/(1+np.exp((np.log(cahalf/ca))/k))

    zgate = moose.element( K_SK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = inf/tau*1e3
    zgate.tableB = 1.0/tau*1e3

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

    cai = 5.e-5
    d1 = .84
    d2 = 1.
    k1 = .48e-3
    k2 = .13e-6
    abar = .28
    bbar = .48
    st=1
    tfactor = 1
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
        tblB[i] = np.add(tblA[i],[bet(x, y) for y in ca])/tfactor
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
    cainf    = 50e-6
    # B = 28789637.7

    Ca.tau = taur/7
    Ca.Ca_base = cainf
    Ca.diameter = 500e-6
    Ca.length = 1000e-6
    Ca.thick = 177.9e-6
    return Ca
