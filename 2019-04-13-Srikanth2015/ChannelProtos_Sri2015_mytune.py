#exec(open('Srikanth2015/ChannelProtos_Sri2015_mytune.py').read())
#This is after tuning by Srikanth's original base code. During mytuning, this file is not altered
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
Camin = 1e-12
Camax = 3e-3
Cadivs = 3000
dCa = (Camax-Camin)/Cadivs
ca = np.arange(Camin,Camax+dCa, dCa)

def Na_SChan(name):
    Na_S = moose.HHChannel( '/library/' + name )
    Na_S.Ek = ENa
    Na_S.Gbar = 300.0*SOMA_A
    Na_S.Gk = 0.0
    Na_S.Xpower = 3.0
    Na_S.Ypower = 1.0
    Na_S.Zpower = 1

    vhalfm  =  -29.758757
    qa   = 7.2
    Ra   = 0.4
    Rb   = 0.124
    vhalfh  = -44.151697
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
    vhalfs=-59.510922
    a0s=0.0003
    zetas=12
    gms=0.2
    smax=10
    vvh=-58
    vvs=2
    ar2=1
    stfactor = 0.986498
    htfactor = 0.816827
    mtfactor = 1.084830

    def trap0(v,th,a,q):
        if np.abs(v-th) > 1e-6:
            return a * (v*1e3 - th) / (1 - np.exp(-(v*1e3 - th)/q))
        else:
            return a * q

    qt=q10**((celsius-24)/10)
    a = np.array([trap0(vm,vhalfm,Ra,qa) for vm in v])
    b = np.array([trap0(-vm,-vhalfm,Rb,qa) for vm in v])
    mtau = mtfactor*1/(a+b)/qt*1e-3
    mtau[mtau<mmin*1e-3] = mmin*1e-3
    minf = a/(a+b)

    a = np.array([trap0(vm,vhalfh,Rd,qd) for vm in v])
    b = np.array([trap0(-vm,-vhalfh,Rg,qg) for vm in v])
    htau =  htfactor*1/(a+b)/qt *1e-3
    htau[htau<hmin*1e-3] = hmin*1e-3
    hinf = 1/(1+np.exp((v*1e3-thinf)/qinf))

    c = 1/(1+np.exp((v*1e3-vvh)/vvs))
    sinf = c+ar2*(1-c)
    alps = np.exp(1.e-3*zetas*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    bets = np.exp(1.e-3*zetas*gms*(v*1e3-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
    taus = stfactor*bets/(a0s*(1+alps))*1e-3
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

def KDR_SChan(name):
    KDR_S = moose.HHChannel( '/library/' + name )
    KDR_S.Ek = EK
    KDR_S.Gbar = 300.0*SOMA_A
    KDR_S.Gk = 0.0
    KDR_S.Xpower = 1.0
    KDR_S.Ypower = 0.0
    KDR_S.Zpower = 0.0

    vhalfn=10.443284
    a0n=0.02
    zetan=-3
    gmn=0.7
    nmax=2
    q10=1
    ntfactor = 0.626055

    qt=q10**((celsius-24)/10)
    alpn = np.exp(1.e-3*zetan*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    betn = np.exp(1.e-3*zetan*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1+alpn)
    taun = ntfactor*betn/(qt*a0n*(1+alpn))*1e-3
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

    vhalfn=8.294468
    vhalfl=-58.836023
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
    ntfactor = 1.790612
    ltfactor = 10.678116

    qt = q10**((celsius-24)/10)
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    alpn = np.exp(1.e-3*zeta*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    zeta = zetan+pw/(1+np.exp((v*1e3-tq)/qq))
    betn = np.exp(1.e-3*zeta*gmn*(v*1e3-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
    ninf = 1/(1 + alpn)
    taun = ntfactor*betn/(qt*a0n*(1+alpn))*1e-3
    taun[taun<nmin*1e-3] = nmin*1e-3

    alpl = np.exp(1.e-3*zetal*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    betl = np.exp(1.e-3*zetal*gml*(v*1e3-vhalfl)*9.648e4/(8.315*(273.16+celsius)))
    linf = 1/(1+ alpl)
    taul = ltfactor*0.26*(v*1e3+50)/qtl*1e-3
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

def KM_Chan(name):
    KM = moose.HHChannel( '/library/' + name )
    KM.Ek = EK
    KM.Gbar = 300.0*SOMA_A
    KM.Gk = 0.0
    KM.Xpower = 1.0
    KM.Ypower = 0.0
    KM.Zpower = 0.0

    vhalfl=-44.607246
    kl=-10
    vhalft=-42
    a0t=0.009
    zetat=7
    gmt=.4
    q10=5
    b0=60
    st=1
    tfactor = 1.442796

    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))
    qt=q10**((celsius-35)/10)
    inf = (1/(1 + np.exp((v*1e3-vhalfl)/kl)))
    a = alpt
    tau = tfactor*(b0 + bett)/(a0t*(1+a))*1e-3

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

    vhalfl=-77.553235
    kl=-8
    vhalft=-75
    a0t=0.011
    zetal=4
    zetat=2.2
    gmt=.4
    q10=4.5
    qtl=1
    ltfactor=1.133596

    alpt = np.exp(0.0378*zetat*(v*1e3-vhalft))
    bett = np.exp(0.0378*zetat*gmt*(v*1e3-vhalft))
    qt=q10**((celsius-33)/10)
    a = alpt
    linf = 1/(1 + np.exp(-(v*1e3-vhalfl)/kl))
    taul = ltfactor*bett/(qtl*qt*a0t*(1+a))*1e-3

    xgate = moose.element( h.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = linf/taul
    xgate.tableB = 1.0/taul
    return h

def CaT_Chan(name):
    CaT = moose.HHChannel( '/library/' + name )
    CaT.Ek = ECa
    CaT.Gbar = 300.0*SOMA_A
    CaT.Gk = 0.0
    CaT.Xpower = 2.0
    CaT.Ypower = 1.0
    CaT.Zpower = 0.0

    celsius = 25
    gcatbar=.003
    cao = 2
    q10 = 5
    mmin=0.2
    hmin=10
    a0h =0.015
    zetah = 3.5
    vhalfh = -70.018017
    gmh=0.6
    a0m =0.04
    zetam = 2
    vhalfm = -17.868123
    gmm=0.61
    vhm=-60
    vhh=-85
    mtfactor=1.552435
    htfactor=0.963171

    alph = np.exp(0.0378*zetah*(v*1e3-vhalfh))
    beth = np.exp(0.0378*zetah*gmh*(v*1e3-vhalfh))
    alpmt = np.exp(0.0378*zetam*(v*1e3-vhalfm))
    betmt = np.exp(0.0378*zetam*gmm*(v*1e3-vhalfm))
    qt=q10**((celsius-25)/10)

    minf = (1/(1 + np.exp(-(v*1e3-vhm)/10)))
    mtau = mtfactor*betmt/(qt*a0m*(1+alpmt))*1e-3
    mtau[mtau<mmin*1e-3] = mmin*1e-3

    hinf = (1/(1 + np.exp((v*1e3-vhh)/10)))
    htau =  htfactor*beth/(a0h*(1+alph))*1e-3
    htau[htau<hmin*1e-3] = hmin*1e-3

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

    addmsg1 = moose.Mstring( CaT.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return CaT

def CaR_SChan(name):
    CaR_S = moose.HHChannel( '/library/' + name )
    CaR_S.Ek = ECa
    CaR_S.Gbar = 300.0*SOMA_A
    CaR_S.Gk = 0.0
    CaR_S.Xpower = 3.0
    CaR_S.Ypower = 1.0
    CaR_S.Zpower = 0.0

    q10  = 3
    taum_exp = 0.92
    z = 2
    vhalfm = 0.058617
    vhalfh = -42.085754
    mtfactor = 1.482045
    htfactor = 0.615933

    qt=q10**((celsius-22)/10)
    minf = 1/(1+np.exp(-(v*1e3- vhalfm)/8.3))
    taum = (mtfactor*qt*taum_exp + 0*v)*1e-3
    hinf = 1/(1+np.exp( (v*1e3-vhalfh)/9.2))
    tauh = (htfactor*qt*53 + 0*v)*1e-3


    xgate = moose.element( CaR_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum
    xgate.tableB = 1.0/taum

    ygate = moose.element( CaR_S.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh
    ygate.tableB = 1.0/tauh

    addmsg1 = moose.Mstring( CaR_S.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return CaR_S

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
    vhalfa = -29.743401
    mtfactor = 1.840658

    alpm = 0.055*(vhalfa - v*1e3)/(np.exp((vhalfa-v*1e3)/3.8) - 1)
    betm = 0.94*np.exp((-63.01-v*1e3)/17)
    a = alpm
    taum = mtfactor*1/(5*(a+betm))*1e-3
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

def CaN_SChan(name):
    CaN_S = moose.HHChannel( '/library/' + name )
    CaN_S.Ek = ECa
    CaN_S.Gbar = 300.0*SOMA_A
    CaN_S.Gk = 0.0
    CaN_S.Xpower = 2.0
    CaN_S.Ypower = 1.0
    CaN_S.Zpower = 1.0
    CaN_S.useConcentration = 1
    CaN_S.instant = 4

    celsius= 34
    ki=.001
    cai = 50e-6
    cao = 10
    mtfactor = 1.495491
    htfactor = 1.760253
    vhalfm = 23.916356
    vhalfh = 39.676957

    alph = 1.6e-4*np.exp(-v*1e3/48.4)
    beth = 1/(np.exp((-v*1e3+vhalfh)/10.)+1.)
    alpm = 0.1967*(-1.0*v*1e3+vhalfm)/(np.exp((-1.0*v*1e3+vhalfm)/10.0)-1.0)
    betm = 0.046*np.exp(-v*1e3/20.73)
    a = alpm
    taum = mtfactor*1/(a + betm)*1e-3
    minf = a*mtfactor*1/(a + betm)
    a = alph
    tauh = htfactor*1/(a + beth)*1e-3
    hinf = a*htfactor*1/(a + beth)

    xgate = moose.element( CaN_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = minf/taum
    xgate.tableB = 1.0/taum

    ygate = moose.element( CaN_S.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = hinf/tauh
    ygate.tableB = 1.0/tauh

    zgate = moose.element( CaN_S.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = ki/(ki+ca)
    zgate.tableB = ki/(ki+ca)

    addmsg2 = moose.Mstring( CaN_S.path + '/addmsg2' )
    addmsg2.value = '.    IkOut    ../Ca_conc    current'

    addmsg4 = moose.Mstring( CaN_S.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return CaN_S

def KSK_Chan(name):
    KSK = moose.HHChannel( '/library/' + name )
    KSK.Ek = EK
    KSK.Gbar = 300.0*SOMA_A
    KSK.Gk = 0.0
    KSK.Xpower = 0.0
    KSK.Ypower = 0.0
    KSK.Zpower = 3.0
    KSK.useConcentration = 1

    n=4
    cai=50.e-6
    a0=1.3e13
    b0=.5e-2
    tfactor = 0.522782
    cahalf=0.000161933
    k=0.10857

    alp = a0*ca**n
    a = alp
    tau = tfactor*1/(a + b0)*1e-3
    inf=1/(1+np.exp((np.log(cahalf/ca))/k))

    zgate = moose.element( KSK.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zgate.tableA = inf/tau
    zgate.tableB = 1.0/tau

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

    cai = 5.e-5
    d1 = .84
    d2 = 1.
    k1 = 0.000588
    k2 = 8.76481e-08
    abar = .28
    bbar = .48
    st=1
    tfactor = 1.279493
    F_KC = F/1000.0

    def exp1(k,d,v):
        return k*np.exp(-2*d*F_KC*v*1e3/R/(273.15 + celsius))
    def alp(v,c):
        return c*abar/(c + exp1(k1,d1,v))
    def bet(v,c):
        return bbar/(1 + c/exp1(k2,d2,v))

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
    tblA = np.zeros([len(v),len(ca)])
    tblB = np.zeros([len(v),len(ca)])
    for i in np.arange(len(v)):
        tblA[i] = np.array([alp(x,y) for y in ca])
        tblB[i] = np.add(tblA[i],[bet(x, y) for y in ca])/tfactor
        x = x + dV
        print(i, end='\r')

    zgate.tableA = tblA*1e3
    zgate.tableB = tblB*1e3

    addmsg4 = moose.Mstring( KBK.path + '/addmsg4' )
    addmsg4.value = '../Ca_conc    concOut    . concen'
    return KBK

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
