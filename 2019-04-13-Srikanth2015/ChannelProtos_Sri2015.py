#'Srikanth2015/ChannelProtos_Srikanth2015.py'
#Hoffman1997, Only Na, KDR, and KA done.
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

def Na_SChan(name): # Hoffman et al.(1997) & Migliore et al. (1999) # Currently from Hoof1997 as both are different. Its probably Migliore.
    #gbarNa = 4.2 #some units
    Na_S = moose.HHChannel( '/library/' + name )
    Na_S.Ek = ENa
    Na_S.Gbar = 300.0*SOMA_A
    Na_S.Gk = 0.0
    Na_S.Xpower = 3.0
    Na_S.Ypower = 1.0
    Na_S.Zpower = 0

    alpham = 0.182*(v*1e3+25)/(1-np.exp(-(32.5+v*1e3)/4.5))
    betam = 0.124*(-v*1e3-32.5)/(1-np.exp((v*1e3+32.5)/4.5))
    taum = 0.8/(aplham+betam)*1e-3
    infm = alpham/(alpham+betam)

    alphah = 0.08*(v*1e3+40)/(1-np.exp(-(v*1e3+40)/3))
    betah = 0.0005*(−v*1e3 −10)/(1−np.exp((v*1e3 + 10)/5))
    infh = 1/(1+np.exp((v*1e3+58)/5))
    tauh = 1/(alphah+betah)*1e-3

    xgate = moose.element( Na_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = infm/taum
    xgate.tableB = 1.0/taum

    ygate = moose.element( Na_S.path + '/gateY' )
    ygate.min = Vmin
    ygate.max = Vmax
    ygate.divs = Vdivs
    ygate.tableA = infh/tauh
    ygate.tableB = 1.0/tauh

    return Na_S

def KDR_SChan(name): # Hoffman et al.(1997) & Migliore et al. (1999) # Currently from Hoof1997 as both are different. Its probably Migliore.
    KDR_S = moose.HHChannel( '/library/' + name )
    KDR_S.Ek = -0.077
    KDR_S.Gbar = 300.0*SOMA_A
    KDR_S.Gk = 0.0
    KDR_S.Xpower = 4.0
    KDR_S.Ypower = 0.0
    KDR_S.Zpower = 0.0

    alpham = −0.0035*(v*1e3 + 30)/(np.exp((v*1e3 + 30)/−13) − 1)
    betam = 0.0035*(v*1e3 + 30)/(exp((v*1e3 + 30)/13) − 1)
    ninf = alpham/(alpham+betam)
    taun = 1.8e-3 + 0*v

    xgate = moose.element( KDR_S.path + '/gateX' )
    xgate.min = Vmin
    xgate.max = Vmax
    xgate.divs = Vdivs
    xgate.tableA = ninf/taun
    xgate.tableB = 1/taun
    return KDR_S

def KA_SChan(name): # Hoffman et al.(1997) & Migliore et al. (1999) # Currently from Hoof1997 as both are different. Its probably Migliore.
    KA_S = moose.HHChannel( '/library/' + name )
    KA_S.Ek = EK
    KA_S.Gbar = 300.0*SOMA_A
    KA_S.Gk = 0.0
    KA_S.Xpower = 4.0
    KA_S.Ypower = 1.0
    KA_S.Zpower = 0.0

    alphan = −0.01*(v*1e3 + 21.3)/(exp((v*1e3 + 21.3)/−35)−1)
    betan = 0.01*(v*1e3 + 21.3)/(exp((v*1e3 + 21.3)/35)−1)
    ninf = alphan/(alphan+betan)
    taun = 0.2e-3 + 0*v

    alphal = −0.01*(v*1e3 + 58)/(exp((v*1e3 + 58)/−8.2) − 1)
    betal = 0.01*(v*1e3 + 58)/(exp((v*1e3 + 58)/−8.2)−1)
    linf = alphal/(alphal+betal)
    taul = (5 + 2.6*(v*1e3 + 20)/10)*1e-3
    taul[v<-20*1e-3] = 5*1e-3

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

def h_Chan(name): # Magee (1998) and Poolos et al. (2002)
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

def CaT_Chan(name): # Shah et al.(2011)
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
    cai = 5.e-5
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

def CaR_SChan(name): # Magee and Johnston (1995) and Poirazi et al. (2003)
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
