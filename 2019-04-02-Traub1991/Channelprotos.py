#Author: Anal Kumar
#Defines the dynamics of NaChan, CaChan, KdrChan, KaChan, KahpChan, KcCHan, Lchan
#LChan

import moose
import rdesigneur as rd
import numpy as np

SOMA_A = 3.32e-9
EREST = -0.06 #EREST here is not the resting membrane potential. It's actually the Em of the membrane.
# EREST = 0
Ek_Na = 0.115 + EREST
Ek_K = -0.015 + EREST
Ek_Ca = 0.140 + EREST
vmin = EREST-0.050
vmax = EREST+0.150
vdivs = 3000.0
Camin = 0.0
Camax = 600.0 #changed from 0.02
# Camax = 500/25000
Cadivs = 3000.0
F = 96485.3329
Ca_scale = 1.0

def NaChan( name ):
    Na = moose.HHChannel( '/library/' + name )
    Na.Ek = Ek_Na
    Na.Gbar = 300.0*SOMA_A
    Na.Gk = 0.0
    Na.Xpower = 2.0
    Na.Ypower = 1.0
    Na.Zpower = 0.0
    
    xgate = moose.element( Na.path + '/gateX' )
    xgate.alpha = np.array([0.32*13.1e3 + 0.32e6*EREST, -0.320e6, -1.0, -13.1e-3-EREST, -4e-3])
    xgate.beta = np.array([-40.1*0.28e3-0.28e6*EREST, 0.28e6, -1.0, -40.1e-3-EREST, 5e-3])
    xgate.min = vmin
    xgate.max = vmax
    xgate.divs = vdivs
    
    ygate = moose.element( Na.path + '/gateY' )
    ygate.alpha = np.array([0.128e3, 0.0, 0.0, -17e-3-EREST, 18e-3])
    ygate.beta = np.array([4e3, 0.0, 1.0, -40e-3-EREST, -5e-3])
    ygate.min = vmin
    ygate.max = vmax
    ygate.divs = vdivs
    return Na

def CaChan ( name ):
    Ca = moose.HHChannel( '/library/' + name )
    Ca.Ek = Ek_Ca
    Ca.Gbar = 40.0 * SOMA_A
    Ca.Gk = 0.0
    Ca.Xpower = 2.0
    Ca.Ypower = 1.0
    Ca.Zpower = 0.0
    
    xgate = moose.element( Ca.path + '/gateX' )
    xgate.alpha = np.array([1.6e3, 0.0, 1.0, -65e-3-EREST, -13.89e-3])
    xgate.beta = np.array([-51.1*0.02e3-0.02e6*EREST, 0.02e6, -1.0, -51.1e-3-EREST, 5e-3])
    xgate.min = vmin
    xgate.max = vmax
    xgate.divs = vdivs
    
    ygate = moose.element( Ca.path + '/gateY' )
    ygate.min = vmin
    ygate.max = vmax
    ygate.divs = vdivs
    yA = np.zeros( (ygate.divs + 1), dtype=float)
    yB = np.zeros( (ygate.divs + 1), dtype=float)
    dx = (ygate.max - ygate.min)/ygate.divs
    
    x = ygate.min
    for i in range( ygate.divs + 1 ):
        if x<EREST:
            yA[i] = 5.0 #changed from yA[i] = 0
            # yA[i] = 0
        else:
            yA[i] = np.exp(-50.0*x + 50.0*EREST)/200e-3
        yB[i] = 5.0
        x = x+dx
    ygate.tableA = yA
    ygate.tableB = yB
    
    addmsg1 = moose.Mstring( Ca.path + '/addmsg1' )
    addmsg1.value = '.    IkOut    ../Ca_conc    current'
    return Ca
    
def KdrChan( name ):
    Kdr = moose.HHChannel( '/library/' + name )
    Kdr.Ek = Ek_K
    Kdr.Gbar = 150.0*SOMA_A
    Kdr.Gk = 0.0
    Kdr.Xpower = 1.0
    Kdr.Ypower = 0.0
    Kdr.Zpower = 0.0
    
    xgate = moose.element( Kdr.path + '/gateX' )
    xgate.alpha = np.array([0.016*35.1e3+0.016e6*EREST, -0.016e6, -1.0, -35.1e-3-EREST, -5e-3])
    xgate.beta = np.array([0.25e3, 0.0, 0.0, -20e-3-EREST, 40e-3])
    xgate.min = vmin
    xgate.max = vmax
    xgate.divs = vdivs  
    return Kdr

def KaChan ( name ):
    Ka = moose.HHChannel( '/library/' + name )
    Ka.Ek = Ek_K
    Ka.Gbar = 50.0*SOMA_A
    Ka.Gk = 0.0
    Ka.Xpower = 1.0
    Ka.Ypower = 1.0
    Ka.Zpower = 0.0
    
    xgate = moose.element( Ka.path + '/gateX' )
    xgate.alpha = np.array([0.02*13.1e3+0.02e6*EREST, -0.02e6, -1.0, -13.1e-3-EREST, -10e-3])
    xgate.beta = np.array([-0.0175*40.1e3-0.0175e6*EREST, 0.0175e6, -1.0, -40.1e-3-EREST, 10e-3 ])
    xgate.min = vmin
    xgate.max = vmax
    xgate.divs = vdivs
    
    ygate = moose.element( Ka.path + '/gateY' )
    ygate.alpha = np.array([0.0016e3, 0.0, 0.0, 13e-3-EREST, 18e-3])
    ygate.beta = np.array([0.05e3, 0.0, 1.0, -10.1e-3-EREST, -5e-3])
    ygate.min = vmin
    ygate.max = vmax
    ygate.divs = vdivs
    return Ka

def KahpChan( name ):
    Kahp = moose.HHChannel( '/library/' + name )
    Kahp.Ek = Ek_K
    Kahp.Gbar = 8.0*SOMA_A
    Kahp.Gk = 0.0
    Kahp.Xpower = 0.0
    Kahp.Ypower = 0.0
    Kahp.Zpower = 1.0
    Kahp.useConcentration = 1 #changed from Kahp.useConcentration = 0
    
    zgate = moose.element( Kahp.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zA = np.zeros( (zgate.divs + 1), dtype=float)
    zB = np.zeros( (zgate.divs + 1), dtype=float)
    dx = (zgate.max - zgate.min)/zgate.divs
    x = zgate.min
    for i in range( zgate.divs + 1 ):
        # zA[i] = min( 500.00 * x * Ca_scale, 10 ) #in paper, its zA[i] = min( 0.02 * x, 10 )
        zA[i] = min( 0.02 * x, 10 ) #changed from zA[i] = min( 250.00 * x, 10 )
        # zA[i] = min( 0.02 * x, 10.0 )
        zB[i] = 1.0 + zA[i] #changed from zB[i] = 1.0
        x = x + dx
    
    zgate.tableA = zA
    zgate.tableB = zB
    
    addmsg1 = moose.Mstring( Kahp.path + '/addmsg1' )
    addmsg1.value = '../Ca_conc    concOut    . concen'
    return Kahp
    
def KahpChanKO( name ):
    Kahp = moose.HHChannel( '/library/' + name )
    Kahp.Ek = Ek_K
    Kahp.Gbar = 8.0*SOMA_A
    Kahp.Gk = 0.0
    Kahp.Xpower = 0.0
    Kahp.Ypower = 0.0
    Kahp.Zpower = 1.0
    Kahp.useConcentration = 1 #changed from Kahp.useConcentration = 0
    
    zgate = moose.element( Kahp.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zA = np.zeros( (zgate.divs + 1), dtype=float)
    zB = np.zeros( (zgate.divs + 1), dtype=float)
    dx = (zgate.max - zgate.min)/zgate.divs
    x = zgate.min
    for i in range( zgate.divs + 1 ):
        # zA[i] = min( 500.00 * x * Ca_scale, 10 ) #in paper, its zA[i] = min( 0.02 * x, 10 )
        zA[i] = min( 0.02 * x/2, 10 ) #changed from zA[i] = min( 250.00 * x, 10 ) #2 times decrease in sensitivity
        # zA[i] = min( 0.02 * x, 10.0 )
        zB[i] = 1.0 + min( 0.02 * x, 10 ) #changed from zB[i] = 1.0
        x = x + dx
    
    zgate.tableA = zA
    zgate.tableB = zB
    
    addmsg1 = moose.Mstring( Kahp.path + '/addmsg1' )
    addmsg1.value = '../Ca_conc    concOut    . concen'
    return Kahp

def KcChan ( name ):
    Kc = moose.HHChannel( '/library/' + name )
    Kc.Ek = Ek_K
    Kc.Gbar = 100.0*SOMA_A
    Kc.Gk = 0.0
    Kc.Xpower = 1.0
    Kc.Ypower = 0.0
    Kc.Zpower = 1.0
    Kc.instant = 4
    Kc.useConcentration = 1
    
    xgate = moose.element( Kc.path + '/gateX' )
    xgate.min = vmin
    xgate.max = vmax
    xgate.divs = vdivs
    xA = np.zeros( (xgate.divs + 1), dtype=float)
    xB = np.zeros( (xgate.divs + 1), dtype=float)
    dx = (xgate.max - xgate.min)/xgate.divs
    x = xgate.min
    for i in range( xgate.divs + 1 ):
        if (x < EREST + 0.05):
            xA[i] = np.exp((x - 10e-3 - EREST)/11e-3 - (x - 6.5e-3 - EREST)/27e-3)/18.975e-3
            xB[i] = 2000.0*np.exp((6.5e-3 + EREST - x)/27e-3)
        else:
            xA[i] = 2000.0*np.exp((6.5e-3 + EREST - x)/27e-3)
            xB[i] = 2000.0*np.exp((6.5e-3 + EREST - x)/27e-3)
        x = x + dx
    xgate.tableA = xA
    xgate.tableB = xB
    
    zgate = moose.element( Kc.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zA = np.zeros( (zgate.divs + 1), dtype=float)
    zB = np.zeros( (zgate.divs + 1), dtype=float)
    dx = ( zgate.max -  zgate.min)/ zgate.divs
    x = zgate.min
    for i in range( zgate.divs + 1 ):
        # zA[i] = min( 1.0, x*Ca_scale/0.02 ) #in the paper, its zA[i] = min( 1.0, x/250 )
        zA[i] = min( 1.0, x/250.0 )
        # zA[i] = min( 1.0, x/250 )
        zB[i] = 1.0
        x += dx
    zgate.tableA = zA
    zgate.tableB = zB
    
    addmsg1 = moose.Mstring( Kc.path + '/addmsg1' )
    addmsg1.value = '../Ca_conc    concOut    . concen'  
    return Kc

#Leak channel # Not used in the code    
def LChan ( name ):
    L = moose.Leakage( '/library/' + name )
    L.Ek = EREST
    L.Gbar = 1*SOMA_A
    L.Gk = 1*SOMA_A
    return L
    
def CaConc ( name ):
    conc = moose.CaConc( '/library/' + name )
    conc.tau = 0.013333
    conc.Ca_base = 0.00000
    conc.thick = 1.79483818751e-10 #Changed from conc.B  = 17.402e12
    #conc length set up same as the length of the compartment
    return conc
    
def QUISrec ( name ) :
    
    return quis
    
# Simplest HHchannel #not used in the code
def simpChan(name):
    simp = moose.HHChannel('/library/' + name )
    simp.Ek = Ek_K
    simp.Gbar = 150*SOMA_A
    simp.Gk = 0
    simp.Xpower = 0
    simp.Ypower = 0
    simp.Zpower = 1
    simp.instant = 4
    simp.useConcentration = 1
    
    zgate = moose.element( simp.path + '/gateZ' )
    zgate.min = Camin
    zgate.max = Camax
    zgate.divs = Cadivs
    zA = np.zeros( (zgate.divs + 1), dtype=float)
    zB = np.zeros( (zgate.divs + 1), dtype=float)
    dx = ( zgate.max -  zgate.min)/ zgate.divs
    x = zgate.min
    for i in range( zgate.divs + 1 ):
        zA[i] = min( 1.0, x/zgate.max ) #in the paper, its zA[i] = min( 1.0, x/250 )
        # zA[i] = min( 1.0, x/250 )
        zB[i] = 1.0
        x += dx
    zgate.tableA = zA
    zgate.tableB = zB
    
    addmsg1 = moose.Mstring( simp.path + '/addmsg1' )
    addmsg1.value = '../Ca_conc    concOut    . concen'
    return simp
    
