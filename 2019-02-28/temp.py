import numpy as np
import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

####################################################################

vmin = -0.100
vmax = 0.100
vdivs = 3000

V12mNaT = -0.025
kmNaT = 0.0072
dV = (vmax-vmin)/vdivs    
AlphamNaT = np.zeros((vdivs+1), dtype=float)
BetamNaT = np.zeros((vdivs+1), dtype=float)
TaumNaT = np.zeros((vdivs+1), dtype=float)
V = vmin
for i in range(vdivs+1):
if abs(V-V12mNaT)>1e-9:
AlphamNaT[i] = 0.4e3*(V-V12mNaT)/(1-np.exp((V12mNaT-V)/kmNaT))
BetamNaT[i] = -0.124e3*(V-V12mNaT)/(1-np.exp((V-V12mNaT)/kmNaT))
else:
AlphamNaT[i] = 0.4e3*kmNaT
BetamNaT[i] = 0.124e3*kmNaT
V = V+dV

TaumNaT = 3e-3/(AlphamNaT+BetamNaT)
TaumNaT[TaumNaT<0.02e-3] = 0.02e-3

##################################################################
try:
[moose.delete(x) for x in ['/model', '/library']]
except:
pass
SOMA_A = 3.32e-9
Temp = 307.15
dt = 0.05e-3
ENa = 0.050
EK = -0.080
Eh = -0.010
ECa = 0.140
Em = -0.070
vmin = -0.100
vmax = 0.100
vdivs = 3000
dV = (vmax-vmin)/vdivs
moose.Neutral('/library')
Chan = KDR_SChan('chan')
gate = moose.element('/library/chan/gateX')

plt.plot(np.arange(vmin,vmax+dV,dV),gate.tableA/gate.tableB, label='Inf')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the X gate of KDR channel')
plt.legend()
plt.show()
plt.plot(np.arange(vmin,vmax+dV,dV),np.ones(len(gate.tableB))/gate.tableB, label='tau')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Time constant (s)')
plt.title('tau for the X gate of KDR channel')
plt.legend()
plt.show()

#########Checking if everything matches##########################################
vmin = -100
vmax = 100
vdivs = 3000
dV = (vmax-vmin)/vdivs

#######KDR_S NEURON#####
vhalfn=13
a0n=0.02
zetan=-3
gmn=0.7
nmax=2
q10=1

v = np.arange(vmin,vmax+dV, dV)
alpn = np.exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
betn = np.exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius)))
qt=q10**((celsius-24)/10)
a = alpn
ninf = 1/(1+a)
taun = betn/(qt*a0n*(1+a))
taun[taun<nmax]=nmax

plt.plot(v, ninf, label='Inf Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the X gate of KDR_S channel')
plt.legend()
plt.show()

plt.plot(v, taun, label='tau Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Tau (ms)')
plt.title('Tau for the X gate of KDR_S channel')
plt.legend()
plt.show()

#########Na_S NEURON########
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
q10=3
Rg   = 0.01
Rd   = .03
qq   = 10
tq   = -55
thinf  = -50
qinf  = 2
vhalfs=-60
a0s=0.0003
zetas=12
gms=0.2
smax=10
vvh=-58
vvs=2
ar2=1
celsius = 34

v = np.arange(vmin,vmax+dV, dV)
alpv = 1/(1+np.exp((v-vvh)/vvs))
alps = np.exp(1.e-3*zetas*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))
bets = np.exp(1.e-3*zetas*gms*(v-vhalfs)*9.648e4/(8.315*(273.16+celsius)))

def trap0(v,th,a,q):
    if abs(v-th)>1e-6:
        return a * (v - th) / (1 - np.exp(-(v - th)/q))
    else:
        return a * q
        
qt=q10**((celsius-24)/10)
a = np.array([trap0(vm,tha,Ra,qa) for vm in v])
b = np.array([trap0(-vm,-tha,Rb,qa) for vm in v])
mtau = 1/(a+b)/qt
mtau[mtau<mmin]=mmin
minf = a/(a+b)

a = np.array([trap0(vm,thi1,Rd,qd) for vm in v])
b = np.array([trap0(-vm,-thi2,Rg,qg) for vm in v])
htau =  1/(a+b)/qt
htau[htau<hmin] = hmin
hinf = 1/(1+np.exp((v-thinf)/qinf))
c=alpv
sinf = c+ar2*(1-c)
taus = bets/(a0s*(1+alps))
taus[taus<smax] = smax

plt.plot(v, minf, label='Inf Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the X gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v, mtau, label='tau Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Tau (ms)')
plt.title('Tau for the X gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v, hinf, label='Inf Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the Y gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v, htau, label='tau Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Tau (ms)')
plt.title('Tau for the Y gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v, sinf, label='Inf Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the Z gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v, taus, label='tau Bianchi')
plt.xlabel('Membrane potential (mV)')
plt.ylabel('Tau (ms)')
plt.title('Tau for the Z gate of Na_S channel')
plt.legend()
plt.show()

#########

#####KDR MOOSE###########
exec(open('/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work/Bianchi2002_WT_CA1/CA1_WT.py').read())

vmin = -0.100
vmax = 0.100
vdivs = 3000
dV = (vmax-vmin)/vdivs
v = np.arange(vmin,vmax+dV,dV)

gateX = moose.element('/library/KDR_Schan/gateX')
plt.plot(v,gateX.tableA/gateX.tableB, label='Inf')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the X gate of KDR channel')
plt.legend()
plt.show()

plt.plot(v,np.ones(len(gateX.tableB))/gateX.tableB, label='tau')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Time constant (s)')
plt.title('tau for the X gate of KDR channel')
plt.legend()
plt.show()

###Na_S MOOSE###########
vmin = -0.100
vmax = 0.100
vdivs = 3000
dV = (vmax-vmin)/vdivs
v = np.arange(vmin,vmax+dV,dV)

gateX = moose.element('/library/Na_Schan/gateX')
plt.plot(v,gateX.tableA/gateX.tableB, label='Inf')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the X gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v,np.ones(len(gateX.tableB))/gateX.tableB, label='tau')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Time constant (s)')
plt.title('tau for the X gate of Na_S channel')
plt.legend()
plt.show()

gateY = moose.element('/library/Na_Schan/gateY')
plt.plot(v,gateY.tableA/gateY.tableB, label='Inf')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the Y gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v,np.ones(len(gateY.tableB))/gateY.tableB, label='tau')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Time constant (s)')
plt.title('tau for the Y gate of Na_S channel')
plt.legend()
plt.show()

gateZ = moose.element('/library/Na_Schan/gateZ')
plt.plot(v,gateZ.tableA/gateZ.tableB, label='Inf')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Inf (unitless)')
plt.title('Inf for the Z gate of Na_S channel')
plt.legend()
plt.show()

plt.plot(v,np.ones(len(gateZ.tableB))/gateZ.tableB, label='tau')
plt.xlabel('Membrane potential (V)')
plt.ylabel('Time constant (s)')
plt.title('tau for the Z gate of Na_S channel')
plt.legend()
plt.show()