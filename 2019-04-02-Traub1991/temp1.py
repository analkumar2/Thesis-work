import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

# rdes = rd.rdesigneur(
	# chanProto = [
		# ['make_K_AHP()', 'K_AHP'],
		# ['make_Ca_conc()', 'Ca_conc']
	# ],
	# plotList = [
		# ['soma', '1', '.', 'Vm', 'Membrane potential']
	# ]
	
	
# )

# rdes.buildModel()
# moose.reinit()
# moose.start( 2 )
# rdes.display()


# solve calcium current to concentration

dt = 0.001
B = 17.402e12/40000
tau = 0.01333
I = 0.6e-9

Calist = [0]
Ca = 0
for i in np.arange(dt,2,dt):
	Ca = dt*(B*I - Ca/tau) + Ca
	Calist.append(Ca)

plt.plot(np.arange(0,2,dt), Calist)
plt.show()



import moose
import pylab
import rdesigneur as rd
rdes = rd.rdesigneur(
    cellProto = [['somaProto', 'soma', 12e-6, 12e-6]],
    chanProto = [
        ['make_Na()', 'Na'],
        ['make_K_DR()', 'K_DR'],
        ['make_K_A()', 'K_A' ],
        ['make_Ca()', 'Ca' ],
        ['make_Ca_conc()', 'Ca_conc' ]
    ],
    # Some changes to the default passive properties of the cell.
    passiveDistrib = [['soma', 'CM', '0.03', 'Em', '-0.06']],
    chanDistrib = [
        ['Na', 'soma', 'Gbar', '300' ],
        ['K_DR', 'soma', 'Gbar', '250' ],
        ['K_A', 'soma', 'Gbar', '200' ],
        ['Ca_conc', 'soma', 'tau', '0.0333' ],
        ['Ca', 'soma', 'Gbar', '40' ]
    ],
    # Give a + pulse from 5 to 7s, and a - pulse from 20 to 21.
    stimList = [['soma', '1', '.', 'inject', '((t>1 && t<1.2) * 1.0e-11' ]],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Membrane potential'],
		# ['soma', '1', 'Ca', 'Ik', 'Calcium current'],
		['soma', '1', 'Ca_conc', 'Ca', 'Calcium concentration'],
		['soma', '1', 'Ca', 'Gk', 'Calcium conductance']
    ],
)

rdes.buildModel()
moose.reinit()
moose.start( 2 )

rdes.display()



#plot state variable rates
EREST = -0.06
def plotstatevarrate(A,B,C,D,E):
	V = np.linspace(EREST-0.05, EREST+0.15, 3000)
	rate = (A + B*V)/(C + np.exp((V + D)/E))
	return plt.plot(V,rate)
	plt.show()
    
ratecon = moose.element('/model/elec/soma/Ca_Chan/gateX')
[moose.delete(x) for x in ['/model', '/library']]
execfile('traub1991.py')




#To find how channels work

    
rdes = rd.rdesigneur(
    chanProto = [['Channelprotos.NaChan()','simp_Chan']],
    chanDistrib = [['simp_Chan', 'soma', 'Gbar', '1']],
    plotList = [['soma', '1', '.', 'Vm', 'Membrane potential'],
        ['soma', '1', 'simp_Chan', 'Gk', 'simple channel conductance']],
    )
    
rdes.buildModel()
moose.reinit()

data = moose.Neutral('/data')
xgateA = moose.Table('/data/xgateA')
xgateB = moose.Table('/data/xgateB')

xgate = moose.element('/model/elec/soma/simp_Chan/gateX')

moose.connect(xgateA, 'requestOut', xgate, 'getA')
moose.connect(xgateB, 'requestOut', xgate, 'getB')

moose.start( 2 )

plt.figure()
plt.plot(xgateA.vector)
plt.plot(xgateB.vector)
plt.show()

rdes.display()


#Checking useCOncentration

def findThickness(R, surfA, B):
    thick = R - np.sqrt(R**2 - 2*R/(B*F*surfA))
    return(str(thick))

def surfA(Diameter, length):
    return np.pi*Diameter*length

[moose.delete(x) for x in ['/model', '/library']]

rdes = rd.rdesigneur(
    chanProto = [
        ['Channelprotos.simpChan()', 'simp_Chan'],
        ['Channelprotos.CaConc()', 'Ca_conc'],
        ],
    chanDistrib = [
        ['Ca_conc', 'soma', 'thick', findThickness(2.5e-4, surfA(0.0005, 0.0005), 17.402e12)],
        ['simp_Chan', 'soma', 'Gbar', '1'],
        ],
    stimList = [
        ['soma', '1', '.', 'inject', '(t>1 && t<=1.2) ? 3e-9 : 0' ],
    ],
    plotList = [
        ['soma', '1', '.', 'Vm', 'Membrane potential'],
        ['soma', '1', 'simp_Chan', 'Gk', 'Kc channel conductance'],
        ['soma', '1', 'Ca_conc', 'Ca', 'Calcium concentration'],
        ],
    )
rdes.buildModel()
moose.reinit()
moose.start( 2 )
rdes.display()


#Calcium channel current #Gbar*s**2*r*(V-Vca)
def alphas(V):
    return 1.6e3/(1 + np.exp((V+0.06-65e-3)*-72))
def betas(V):
    return 20*(V+0.06-51.1e-3)/(-1 + np.exp((V+0.06-51.1e-3)/5e-3))

def alphar(V):
    if V<=-0.06:
        return 5
    else:
        return 5*np.exp(-50*(0.06-V))
def betar(V):
    if V<=-0.06:
        return 0
    else:
        return 5-alphar(V)
    
        

dt = 1e-6
dx = 1e-6
Gbar = 40*3320e-12
def s(prevs, volt):
    return prevs + dt*( alphas(volt)*(1-prevs) - betas(volt)*prevs )

def r(prevr, volt):
    return prevr + dt*( alphar(volt)*(1-prevr) - betar(volt)*prevr )
    

prevs = 0
prevr = 0    
for v in [-0.075, -0.06, -0.02, 0, 0.02, 0.05]:
    curr = []
    for i in range(3000000):
        curr.append( Gbar*s(prevs, v)**2*r(prevr, v)*(v-0.140) )
        prevs = s(prevs, v)
        prevr = r(prevr, v)
    print np.min(curr)


#PLotting fI curves
plt.plot([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], [273.0/15, 393.0/15, 476.0/15, 598.0/15, 702.0/15, 833.0/15], 'o-')
plt.title('fI curve for action potentials')         393.0/15, 476.0/1
plt.ylabel('Action potential frequency (Hz)')
plt.xlabel('Somatic injected current (nA)')
plt.show()

plt.plot([0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.275, 0.3], [4.0/15, 5.0/15, 6.0/15, 7.0/15, 8.0/15, 9.0/15, 10.0/15, 11.0/15, 11.0/15], 'o-')
plt.title('fI curve for bursting')
plt.ylabel('Burst frequency (Hz)')
plt.xlabel('SOmatic injected current (nA)')
plt.ylim([0.0, 1.0])
plt.show()

