# exec(open('Izhitemplate.py').read())
# Tried using USEION to communicate. Works but the weights and amp has to be setup to very high values. Probably just use the pointer version
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

v_init=-70

dummyE1 = h.Section(name='dummyE1') #dummy section to insert Izhikevich

izhiE1 = h.Izhi2003a(0.5,sec=dummyE1) #insert Izhikevich
izhiE1.V = v_init

AMPAR = h.Exp2SynAMPAR(dummyE1(0.5)) #insert AMPAR synapse
AMPAR.tau1 = 3.5
AMPAR.tau2 = 21
AMPAR.e = 0

GABAR = h.Exp2SynGABAR(dummyE1(0.5)) #insert GABAR synapse
GABAR.tau1 = 6.4
GABAR.tau2 = 40
GABAR.e = -80

stim = h.IClampizhi(dummyE1(0.5)) #insert IClamp izhikevich version
stim.delay = 2000
stim.dur = 500
stim.amp = 10000

ArtinputE = h.NetStim() #Excitatory neuron input
ArtinputE.interval = 20
ArtinputE.number = 5
ArtinputE.start = 500

ArtinputI = h.NetStim() #Inhibitory neuron input
ArtinputI.interval = 20
ArtinputI.number = 5
ArtinputI.start = 1000

con_ArtinputE_AMPAR = h.NetCon(ArtinputE, AMPAR) #Connecting excitatory input to AMPAR
con_ArtinputE_AMPAR.weight[0] = 150
con_ArtinputI_GABAR = h.NetCon(ArtinputI, GABAR) #Connecting inhibitory input to GABAR
con_ArtinputI_GABAR.weight[0] = 100

h.setpointer(izhiE1._ref_V, 'Vizhi', AMPAR) #So that AMPAR uses Izhikevich voltage V instead of the cell's voltage
h.setpointer(izhiE1._ref_V, 'Vizhi', GABAR) #So that GABAR uses Izhikevich voltage V instead of the cell's voltage
# h.setpointer(stim._ref_icla, 'Ic', izhiE1) #So that
# h.setpointer(AMPAR._ref_iam, 'Iampar', izhiE1)
# h.setpointer(GABAR._ref_iga, 'Igabar', izhiE1)

v_vec = h.Vector()             # Membrane potential vector
t_vec = h.Vector()             # Time stamp vector
v_vec.record(izhiE1._ref_V)
t_vec.record(h._ref_t)

h.finitialize(v_init)
h.tstop = 3000
h.v_init = v_init
h.run()

plt.plot(t_vec,v_vec)
plt.show()