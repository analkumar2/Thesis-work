# exec(open('Izhitemplate.py').read())
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

v_init=-70

dummyE1 = h.Section(name='dummyE1')
dummyE1.v = v_init

izhiE1 = h.Izhi2003a(0.5,sec=dummyE1) #insert Izhikevich
izhiE1.V = v_init

AMPAR = h.ExpSyn(dummyE1(0.5)) #insert AMPAR synapse
AMPAR.tau = 2
AMPAR.e = 0

GABAR = h.ExpSyn(dummyE1(0.5)) #insert GABAR synapse
GABAR.tau = 50
GABAR.e = -80

stim = h.IClamp(dummyE1(0.5)) #inser IClamp
stim.delay = 2000
stim.dur = 500
stim.amp = 5

ArtinputE = h.NetStim()
ArtinputE.interval = 1
ArtinputE.number = 5
ArtinputE.start = 500

ArtinputI = h.NetStim()
ArtinputI.interval = 1
ArtinputI.number = 5
ArtinputI.start = 1000

con_ArtinputE_AMPAR = h.NetCon(ArtinputE, AMPAR)
con_ArtinputE_AMPAR.weight[0] = 1
con_ArtinputE_GABAR = h.NetCon(ArtinputI, GABAR)
con_ArtinputE_GABAR.weight[0] = 1

h.setpointer(stim._ref_i, 'Ic', izhiE1)
h.setpointer(AMPAR._ref_i, 'Iampar', izhiE1)
h.setpointer(GABAR._ref_i, 'Igabar', izhiE1)

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
