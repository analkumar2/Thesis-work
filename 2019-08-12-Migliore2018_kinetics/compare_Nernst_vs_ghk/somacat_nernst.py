from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

v_init=-70

soma = h.Section(name='soma')
soma.insert('pas')
soma.insert('cat_nernst')

stim = h.IClamp(soma(0.5))
stim.delay = 100
stim.dur = 100
stim.amp = 100

v_vec = h.Vector()             # Membrane potential vector
t_vec = h.Vector()             # Time stamp vector
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)

h.tstop = 500
h.v_init = v_init
h.run()

plt.plot(t_vec,v_vec,label = 'nernst T-type calcium channel')
# plt.show()
