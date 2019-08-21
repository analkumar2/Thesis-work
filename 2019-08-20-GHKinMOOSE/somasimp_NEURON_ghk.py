# exec(open('somasimp_NEURON_ghk.py').read())
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

v_init=-70

soma = h.Section(name='soma')
soma.insert('pas')
soma.insert('simp')
# soma.insert('cacum')
# CaD = h.cacum(soma(0.5))

stim = h.IClamp(soma(0.5))
stim.delay = 1000
stim.dur = 500
stim.amp = 100

v_vec = h.Vector()             # Membrane potential vector
t_vec = h.Vector()             # Time stamp vector
cai_vec = h.Vector()
v_vec.record(soma(0.5)._ref_v)
t_vec.record(h._ref_t)
cai_vec.record(soma(0.5)._ref_cai)

h.tstop = 2000
h.v_init = v_init
h.run()

plt.plot(t_vec,v_vec,label = 'ghk simple channel')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.legend()
plt.title('ghk simple channel')
plt.show()
