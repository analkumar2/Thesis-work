# exec(open('somasimp_NEURON_ghk.py').read())
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

v_init = -70

soma = h.Section(name='soma')
soma(0.5).v = -70
soma.diam = 20
soma.L = 20
soma.insert('pas')
soma.insert('simp')
soma.insert('cacum')
h.cai0_ca_ion = 0.05e-3
h.cao0_ca_ion = 2
h.celsius = 32


stim = h.IClamp(soma(0.5))
stim.delay = 1000
stim.dur = 500
stim.amp = 1

t_vec = h.Vector()             # Time stamp vector
v_vec = h.Vector()             # Membrane potential vector
cai_vec = h.Vector()
isimp_vec = h.Vector()
m_vec = h.Vector()

t_vec.record(h._ref_t)
v_vec.record(soma(0.5)._ref_v)
cai_vec.record(soma(0.5)._ref_cai)
isimp_vec.record(soma(0.5).simp._ref_ica)
m_vec.record(soma(0.5).simp._ref_m)


h.tstop = 2000
h.v_init = v_init #Initialization screws everything as it sets the initial minf to set the v to v_init. Thus, gives diff results than MOOSE
h.finitialize()
h.run()

fig1, ax1 = plt.subplots()
ax1.plot(t_vec,v_vec,label = 'ghk simple channel')
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Membrane potential (mV)')
ax1.legend()
ax1.set_title('ghk simple channel')

# fig2, ax2 = plt.subplots()
# ax2.plot(t_vec,cai_vec, label='internal calcium conc')
# ax2.set_xlabel('Time (ms)')
# ax2.set_ylabel('Concentration (mM)')
# ax2.legend()
# ax2.set_title('ghk simple channel')
#
# fig3, ax3 = plt.subplots()
# ax3.plot(t_vec,isimp_vec, label='simp channel current')
# ax3.set_xlabel('Time (ms)')
# ax3.set_ylabel('Current (nA)')
# ax3.legend()
# ax3.set_title('ghk simple channel')
plt.show()
