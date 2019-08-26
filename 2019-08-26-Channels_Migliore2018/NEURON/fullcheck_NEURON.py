# exec(open('fullcheck_NEURON.py').read())
from neuron import h, gui
import numpy as np
import matplotlib.pyplot as plt

Chan = 'na3'

v_init = -70

soma = h.Section(name='soma')
soma(0.5).v = -70
soma.diam = 20
soma.L = 20
soma.insert('pas')
soma.insert(Chan)
soma.insert('cacum')
soma(0.5).cacum.cai0 = 0.05e-3
h.cai0_ca_ion = 0.05e-3
h.cao0_ca_ion = 2
h.celsius = 32


stim = h.IClamp(soma(0.5))
stim.delay = 1000
stim.dur = 500
stim.amp = 0.5

# stim = h.SEClamp(soma(0.5))
# stim.rs = 0.05
# stim.dur1 = 1000
# stim.dur2 = 500
# stim.dur3 = 500
# stim.amp1 = -65
# stim.amp2 = 10
# stim.amp3 = -65

t_vec = h.Vector()             # Time stamp vector
v_vec = h.Vector()             # Membrane potential vector
cai_vec = h.Vector()
ichan_vec = h.Vector()
ihold_vec = h.Vector()
gate_vec = h.Vector()

t_vec.record(h._ref_t)
v_vec.record(soma(0.5)._ref_v)
cai_vec.record(soma(0.5)._ref_cai)
ichan_vec.record(soma(0.5).na_ion._ref_ina) # change needed here for different channels
# ihold_vec.record(stim._ref_i)
# gate_vec.record(soma(0.5).cal._ref_o)


h.tstop = 2000
h.v_init = v_init #Initialization screws everything as it sets the initial minf to set the v to v_init. Thus, gives diff results than MOOSE
h.finitialize()
h.run()

fig1, ax1 = plt.subplots()
ax1.plot(t_vec,v_vec,label = Chan + ' channel')
ax1.set_xlabel('Time (ms)')
ax1.set_ylabel('Membrane potential (mV)')
ax1.legend()
ax1.set_title(Chan + ' channel NEURON')

# fig2, ax2 = plt.subplots()
# ax2.plot(t_vec,cai_vec, label='internal calcium conc')
# ax2.set_xlabel('Time (ms)')
# ax2.set_ylabel('Concentration (mM)')
# ax2.legend()
# ax2.set_title(Chan + ' channel NEURON')

fig3, ax3 = plt.subplots()
ichan_vec = np.array(ichan_vec)
ax3.plot(t_vec,ichan_vec*soma(0.5).area()*1e-2, label='channel current') #soma(0.5).area()*1e-2 beacuse ica is in current density in mA/cm^2
ax3.set_xlabel('Time (ms)')
ax3.set_ylabel('Current (nA)')
ax3.legend()
ax3.set_title(Chan + ' channel NEURON')

# fig4, ax4 = plt.subplots()
# ax4.plot(t_vec,ihold_vec,label = 'Holding current')
# ax4.set_xlabel('Time (ms)')
# ax4.set_ylabel('Holding current (nA)')
# ax4.legend()
# # ax4.set_ylim(-2, 2)
# ax4.set_title(Chan + ' channel NEURON')

plt.show()
