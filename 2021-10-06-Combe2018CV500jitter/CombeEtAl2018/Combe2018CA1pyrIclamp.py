## exec(open('Migliore2018CA1pyrIclamp.py').read())

from neuron import h,gui
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
sbn.set()

h('load_file("simplestim.hoc")')

v_vec = h.Vector()             # Membrane potential vector
t_vec = h.Vector()             # Time stamp vector
v_vec.record(h.soma[0](0.5)._ref_v)
t_vec.record(h._ref_t)

stim = h.IClamp(h.soma[0](0.5))
stim.delay = 1000
stim.amp = 0.150
stim.dur = 500
h.finitialize()
h.tstop = 2000
h.run()

plt.plot(np.array(t_vec)*1e-3,np.array(v_vec)*1e-3, label=f'{stim.amp*1000}pA')
plt.show()

stim.amp = 0.120
h.finitialize()
h.tstop = 2000
h.run()
plt.plot(np.array(t_vec)*1e-3,np.array(v_vec)*1e-3, label=f'{stim.amp*1000}pA')

stim.amp = 0.200
h.finitialize()
h.tstop = 2000
h.run()
plt.plot(np.array(t_vec)*1e-3,np.array(v_vec)*1e-3, label=f'{stim.amp*1000}pA')

plt.legend()
plt.title('Combe 2018 in-silico CA1 pyramidal neuron current clamp')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (V)')
plt.xlim(0.7,1.7)
plt.tight_layout()
plt.show()
