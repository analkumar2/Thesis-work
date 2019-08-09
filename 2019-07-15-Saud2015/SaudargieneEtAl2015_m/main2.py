# exec(open('main2.py').read())

from neuron import h, gui
from matplotlib import pyplot as plt

h.load_file("nrngui.hoc")

h.load_file("BasalPath.hoc")
h.load_file("ObliquePath.hoc")
h.load_file("randomLocation.hoc")

h.load_file("pyramidalNeuron.hoc")

# objectvar cell
cell = h.PyramidalCell()

# // access cell.soma
# // objref stim
# // stim = new IClamp(cell.soma[0](0.5))


stim = h.IClamp(0.5,cell.soma[1])
stim.amp = 0.400
stim.delay = 1000
stim.dur = 500

v_vec = h.Vector()
v_vec.record(cell.soma[1](0.5)._ref_v)
t_vec = h.Vector()
t_vec.record(h._ref_t)


h.tstop = 2000.0
h.v_init = -70
h.run()

plt.plot(t_vec,v_vec)
plt.title('400pA current injection to soma')
plt.text(0.5,0.5, 'Current injected from 1000ms to 1500ms')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane potential (mV)')
plt.show()
