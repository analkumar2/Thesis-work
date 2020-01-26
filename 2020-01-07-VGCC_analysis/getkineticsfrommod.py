## exec(open('getkineticsfrommod.py').read())
## sudo nrnivmodl '../../Compilations/ExistingModels/Poolos2002/lamotrigine/na3n.mod'

from neuron import h,gui
import numpy as np
import matplotlib.pyplot as plt
import brute_curvefit

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
	# alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

h.celsius = 32
soma = h.Section(name='soma')
soma.insert('pas')
# soma.insert('kdr')
soma.insert('cal')
soma.insert('can')
soma.insert('cat')
# soma.insert('calH')
# soma.insert('car')
# soma.insert('cat')
# soma.insert('car_mag')
# soma.insert('Nap_Et2')
# soma.insert('NaTg')
# soma.insert('nax')
soma.insert('hd')
soma.insert('kmb')
soma.insert('cagk')

h.cai0_ca_ion = 100e-5

clamp = h.SEClamp(soma(0.5))
clamp.rs = 1e-3 # series resistance should be much smaller than input resistance of the cell
clamp.dur1 = 1e9

cmd = h.Vector(np.linspace(-100,100,3000))
cmd.play(clamp._ref_amp1, 2000/3000)
h.v_init = -100


minf_vec = h.Vector()
mtau_vec = h.Vector()             # Membrane potential vector
v_vec = h.Vector()
t_vec = h.Vector()             # Time stamp vector
v_vec.record(soma(0.5)._ref_v)
minf_vec.record(soma(0.5).cagk._ref_oinf)
mtau_vec.record(soma(0.5).cagk._ref_tau)
# minf_vec.record(soma(0.5).cat._ref_hinf)
# mtau_vec.record(soma(0.5).cat._ref_tauh)
# minf_vec.record(soma(0.5).car._ref_inf[1])
# mtau_vec.record(soma(0.5).car._ref_tau[1])
# minf_vec.record(soma(0.5).nax._ref_minf)
# mtau_vec.record(soma(0.5).nax._ref_mtau)
t_vec.record(h._ref_t)


h.finitialize()
h.tstop = 2000
h.run()

plt.plot(t_vec, v_vec)
plt.show()

def wr_ChanGate_inf(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
	return ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F)[0]

def wr_ChanGate_tau(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
	return ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F)[1]

paramfitted_inf, error = brute_curvefit.brute_scifit(wr_ChanGate_inf, np.array(v_vec)*1e-3, np.array(minf_vec), restrict=[[-0.1,-0.1, -0.1,0,0,0,0,0],[0.1,0.1, 0.1,0.1,0.1,0.5,0.5,1]])

plt.plot(v_vec, minf_vec, label='ori')
plt.plot(v_vec, wr_ChanGate_inf(np.array(v_vec)*1e-3, *paramfitted_inf), label='fitted')
print(paramfitted_inf)
plt.legend()
plt.show()

# paramfitted_tau, error = brute_curvefit.brute_scifit(wr_ChanGate_tau, np.array(v_vec)*1e-3, np.array(mtau_vec)*1e-3, restrict=[[-0.1,-0.1, -0.1,0,0,0,0,0],[0.1,0.1, 0.1,0.5,0.5,0.5,0.5,1]], ntol = 10000, returnnfactor = 0.001)

# plt.plot(v_vec, np.array(mtau_vec)*1e-3, label='ori')
# plt.plot(v_vec, wr_ChanGate_tau(np.array(v_vec)*1e-3, *paramfitted_tau), label='fitted')
# print(paramfitted_tau)
# plt.legend()
# plt.show()