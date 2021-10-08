## exec(open('Migliore2018CA1pyrIclamp.py').read())

from neuron import h,gui
import numpy as np
import matplotlib.pyplot as plt

h('load_file("CombeEtAl2018/simplestim.hoc")')

v_vec = h.Vector()             # Membrane potential vector
t_vec = h.Vector()             # Time stamp vector
v_vec.record(h.soma[0](0.5)._ref_v)
t_vec.record(h._ref_t)

stim = h.IClamp(h.soma[0](0.5))

samprate = 20000
totalsec = 2
curr = np.zeros(int(samprate*totalsec))
curr[int(1*samprate):int(1.9*samprate)] = 150e-12
noise = np.random.normal(0,20e-12,int(samprate*totalsec))
curr = curr + noise

currh = h.Vector(curr*1e9)
th = h.Vector(np.linspace(0,totalsec, int(samprate*totalsec)) * 1e3)

stim.delay = 0
stim.dur = 1e9
currh.play(stim._ref_amp, th, 1)

h.finitialize()
h.tstop = 2000
h.run()
plt.plot(np.array(t_vec)*1e-3,np.array(v_vec)*1e-3, label=f'{stim.amp*1000}pA')


def main(Channame):
	h('load_file("CombeEtAl2018/simplestim.hoc")')
    gbarratio = [0,0.1,0.2,0.5,0.75,0.9,1,1.1,1.5,2,3,5,10]
    # gbarratio = [0.9,1,1.1]

    #####
    if Channame == 'Na_Chan':
        CV500_Na = []
        jit_Na = []

        for i in range(len(gbarratio)):
            print('Na_Chan', end='\t')
            print(gbarratio[i], end='\t')
            CV500, jit = calcCVjit()
            CV500_Na.append(CV500)
            jit_Na.append(jit)
            ff = open('CVvsMean_True_channels4_jitCV.py', 'a+')
            ff.write(f'jit_Na = {jit_Na} \n')
            ff.write(f'CV500_Na = {CV500_Na} \n \n')
            ff.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('Chan', type=str)
    args = parser.parse_args()
    main(args.Chan)
    #fig = pickle.load(open('CVvsMean.pkl', 'rb'))
    #plt.show()
