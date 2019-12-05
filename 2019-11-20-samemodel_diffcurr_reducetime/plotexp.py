# exec(open('plotexp.py').read())


import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os

def expdata(Address, currentlevel):
    CurrAmp = int((currentlevel - (-100e-12))/25e-12)
    reader = AxonIO(filename=Address)
    Samprate = reader.get_signal_sampling_rate()
    seg = reader.read_block().segments[CurrAmp]
    Tdur = np.array(seg.t_stop - seg.t_start)
    return [np.linspace(0,Tdur+0,Samprate*Tdur), np.ravel(seg.analogsignals[0])*1e-3]

def plotexp(Address, currentlevel):
    T, V = expdata(Address, currentlevel)
    plt.plot(T,V, label=Address)
    plt.legend()
    plt.title(f'Current clamp at {currentlevel}A')
    plt.xlabel('Time (s)')
    plt.ylabel('Membrane potential (V)')
    plt.show()

if __name__ == '__main__':
    plotexp('Experimental recordings/cell 4 of 61016.abf', 150e-12)
