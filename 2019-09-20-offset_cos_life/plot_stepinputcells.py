# exec(open('plot_stepinputcells.py').read())

# Plots all the data one by one
# -100pA to 400pA
# Change CurrAmp

import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os


Direc = '../../Raw_data/Deepanjali_data/WT step input cells/' #Directory where all the files are stored
CurrAmp = 10 #0 is for -100pA, 1 for -75pA, and so on till 20 which is 400pA. 10 for 150pA
# flist = os.listdir(Direc)
flist = ['cell 4 of 61016.abf']

i = 0
print(os.listdir(Direc))
for filename in flist:
    i = i+1
    print(i, end = '\r')
    try: #Try is used to skip unsupported files that may be there in the folder
        reader = AxonIO(filename=Direc+filename)
        Samprate = reader.get_signal_sampling_rate()
        seg = reader.read_block().segments[CurrAmp]
        Tdur = np.array(seg.t_stop - seg.t_start)
        plt.plot(np.linspace(0,Tdur+0,Samprate*Tdur), seg.analogsignals[0], label = f'{filename}')
        # plt.title('WT 150pA current injection for 0.5s')
        # plt.xlim([0,1])
        plt.ylim([-75,50])
        plt.xlabel('Time (s)')
        plt.ylabel('Membrane potential (mV)')
        plt.legend()
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())
        plt.show()
    except:
        pass

# plt.show()
