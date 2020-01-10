jarjar = input('haha1')
import sys
import numpy as np
import quantities as pq
# import matplotlib.pyplot as plt
# from neo.io import AxonIO
import os
jarjar = input('haha2')

# def expdata(Address, currentlevel):
#     CurrAmp = int((currentlevel - (-100e-12))/25e-12)
#     reader = AxonIO(filename=Address)
#     Samprate = reader.get_signal_sampling_rate()
#     seg = reader.read_block(signal_group_mode='split-all').segments[CurrAmp]
#     Tdur = np.array(seg.t_stop - seg.t_start)
#     return [np.linspace(0,Tdur+0,int(Samprate*Tdur)), np.array(np.ravel(seg.analogsignals[0]))*1e-3]


# def plotexp(Address, currentlevel):
#     T, V = expdata(Address, currentlevel)
#     plt.plot(T,V, label=Address)
#     plt.legend()
#     plt.title(f'Current clamp at {currentlevel}A')
#     plt.xlabel('Time (s)')
#     plt.ylabel('Membrane potential (V)')
#     plt.show()

# jarjar = input('haha3')
# if len(sys.argv)==1:
# 	print(sys.argv)
# 	currentlevel = float(input('currentlevel'))
# 	print(currentlevel)
# 	plotexp(sys.argv[0], int(currentlevel))
# elif len(sys.argv)==2:
# 	print(sys.argv)
# 	currentlevel = float(input('currentlevel'))
# 	print(currentlevel)
# 	plotexp(sys.argv[1], currentlevel)
