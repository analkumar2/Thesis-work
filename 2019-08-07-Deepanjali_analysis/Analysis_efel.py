# exec(open('Deepanjali data/Analysis_efel.py').read())

# Extracts single feature using efel package


import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import os
import efel

efel.getFeatureNames()
filename = 'Cell 6 of 171117.abf'
seg_no = 17 #0 is -100pA, 4 is 0pA, 20 is 400pA. Now extrapolate
stim_start = 81.4 #in ms
stim_end = 581.4 #in ms
#Spike count threshold is set to -20 by default

Actualstim_start = seg_no*2000+stim_start
Actualstim_end = seg_no*2000+stim_end
data = efel.io.load_neo_file("Deepanjali data/WT step input cells/"+filename, stim_start=Actualstim_start, stim_end=Actualstim_end)
features = efel.getFeatureValues(data[0][seg_no], ['voltage_base', 'AP1_width', ])
print(features)

plt.plot(data[0][seg_no][0]['T'], data[0][seg_no][0]['V'])
plt.show()

###########Bulk##############
