# exec(open('Bianchi2012/Analysis.py').read())
import numpy as np
import quantities as pq
import matplotlib.pyplot as plt
from neo.io import AxonIO
import csv
import os
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

filename = '300pA_0.0BK.csv'
v = []
t = []
with open(f'Bianchi2012/{filename}') as file:
    reader = csv.reader(file,delimiter = '\t')
    a = next(reader)
    a = next(reader)
    for row in reader:
        v.append(float(row[1]))
        t.append(float(row[0]))

plt.plot(t,v, label=filename)
plt.xlabel('time (ms)')
plt.xlim(750,1750)
plt.ylabel('Potential (mV)')
plt.title(filename)
plt.legend()
plt.show()
