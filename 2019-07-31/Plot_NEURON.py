#exec(open('Plot_NEURON.py').read())

import csv
import numpy as np
import matplotlib.pyplot as plt

Vtrace = []
time = []

with open('Combe2018/curr150.txt') as file:
    reader = csv.reader(file, delimiter = '\t')
    a = next(reader)
    a = next(reader)
    for row in reader:
        Vtrace.append(float(row[1]))
        time.append(float(row[0]))

plt.plot([t*1e-3 for t in time], [V*1e-3 for V in Vtrace])
plt.title('Combe2018 modelDB version current clamp of 150pA from 1 to 1.5s')
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.show()
