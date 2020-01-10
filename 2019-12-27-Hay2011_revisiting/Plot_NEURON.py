#exec(open('Plot_NEURON.py').read())

import csv
import numpy as np
import matplotlib.pyplot as plt
import os

foldername=os.path.basename(os.getcwd())
filename = '150pAekm75.txt'
Vtrace = []
time = []

with open(f'../../Output/{foldername}/{filename}') as file:
    reader = csv.reader(file, delimiter = '\t')
    a = next(reader)
    a = next(reader)
    for row in reader:
        Vtrace.append(float(row[1]))
        time.append(float(row[0]))

plt.plot([t*1e-3 for t in time], [V*1e-3 for V in Vtrace], label='Hay2011')
plt.title('Hay2011 L5PC model vs CA1PC experiment. 150pA injection')
# plt.title('500pA vs 1000pA in Hay2011 L5PC model')
plt.xlabel('Time (s)')
plt.ylabel('Membrane potential (V)')
plt.legend()
# plt.show()
