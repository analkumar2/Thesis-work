#exec(open('Combe2018/Plot_NEURON.py').read())

import csv
import numpy as np
import matplotlib.pyplot as plt

Vtrace = []
time = []

with open('Combe2018/110pA_combe2018.txt') as file:
    reader = csv.reader(file, delimiter = '\t')
    a = next(reader)
    a = next(reader)
    for row in reader:
        Vtrace.append(float(row[1]))
        time.append(float(row[0]))

plt.plot(time, Vtrace)
plt.show()
