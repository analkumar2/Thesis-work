#exec(open('temp1.py').read())

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
mpl.rcParams["savefig.directory"] = '/mnt/c/Analkumar2/Study/Biology/Neuroscience/2018 - 23 PhD Thesis/Thesis work'

T = 33+273.15
z = 1
R = 8.314
F = 96485
Ko = 140
v = np.linspace(-0.100,0.100, 10000)
Er1 = -98.82e-3
Er2 = -98.79e-3

fbyP1 = z*z*F*F*Ko/R/T*v*(np.exp(z*F/R/T*(v-Er1))-1)/(np.exp(z*F*v/R/T)-1)
fbyP2 = z*z*F*F*Ko/R/T*v*(np.exp(z*F/R/T*(v-Er2))-1)/(np.exp(z*F*v/R/T)-1)
# Pm = 1.25e-7

plt.plot(v,fbyP1, label='Er=-98.82mV')
plt.plot(v,fbyP2, label='Er=-98.79mV')
plt.legend()
plt.xlabel('Potential')
plt.ylabel('Flux/Permeability')
plt.show()
