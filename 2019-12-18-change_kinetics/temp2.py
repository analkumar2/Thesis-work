# exec(open('temp2.py').read())

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import moose

exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore2018).py').read())
exec(open('../../Compilations/Kinetics/K_DR_Chan_(Migliore2018).py').read())
exec(open('../../Compilations/Kinetics/K_A_Chan_(Migliore2018).py').read())
moose.Neutral('/library')
Na_Chan('Na_Chan')
K_DR_Chan('K_DR_Chan')
K_A_Chan('K_A_Chan')

v = np.linspace(-0.1,0.1,3000)

Mig2018minf_Na = moose.element('/library/Na_Chan/gateX').tableA/moose.element('/library/Na_Chan/gateX').tableB
Mig2018mtau_Na = 1/moose.element('/library/Na_Chan/gateX').tableB

Mig2018hinf_Na = moose.element('/library/Na_Chan/gateY').tableA/moose.element('/library/Na_Chan/gateY').tableB
Mig2018htau_Na = 1/moose.element('/library/Na_Chan/gateY').tableB

Mig2018minf_K_DR = moose.element('/library/K_DR_Chan/gateX').tableA/moose.element('/library/K_DR_Chan/gateX').tableB
Mig2018mtau_K_DR = 1/moose.element('/library/K_DR_Chan/gateX').tableB

Mig2018minf_K_A = moose.element('/library/K_A_Chan/gateX').tableA/moose.element('/library/K_A_Chan/gateX').tableB
Mig2018mtau_K_A = 1/moose.element('/library/K_A_Chan/gateX').tableB

Mig2018hinf_K_A = moose.element('/library/K_A_Chan/gateY').tableA/moose.element('/library/K_A_Chan/gateY').tableB
Mig2018htau_K_A = 1/moose.element('/library/K_A_Chan/gateY').tableB

plt.figure()
ax1 = plt.subplot(111)

ax1.plot(v, Mig2018minf_Na, label='Mig2018minf_Na')
ax1.plot(v, Mig2018hinf_Na, label='Mig2018hinf_Na')
ax1.plot(v, Mig2018minf_K_DR, label='Mig2018minf_K_DR')
ax1.plot(v, Mig2018minf_K_A, label='Mig2018minf_K_A')
ax1.plot(v, Mig2018hinf_K_A, label='Mig2018hinf_K_A')

plt.figure()
ax2 = plt.subplot(111)
ax3 = ax2.twinx()

ax2.plot(v, Mig2018mtau_Na, label='Mig2018mtau_Na')
ax2.plot(v, Mig2018htau_Na, label='Mig2018htau_Na')

ax3.plot(v, Mig2018mtau_K_DR, label='Mig2018mtau_K_DR')

ax3.plot(v, Mig2018mtau_K_A, label='Mig2018mtau_K_A')
ax3.plot(v, Mig2018htau_K_A, label='Mig2018htau_K_A')

plt.legend()
plt.show()


