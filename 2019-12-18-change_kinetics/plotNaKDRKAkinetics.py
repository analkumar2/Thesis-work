# exec(open('plotNaKDRKAkinetics.py').read())

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import moose

####////////////// Migliore 2018 //////////////////////////////////////////////


# exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore2018).py').read())
# exec(open('../../Compilations/Kinetics/K_DR_Chan_(Migliore2018).py').read())
# exec(open('../../Compilations/Kinetics/K_A_Chan_(Migliore2018).py').read())
# moose.Neutral('/library')
# Na_Chan('Na_Chan')
# K_DR_Chan('K_DR_Chan')
# K_A_Chan('K_A_Chan')

# v = np.linspace(-0.1,0.1,3000)

# Mig2018minf_Na = moose.element('/library/Na_Chan/gateX').tableA/moose.element('/library/Na_Chan/gateX').tableB
# Mig2018mtau_Na = 1/moose.element('/library/Na_Chan/gateX').tableB

# Mig2018hinf_Na = moose.element('/library/Na_Chan/gateY').tableA/moose.element('/library/Na_Chan/gateY').tableB
# Mig2018htau_Na = 1/moose.element('/library/Na_Chan/gateY').tableB

# Mig2018minf_K_DR = moose.element('/library/K_DR_Chan/gateX').tableA/moose.element('/library/K_DR_Chan/gateX').tableB
# Mig2018mtau_K_DR = 1/moose.element('/library/K_DR_Chan/gateX').tableB

# Mig2018minf_K_A = moose.element('/library/K_A_Chan/gateX').tableA/moose.element('/library/K_A_Chan/gateX').tableB
# Mig2018mtau_K_A = 1/moose.element('/library/K_A_Chan/gateX').tableB

# Mig2018hinf_K_A = moose.element('/library/K_A_Chan/gateY').tableA/moose.element('/library/K_A_Chan/gateY').tableB
# Mig2018htau_K_A = 1/moose.element('/library/K_A_Chan/gateY').tableB

# plt.figure()
# ax1 = plt.subplot(111)

# ax1.plot(v, Mig2018minf_Na, 'b', label='Mig2018minf_Na')
# ax1.plot(v, Mig2018hinf_Na, 'g', label='Mig2018hinf_Na')
# ax1.plot(v, Mig2018minf_K_DR, 'r', label='Mig2018minf_K_DR')
# ax1.plot(v, Mig2018minf_K_A, 'c', label='Mig2018minf_K_A')
# ax1.plot(v, Mig2018hinf_K_A, 'm', label='Mig2018hinf_K_A')

# ax1.legend()

# plt.figure()
# ax2 = plt.subplot(111)
# ax3 = ax2.twinx()

# ax2.plot(v, Mig2018mtau_Na, 'b', label='Mig2018mtau_Na')
# ax2.plot(v, Mig2018htau_Na, 'g', label='Mig2018htau_Na')

# ax2.plot(v, Mig2018mtau_K_DR, 'r', label='Mig2018mtau_K_DR')

# ax2.plot(v, Mig2018mtau_K_A, 'c', label='Mig2018mtau_K_A')
# ax2.plot(v, Mig2018htau_K_A, 'm', label='Mig2018htau_K_A')

# ax2.legend()
# ax3.legend()

# plt.show()


## ///////////////////////////////////////////////////////////////////////


## //// Hay2011 //////////////////////////////////////////////////////////

exec(open('../../Compilations/Kinetics/Na_T_Chan_(Hay2011).py').read())
exec(open('../../Compilations/Kinetics/K_DR_Chan_(Hay2011).py').read())
moose.Neutral('/library')
Na_T_Chan('Na_T_Chan')
K_DR_Chan('K_DR_Chan')

v = np.linspace(-0.1,0.1,3000)

Hay2011minf_Na_T = moose.element('/library/Na_T_Chan/gateX').tableA/moose.element('/library/Na_T_Chan/gateX').tableB
Hay2011mtau_Na_T = 1/moose.element('/library/Na_T_Chan/gateX').tableB

Hay2011hinf_Na_T = moose.element('/library/Na_T_Chan/gateY').tableA/moose.element('/library/Na_T_Chan/gateY').tableB
Hay2011htau_Na_T = 1/moose.element('/library/Na_T_Chan/gateY').tableB

Hay2011minf_K_DR = moose.element('/library/K_DR_Chan/gateX').tableA/moose.element('/library/K_DR_Chan/gateX').tableB
Hay2011mtau_K_DR = 1/moose.element('/library/K_DR_Chan/gateX').tableB


plt.figure()
ax1 = plt.subplot(111)

ax1.plot(v, Hay2011minf_Na_T, 'b', label='Hay2011minf_Na_T')
ax1.plot(v, Hay2011hinf_Na_T, 'g', label='Hay2011hinf_Na_T')
ax1.plot(v, Hay2011minf_K_DR, 'r', label='Hay2011minf_K_DR')

ax1.legend()

plt.figure()
ax2 = plt.subplot(111)
ax3 = ax2.twinx()

ax2.plot(v, Hay2011mtau_Na_T, 'b', label='Hay2011mtau_Na_T')
ax2.plot(v, Hay2011htau_Na_T, 'g', label='Hay2011htau_Na_T')

ax2.plot(v, Hay2011mtau_K_DR, 'r', label='Hay2011mtau_K_DR')

ax2.legend()
ax3.legend()

plt.show()

## ///////////////////////////////////////////////////////////////////////