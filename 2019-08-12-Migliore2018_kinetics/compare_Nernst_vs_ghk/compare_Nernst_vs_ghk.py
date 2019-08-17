## THe problem with direct comparison is that ghk uses current flux. Nernst uses current.Incomplete code.



import matplotlib.pyplot as plt
import numpy as np

# Supoose G = 1S constant, we plot I/g vs grpah for cai = 0.05e-3, 0.1e-3, 0.5e-3 and 1e-3.
z = 2
F = 96485
T = 306
R = 8.314
cao = 2
area =


cai = [0.05e-3,0.1e-3,0.5e-3,1e-3]
Vm = np.linspace(-0.100,0.040,1000)
INernst = (Vm - 0.01318*np.log(cao/0.05e-3))/
Ighk005 = Vm*(z**2)*(F**2/R/T)*(cai[0] - cao*np.exp(-z*Vm*F/R/T))/(1-np.exp(-z*Vm*F/R/T))
Ighk010 = Vm*(z**2)*(F**2/R/T)*(cai[1] - cao*np.exp(-z*Vm*F/R/T))/(1-np.exp(-z*Vm*F/R/T))
Ighk050 = Vm*(z**2)*(F**2/R/T)*(cai[2] - cao*np.exp(-z*Vm*F/R/T))/(1-np.exp(-z*Vm*F/R/T))
Ighk100 = Vm*(z**2)*(F**2/R/T)*(cai[3] - cao*np.exp(-z*Vm*F/R/T))/(1-np.exp(-z*Vm*F/R/T))

plt.plot(Vm,INernst, label = "INernst")
plt.plot(Vm, Ighk005, label = "Ighk0.05e-3mM")
plt.plot(Vm, Ighk010, label = "Ighk0.10e-3mM")
plt.plot(Vm, Ighk050, label = "Ighk0.50e-3mM")
plt.plot(Vm, Ighk100, label = "Ighk1.00e-3mM")
plt.legend()
plt.show()
