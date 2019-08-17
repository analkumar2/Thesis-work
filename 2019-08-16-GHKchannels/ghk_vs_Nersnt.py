# The ghk equatioin is taken from migliore2018 cat.mod where gcat = P*z**2*F**2*[S]o/(R*T)
import numpy as np
import matplotlib.pyplot as plt


F = 96485
R = 8.314
T = 306
FRT = F/R/T

cao = 2
cai = np.array([0.05e-3,1e-3]) #ci is usually in the range 0.05e-3mM to 1e-3. But it does not matter much
caz = 2
caE = R*T/(caz*F)*np.log(cao/cai)
nao = 141.5
nai = 4.3
naz = 1
naE = 0.092
ko = 3.3
ki = 140
kz = 1
kE = -0.100

Vm = np.linspace(-0.100,0.100,1000)

nai_by_g_ghk = Vm*(nai/nao*np.exp(Vm*naz*FRT)-1)/(np.exp(Vm*naz*FRT)-1)
nai_by_g_nernst = Vm - naE
ki_by_g_ghk = Vm*(ki/ko*np.exp(Vm*kz*FRT)-1)/(np.exp(Vm*kz*FRT)-1)
ki_by_g_nernst = Vm - kE
cai_by_g_ghk = Vm*(cai[0]/cao*np.exp(Vm*caz*FRT)-1)/(np.exp(Vm*caz*FRT)-1)
cai_by_g_nernst = Vm - caE[0]

# plt.plot(Vm, nai_by_g_ghk, label='Na ghk')
# plt.plot(Vm, ki_by_g_ghk, label='K ghk')
# plt.plot(Vm, cai_by_g_ghk, label='Ca ghk')
# plt.plot(Vm, nai_by_g_nernst, label='Na nernst')
# plt.plot(Vm, ki_by_g_nernst, label='K nernst')
# plt.plot(Vm, cai_by_g_nernst, label='Ca nernst')

# plt.plot(Vm, nai_by_g_ghk/nai_by_g_nernst, label='Na')
# plt.plot(Vm, ki_by_g_ghk/ki_by_g_nernst, label='K')
plt.plot(Vm, cai_by_g_ghk/cai_by_g_nernst, label='Ca')
plt.xlabel('Membrane potential (V)')
# plt.ylabel('Driving force (V)')
plt.ylabel('Driving force ratio (1)')
plt.legend()
# plt.title('Comparing Na, K, and Ca currents/g for ghk and nernst implementations')
plt.title('Driving force ratios')
plt.show()
