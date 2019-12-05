#exec(open('calcinftau.py').read())

import numpy as np
import matplotlib.pyplot as plt
plt.cla()
ax = plt.subplot(111)
inf = []
vfinallist = np.arange(-0.120,0.110,0.010)
maxI_list = []
for vfinal in vfinallist:
    exec(open('testmodel.py').read())
    _I = max(Ivec[(np.abs(np.array(tvec)-0.9)).argmin():(np.abs(np.array(tvec)-1)).argmin()])
    maxI_list.append(_I)

Gbar_leak_c = (maxI_list[1]-maxI_list[0])/(vfinallist[1]-vfinallist[0])
Em_c = -(maxI_list[1]*vfinallist[0] - maxI_list[0]*vfinallist[1])/(maxI_list[0]-maxI_list[1])
Gbar_active_c = ((maxI_list[-1]-Gbar_leak_c*(vfinallist[-1]-Em_c))-(maxI_list[-2]-Gbar_leak_c*(vfinallist[-2]-Em_c)))/(vfinallist[-1]-vfinallist[-2]) #Calculates by subtracting leak current
Ek_c = -((maxI_list[-1]-Gbar_leak_c*(vfinallist[-1]-Em_c))*vfinallist[-2] - (maxI_list[-2]-Gbar_leak_c*(vfinallist[-2]-Em_c))*vfinallist[-1])/((maxI_list[-2]-Gbar_leak_c*(vfinallist[-2]-Em_c))-(maxI_list[-1]-Gbar_leak_c*(vfinallist[-1]-Em_c)))

G_active_c = (maxI_list - Gbar_leak_c*(vfinallist-Em_c))/(vfinallist-Ek_c)
inf = G_active_c/max(G_active_c)

plt.title('Voltage clamp just KM=50 and KDR=50')
plt.ylabel('Holding current (A)')
plt.xlabel('time (s)')
plt.legend()
pickle.dump(ax, open('justKMDR_Vclamp.pickle', 'wb'))
plt.show()

# Vmin = -0.100
# Vmax = 0.100
# Vdivs = 3000
# v = np.linspace(Vmin,Vmax, Vdivs)
# plt.plot(vfinallist, inf, label='KDR=50, KM=5')
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt
# import pickle
# ax = pickle.load(open('justKMDR_Vclamp.pickle', 'rb'))
# plt.show()
