# exec(open('inftau_fit.py').read())

import numpy as np
import brute_curvefit as bf
import matplotlib.pyplot as plt

# def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, F):
# 	# Used by moose. Is bullshit
#     Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
#     Tau = B*(v+A/B)/(C+np.exp((v+D)/F))
#     Tau[Tau<0.00002] = 0.00002
#     return [Inf,Tau]

# def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, F):
# 	# Previously selected model. Is bullshit
#     Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
#     Tau = (A*np.exp(B*v))/(1+C*np.exp(D*v)) + F
#     Tau[Tau<0.00002] = 0.00002
#     return [Inf,Tau]

# def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
# 	# logistic model
#     Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
#     Tau = (C + 1/(1+np.exp((v-A)/-B))) * (D + 1/(1+np.exp((v-A)/E))) * F
#     Tau[Tau<0.00002] = 0.00002
#     return [Inf,Tau]

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
	# alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

# def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
# 	# arctan model
#     Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
#     yl = (v-A)/-B
#     yr = (v-A)/E
#     Tau = (C + (np.pi/2+np.arctan(yl))/np.pi) * (D + (np.pi/2+np.arctan(yr))/np.pi) * F
#     Tau[Tau<0.00002] = 0.00002
#     return [Inf,Tau]

# def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
# 	# generalized logistic model. not in working condition
#     Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
#     yl = (v-A)/-B
#     yr = (v-A)/E
#     Tau = (C + (1 + np.tanh(yl))/2) * (D + (1 + np.tanh(yr))/2) * F
#     Tau[Tau<0.00002] = 0.00002
#     return [Inf,Tau]

# def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
# 	# tanh model
#     Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
#     yl = (v-A)/-B
#     yr = (v-A)/E
#     Tau = (C + (1 + np.tanh(yl))/2) * (D + (1 + np.tanh(yr))/2) * F
#     Tau[Tau<0.00002] = 0.00002
#     return [Inf,Tau]


Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
v = np.linspace(Vmin,Vmax, Vdivs)

exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore2018).py').read())
moose.Neutral('/library')
Na_Chan('Na_Chan')
Mig2018hinf = moose.element('/library/Na_Chan/gateY').tableA/moose.element('/library/Na_Chan/gateY').tableB
Mig2018htau = 1/moose.element('/library/Na_Chan/gateY').tableB

def wr_ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
    return ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F)[1]


params, error = bf.brute_scifit(wr_ChanGate, v, Mig2018htau, restrict=[[-0.1,-0.1, -0.1,0,0,0,0,0],[0.1,0.1, 0.1,0.1,0.1,0.1,0.1,1]], maxfev=1000, ntol=10000, returnnfactor=0.002)

plt.plot(v, Mig2018htau, label='ori')
plt.plot(v, wr_ChanGate(v, *params), label='fitted')
plt.legend()
plt.show()


# params, error = bf.brute_scimin(wr_ChanGate, v, Mig2018hinf, restrict=[[-1000.,-1000.,-0.1,0,0,0,0,0],[1000.,1000.,0.1,0.01,0.012,0.012,0.01,1]], ntol=10000, returnnfactor=0.002, method="Nelder-Mead", jac=None)

# plt.plot(v, Mig2018hinf, label='ori')
# plt.plot(v, wr_ChanGate(v, *params), label='fitted')
# plt.legend()
# plt.show()
### [:int(len(v)/2)]
