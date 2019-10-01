#exec(open('Fitting.py').read())

# Using scipy.optimize import curve_fit to find parameter values
from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt


# inf is of the form 1/(1+np.exp(A*v+B))
# tau for K_DR_chanis of the form C*np.exp(D*v)/(1+np.exp(A*v+B))
# tau for Na_chan is of the form 1e-3/(C*(v+D)/(1-np.exp(-E*v-D*E)) + F*(v+D)/(1-np.exp(E*v+D*E)))
# Hence each gate's Kinetics is defined by 6 parameters

# For K_DR_chan activatioin gate, A=-114.1, B=1.483, C=32, D=34.25, E=114.1 ,F=-1.483
# For Na_chan activatioin gate, A=-138.9, B=-5.34, C=696, D=0.03, E=139 ,F=-108
# For Na_chan inactivation gate,

#Now finding A,B,C,D for Na gates
exec(open('ChannelSlider.py').read())
plt.show()
Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
v = np.linspace(Vmin,Vmax, Vdivs)

def inf(x,A,B):
    return 1/(1+np.exp(A*x+B))

def tau(x,C,D,E,F,G):
    return 1e-3*C*np.exp(D*x)/(1+np.exp(E*x+F)) #For K_DR act
    return 1e-3/(C*(v+D)/(1-np.exp(-E*v-D*E)) + F*(v+D)/(1-np.exp(E*v+D*E))) #For Na act
    return 1e3/(C*x*np.exp(D*x)+E*np.exp(F*x))
    return (C*np.exp(D*x))/(1+E*np.exp(F*x)) + G
    return C+D*x+E*x**2+F*x**3+G*x**4
    return 1/((I*x+J)/(1+C*np.exp(D*x)) + (E*x+F)/(1+G*np.exp(H*x)))
    return 1/((I*x+J)/(1+C*np.exp(x)) + (E*x+F)/(1+np.exp(H*x)))

def Jfind(x,J):
    return tau(x,*bestCDEFGHIJ)+J


# inti_vals_ABCDEF = [-114,1.483,141.18e-3,-79.85,-114,1.483]
# inti_vals_ABCDEF = [-114,1.483,0.4e3,0.03,30/7.2,1.483]
# inti_vals_ABCDEFGHIJ = [-138.9,-5.34, 5.41530968e-05, -4.38728649e+01,  2.82998525e-02, -8.10066283e+01, 1.00000000e+00,  5.00000000e-4,  1.00000000e+00,  1.00000000e+00]
inti_vals_ABCDEFGHIJ = [-138.9,-5.34, 5.41530968e-05, -4.38728649e+01,  2.82998525e-02, -1.10066283e+02, 5.00000000e-5,  5.00000000e-4,  1.00000000e+00,  1.00000000e+00]
xgate = moose.element('/library/Na_chan/gateY')
# xgate = moose.element('/library/K_DR_chan/gateX')
bestAB, covar = curve_fit(inf,v,xgate.tableA/xgate.tableB,p0=inti_vals_ABCDEFGHIJ[:2])
# bestCDEFGHIJ, covar = curve_fit(tau,v,1/xgate.tableB,p0=inti_vals_ABCDEFGHIJ[2:7], maxfev=50000, ftol=1e-10, xtol=1e-10)
# bestJ, covar = curve_fit(Jfind,v,1/xgate.tableB,p0=inti_vals_ABCDEFGHIJ[-1], maxfev=50000, ftol=1e-10, xtol=1e-10)

# bestCDEFGHIJ_list = []
# for i inrange(1000):
#     bestCDEFGHIJ_list.append([1e5*np.random.random(1),40*np.random.random(1)-20,400*np.random.random(1)-200,400*np.random.random(1)-200, 400*np.random.random(1)-200,400*np.random.random(1)-200,400*np.random.random(1)-200,400*np.random.random(1)-200])


plt.plot(v,xgate.tableA/xgate.tableB,label='oriinf')
plt.plot(v,inf(v,*bestAB),label='custominf')
plt.legend()
plt.show()

# plt.plot(v,1/xgate.tableB,label='oritau')
# plt.plot(v,tau(v,*bestCDEFGHIJ),label='customtau')
# plt.legend()
# plt.show()


from lmfit import Model
gmodel = Model(tau)
params = gmodel.make_params()
[C,D,E,F,G] = [ 9.05675818e-09, -3.30670792e+02,  1.91046014e-08, -4.07804573e+02, 5.0e-4]
params = gmodel.make_params(C=C, D=D, E=E, F=F, G=G)
result = gmodel.fit(1/xgate.tableB, params, x=v)
print(result.values)

plt.plot(v,1/xgate.tableB,label='oritau')
plt.plot(v,result.best_fit,label='customtau')
plt.legend()
plt.show()
