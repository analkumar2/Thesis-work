#exec(open('Tutorial.py').read())

from brian2 import *
from brian2modelfitting import *

area=113e-6
Cm=1*ufarad*cm**-2 * area
El=-65*mV
EK=-90*mV
ENa=50*mV
VT=-63*mV

model = '''
dv/dt = (gl*(El-v) - g_na*(m*m*m)*h*(v-ENa) - g_kd*(n*n*n*n)*(v-EK) + I)/Cm : volt
dm/dt = 0.32*(mV**-1)*(13.*mV-v+VT)/
  (exp((13.*mV-v+VT)/(4.*mV))-1.)/ms*(1-m)-0.28*(mV**-1)*(v-VT-40.*mV)/
  (exp((v-VT-40.*mV)/(5.*mV))-1.)/ms*m : 1
dn/dt = 0.032*(mV**-1)*(15.*mV-v+VT)/
  (exp((15.*mV-v+VT)/(5.*mV))-1.)/ms*(1.-n)-.5*exp((10.*mV-v+VT)/(40.*mV))/ms*n : 1
dh/dt = 0.128*exp((17.*mV-v+VT)/(18.*mV))/ms*(1.-h)-4./(1+exp((40.*mV-v+VT)/(5.*mV)))/ms*h : 1
g_na : siemens (constant)
g_kd : siemens (constant)
gl   : siemens (constant)
'''

opt = NevergradOptimizer()
metric = MSEMetric()

fitter = TraceFitter(model=model,
                     input_var='I',
                     output_var='v',
                     input=inp_trace * amp,
                     output=out_trace*mV,
                     dt=0.01*ms,
                     n_samples=10,
                     method='exponential_euler',
                     param_init={'v': -65*mV})

res, error = fitter.fit(n_rounds=2,
                        optimizer=opt,
                        metric=metric,
                        gl=[2*psiemens, 200*nsiemens],
                        g_na=[200*nsiemens, 0.4*msiemens],
                        g_kd=[200*nsiemens, 200*usiemens])

traces = fitter.generate_traces()
