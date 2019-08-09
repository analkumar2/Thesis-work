#exec(open('Variablity/representativeinj.py').read())

import numpy as np
import matplotlib.pyplot as plt

samp_freq = 10000
t = np.linspace(0,3,3*samp_freq)
np.random.seed(123)
noise = [150e-12+k for k in np.random.normal(0,4.3e-12,1*samp_freq)]
i = np.concatenate([np.zeros(1*samp_freq),noise,np.zeros(1*samp_freq)])

plt.plot(t,i)
plt.title('Noisy current of 4.3pA std on top of 150pA current injection')
plt.xlabel('Time (s)')
plt.ylabel('Injected current (A)')
plt.show()
