import numpy as np
import matplotlib.pyplot as plt

# plt.plot(np.arange(-8,11,1), [0,0,0,0,0,200,0,150,300,150,0,200,0,0,0,0,0,0,0], label = 'Na Channel')
# plt.plot(np.arange(-8,11,1), [0,50,50,120,120,120,50,80,40,80,50,170,170,170,100,100,50,50,0], label = 'Ca Channel')
# plt.plot(np.arange(-8,11,1), [0,0,0,0,0,200,0,50,150,50,0,200,0,0,0,0,0,0,0], label = 'K(DR) Channel')
# plt.plot(np.arange(-8,11,1), [0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,0], label = 'K(AHP) Channel')
# plt.plot(np.arange(-8,11,1), [0,50,50,100,100,100,50,200,100,200,50,150,150,150,150,150,50,50,0], label = 'K(C) Channel')
# plt.plot(np.arange(-8,11,1), [0,0,0,0,0,0,0,0,50,0,0,0,0,0,0,0,0,0,0], label = 'K(A) Channel')

plt.plot(np.arange(-8,11,1), [0,0,0,0,0,200,0,150,300,150,0,200,0,0,0,0,0,0,0], label = 'Na Channel')
plt.plot(np.arange(-8,11,1), [0,50,50,70,70,120,50,80,40,80,50,170,70,70,70,50,50,50,0], label = 'Ca Channel')
plt.plot(np.arange(-8,11,1), [0,0,0,0,0,200,50,100,250,100,50,200,0,0,0,0,0,0,0], label = 'K(DR) Channel')
plt.plot(np.arange(-8,11,1), [0,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,0], label = 'K(AHP) Channel')
plt.plot(np.arange(-8,11,1), [0,50,50,50,50,100,50,200,100,200,50,150,50,50,50,50,50,50,0], label = 'K(C) Channel')
plt.plot(np.arange(-8,11,1), [0,0,0,0,0,0,0,0,50,0,0,0,0,0,0,0,0,0,0], label = 'K(A) Channel')

plt.xlabel('Compartment number')
plt.ylabel('Conductance (S*m^-2)')
plt.xticks(np.arange(-8, 10, step=1))
plt.title('Distribution of channel densities across compartments in CA1 pyramidal neurons')
plt.legend()
plt.show()