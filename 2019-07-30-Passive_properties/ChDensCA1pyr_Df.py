# exec(open('Passive properties/ChDensCA1pyr_Df.py').read())

import numpy as np
import pandas as pd
import rdesigneur as rd

# Builds the model to extract length and diameter of each compartment
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

rdes = rd.rdesigneur(
    cellProto = [
        ['Passive properties/Compartments.swc','elec'],
        # ['Morphologies/n123.swc','elec'],
    ],
)
rdes.buildModel()
leng = []
diame = []
Dfindex = []
for el in moose.le('/model/elec'):
    try:
        leng.append(moose.element(el).length)
        diame.append(moose.element(el).diameter)
        Dfindex.append(moose.element(el).path.split('/')[-1])
        print(moose.element(el).path)
        print(moose.element(el).length)
        print(moose.element(el).diameter)
        print('#####')
    except:
        continue

# defines the dataframe
ChDensCA1pyr_Df = pd.DataFrame(columns=['len', 'diam', 'Em', 'RM','RA', 'CM', 'Ca_Basal', 'Ca_tau', 'Ca_B', 'Na_Gbar', 'K_DR_Gbar', 'K_A_Gbar', 'K_M_Gbar', 'h_Gbar', 'Ca_T_Gbar', 'Ca_R_Gbar', 'Ca_L_Gbar', 'Ca_N_Gbar', 'K_SK_Gbar', 'K_BK_Gbar'], index=Dfindex)

# ChDensCA1pyr_Df.len = [3.25e-05, 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.00012, 0.00011, 0.00011, 0.00011, 0.00011, 0.00011, 0.00011, 0.00011, 0.00011]
# ChDensCA1pyr_Df.diam = [3.25e-05, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 5.78e-06, 4.84e-06, 4.84e-06, 4.84e-06, 4.84e-06, 4.84e-06, 4.84e-06, 4.84e-06, 4.84e-06]
ChDensCA1pyr_Df.len = leng
ChDensCA1pyr_Df.diam = diame
noof = len(diame)
ChDensCA1pyr_Df.Em = -0.070*np.ones(noof)
ChDensCA1pyr_Df.RM = 1*np.ones(noof)
ChDensCA1pyr_Df.RA = 1.5*np.ones(noof)
ChDensCA1pyr_Df.CM = 0.0031*np.ones(noof)
ChDensCA1pyr_Df.Ca_Basal = 0.05e-3*np.ones(noof)
ChDensCA1pyr_Df.Ca_tau = 0.029*np.ones(noof)
ChDensCA1pyr_Df.Ca_B = 575792.7*np.ones(noof) # Will be setup after model building
ChDensCA1pyr_Df.Na_Gbar = 0*np.ones(noof) # [300,150,0,200,0,0,0,0,0,0,0,150,0,200,0,0,0,0,0]
ChDensCA1pyr_Df.K_DR_Gbar = 0*np.ones(noof) # [250,100,50,200,0,0,0,0,0,0,0,100,50,200,0,0,0,0,0]
ChDensCA1pyr_Df.K_A_Gbar = 0*np.ones(noof) # [50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ChDensCA1pyr_Df.K_M_Gbar = 0*np.ones(noof) # [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ChDensCA1pyr_Df.h_Gbar = 0*np.linspace(0.6,6,noof) # [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
ChDensCA1pyr_Df.Ca_T_Gbar = 0*0.06*np.ones(noof) # [10,20,10,40,20,20,20,10,10,10,0,20,10,30,20,20,10,10,0]
ChDensCA1pyr_Df.Ca_R_Gbar = 0*np.ones(noof) # [10,20,10,40,20,20,20,10,10,10,0,20,10,30,20,20,10,10,0]
ChDensCA1pyr_Df.Ca_L_Gbar = 0*np.ones(noof) # [10,20,10,40,20,20,20,10,10,10,0,20,10,30,20,20,10,10,0]
ChDensCA1pyr_Df.Ca_N_Gbar = 0*np.ones(noof) # [10,20,10,40,20,20,20,10,10,10,0,20,10,30,20,20,10,10,0]
ChDensCA1pyr_Df.K_SK_Gbar = 0*np.ones(noof) # [8,8,8,8,8,8,8,8,8,8,0,8,8,8,8,8,8,8,0]
ChDensCA1pyr_Df.K_BK_Gbar = 0*0.06*np.ones(noof) # [50,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
