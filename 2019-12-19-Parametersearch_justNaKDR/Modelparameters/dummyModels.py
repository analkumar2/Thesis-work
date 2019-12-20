# exec(open('Modelparameters/dummyModels.py').read())

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, E, F):
    # alge model
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    yl = (v-A)/-B
    yr = (v-A)/E
    Tau = (C + (1 + yl/(np.sqrt(1+yl**2)))/2) * (D + (1 + yr/(np.sqrt(1+yr**2)))/2) * F
    Tau[Tau<0.00002] = 0.00002
    return [Inf,Tau]

Models = {}

Models['Model1'] = {'Error': 5.907811441753849, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.437332664019589e-10, 'Rm': 714010979.6046907, 'Em': -0.026818968916777236}, 'Channels': {'Na_Chan': {'Gbar': 1.1181775516241626e-05, 'Erev':0.092, 'Kinetics': '../../Compilations/Kinetics/Na_Chan_(Migliore2018)', 'gateX':[-0.038,0.0072,-0.0365,0.020,0.0161,0.0547,0.0311,0.00064], 'gateY':[-0.050,-0.004, -0.0456,0.00433,0.01198,0.0262,0.00854,0.039]}, 'K_DR_Chan': {'Gbar': 3.1760459929263415e-07, 'Erev':-0.100, 'Kinetics': '../../Compilations/Kinetics/K_DR_Chan_(Migliore2018)', 'gateX':[0.013,0.0088, 0.0125,0.0173,0,0,0.0341,0.1022]}}}}
