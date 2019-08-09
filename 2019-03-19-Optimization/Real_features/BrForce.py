#exec(open('Optimization/Custom/Real_features/BrForce.py').read())
# Only call using BrForse_inloop.py. WOn't work standalone

#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass


# Em = -0.075
# RM = 3.23
# CM = 0.0083
# Ca_tau = 0.029
# Ca_B = 1/(sm_vol*F*2)
# Na_Gbar = 76
# K_DR_Gbar = 43.4
# K_A_Gbar = 20.7
# K_M_Gbar = 3.88294e-02
# h_Gbar = 0.11
# Ca_T_Gbar = 1.57
# Ca_R_Gbar = 0.07
# Ca_L_Gbar = 0.51
# Ca_N_Gbar = 1.76
# K_SK_Gbar = 4.28006e-02
# K_BK_Gbar = 5.08559e-03

Em = rnd.uniform(-0.080, -0.060) #-0.075
RM = rnd.uniform(0.1,35) #2
CM = rnd.uniform(0.0008,0.1) #0.01
Ca_tau = rnd.uniform(0.001,0.050)
Ca_B = rnd.uniform(1/(sm_vol*F*2)*0.1,1/(sm_vol*F*2)*10)
Na_Gbar = rnd.uniform(5, 1000) #350
K_DR_Gbar = rnd.uniform(4,500) #150
K_A_Gbar = rnd.uniform(2, 250) #5
K_M_Gbar = rnd.uniform(0.003, 40) #10
h_Gbar = rnd.uniform(0.01, 1.5) #0.18
Ca_T_Gbar = rnd.uniform(0.1, 16) #0.5
Ca_R_Gbar = rnd.uniform(0.007, 0.7) #1
Ca_L_Gbar = rnd.uniform(0.05, 5.1) #5
Ca_N_Gbar = rnd.uniform(0.1, 20)
K_SK_Gbar = rnd.uniform(4e-3, 5e-1) #150
K_BK_Gbar = rnd.uniform(5e-4, 5e-2) #2475

rdes = makeModel()
rdes.buildModel()
try:
    moose.element( '/model/elec/soma/vclamp' ).gain *= 0.1
except:
    pass
moose.element('/model/elec/soma/Ca_conc').B = Ca_B
moose.element('/model/elec/soma/Ca_conc').tau = Ca_tau
moose.reinit()
moose.start(runtime)

Vtrace = moose.element('/model/graphs/plot0' ).vector

try:
    features_df = features(Vtrace*1e3,stim_start,stim_end,Inputcurr)
except:
    pass

try:
    tcost = np.nansum(features_df['cost'])
    # print(f"Total cost = {tcost}")
    print('###############################################')
except:
    tcost = 1000
    # print(f'Sorry, No spike')
    print('###############################################')
