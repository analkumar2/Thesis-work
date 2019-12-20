# exec(open('faultycode.py').read())

Models = {}

Models['Model1'] = {'Error': 4.852056981763452, 'parameters': {'notes': '', 'Morphology': {'sm_len': 6.73077545020806e-05, 'sm_diam': 6.73077545020806e-05}, 'Passive': {'Cm': 1.4844902082245464e-10, 'Rm': 697375561.8482424, 'Em': -0.07540097023719322}, 'Channels': {'K_BK_Chan': {'Gbar': 8.130519989635794e-09, 'Kinetics': 'Kinetics/K_BK_Chan_(Migliore2018)', 'Erev': -0.09659549184061027}, 'Ca_T_Chan': {'Gbar': 1.799390010015399e-09, 'Kinetics': 'Kinetics/Ca_T_Chan_(Migliore2018)', 'Erev': 0.13521558061079192}, 'Ca_L_Chan': {'Gbar': 1.263189399144424e-08, 'Kinetics': 'Kinetics/Ca_L_Chan_(Migliore2018)', 'Erev': 0.13521558061079192}, 'Ca_N_Chan': {'Gbar': 8.948093144967204e-09, 'Kinetics': 'Kinetics/Ca_N_Chan_(Migliore2018)', 'Erev': 0.13521558061079192}, 'Na_Chan': {'Gbar': 8.95583549903416e-07, 'Kinetics': 'Kinetics/Na_Chan_(Migliore2018)', 'Erev': 0.07179644483699847}, 'Na_P_Chan': {'Gbar': 7.968172196051026e-09, 'Kinetics': 'Kinetics/Na_P_Chan_(Migliore2018)', 'Erev': 0.07179644483699847}, 'K_DR_Chan': {'Gbar': 1.9998458010877026e-06, 'Kinetics': 'Kinetics/K_DR_Chan_(Migliore2018)', 'Erev': -0.09659549184061027}, 'K_D_Chan': {'Gbar': 4.036044928042262e-09, 'Kinetics': 'Kinetics/K_D_Chan_(Migliore2018)', 'Erev': -0.09659549184061027}, 'K_A_Chan': {'Gbar': 1.3167869978551778e-05, 'Kinetics': 'Kinetics/K_A_Chan_(Migliore2018)_ghk', 'Erev': -0.09659549184061027}, 'K_M_Chan': {'Gbar': 1.3433272421527097e-08, 'Kinetics': 'Kinetics/K_M_Chan_(Migliore2018)', 'Erev': -0.09659549184061027}, 'K_SK_Chan': {'Gbar': 5.141743313994309e-10, 'Kinetics': 'Kinetics/K_SK_Chan_(Migliore2018)', 'Erev': -0.09659549184061027}, 'h_Chan': {'Gbar': 2.2461660603266964e-08, 'Kinetics': 'Kinetics/h_Chan_(Migliore2018)', 'Erev': -0.039230165502240476}}, 'Ca_Conc': {'Ca_B': 31878987667.588173, 'Ca_tau': 0.06528449318193238, 'Ca_base': 0.000619637709405835, 'Kinetics': 'Kinetics/Ca_Conc_(Common)'}}}

import moose
import MOOSEModel_2

try:
    moose.delete('/library')
except:
    pass

MOOSEModel_2.plotModel(Models['Model1'], 150e-12)

moose.delete('/library')

MOOSEModel_2.plotModel(Models['Model1'], 150e-12)

moose.delete('/library/Na_Chan')

MOOSEModel_2.plotModel(Models['Model1'], 150e-12)
