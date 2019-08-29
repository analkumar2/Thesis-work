# exec(open('baseModel.py').read())

import moose
import numpy as np
import matplotlib.pyplot as plt
import rdesigneur as rd

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
    # moose.delete('/library/Ca_conc')
except:
    pass

F = 96485.3329
sm_diam = 60e-6 # so that Cm is 113pF
sm_len = 60e-6 # so that Cm is 113pF
Gin = 7.35e-9
P_K_D_chan_m70 = 98.6e-2
P_h_chan_m70 = 20.2e-2
P_K_M_chan_m70 = 4.74e-2
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
gl = 0.25
mingl = 0.15
maxgl = Gin/sm_area
CM = 0.01
Em = -0.070
Vrest = -0.070
depth = 0.1 # No units. For ca_conc B
elecPlotDt = 0.00005
preStimTime = 1
injectTime = 0.5
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
Injectcurr = 150e-12

def makeModel(runfor=0.8, stimul='Vclamp', RM='4', Em='-0.07', K_M_changbar='1.1', K_D_changbar='0.05', h_changbar='0.25', K_BK_changbar='8'):
    if stimul == 'Vclamp':
        stimuli = ['soma', '1', '.', 'vclamp', f'-0.070 + (t>{preStimTime} && t<{preStimTime+injectTime}) * 0.07' ]
    elif stimul == 'Iclamp':
        stimuli = ['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {Injectcurr} : 0']

    rdes = rd.rdesigneur(
        elecPlotDt = elecPlotDt,
        cellProto = [
            ['somaProto', 'soma', sm_diam, sm_len],
        ],
        chanProto = [
            ['Ca_Conc_(Common).Ca_Conc()', 'Ca_conc'],
            ['Ca_L_Chan_(Migliore2018).Ca_L_Chan()', 'Ca_L_chan'],
            ['Ca_N_Chan_(Migliore2018).Ca_N_Chan()', 'Ca_N_chan'],
            ['Ca_T_Chan_(Migliore2018).Ca_T_Chan()', 'Ca_T_chan'],
            ['h_Chan_(Migliore2018).h_Chan()', 'h_chan'],
            ['K_A_Chan_(Migliore2018).K_A_Chan()', 'K_A_chan'],
            ['K_BK_Chan_(Migliore2018).K_BK_Chan()', 'K_BK_chan'],
            ['K_D_Chan_(Migliore2018).K_D_Chan()', 'K_D_chan'],
            ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
            ['K_M_Chan_(Migliore2018).K_M_Chan()', 'K_M_chan'],
            ['K_SK_Chan_(Migliore2018).K_SK_Chan()', 'K_SK_chan'],
            ['Na_Chan_(Migliore2018).Na_Chan()', 'Na_chan'],
        ],
        passiveDistrib = [
            ['soma', 'RM', RM, 'CM', str(CM), 'initVm', str(Vrest), 'Em', Em],
        ],
        chanDistrib = [
            ['Ca_conc', 'soma', 'CaBasal', str(0.05e-3)],
            ['Ca_L_chan', 'soma', 'Gbar', str(3)],
            ['Ca_N_chan', 'soma', 'Gbar', str(3)],
            ['Ca_T_chan', 'soma', 'Gbar', str(3)],
            ['h_chan', 'soma', 'Gbar', h_changbar],
            ['K_A_chan', 'soma', 'Gbar', str(30)],
            ['K_BK_chan', 'soma', 'Gbar', K_BK_changbar],
            ['K_D_chan', 'soma', 'Gbar', K_D_changbar],
            ['K_DR_chan', 'soma', 'Gbar', str(3)],
            ['K_M_chan', 'soma', 'Gbar', K_M_changbar],
            ['K_SK_chan', 'soma', 'Gbar', str(1)],
            ['Na_chan', 'soma', 'Gbar', str(1000)],
        ],
        stimList = [stimuli
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
            ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
            # ['soma', '1', ',', 'inject', 'Injected current MOOSE'],
            ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
            # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
            # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
        ],
    )

    rdes.buildModel()

    try:
        moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).td = 0
    except:
        pass

    try:
        moose.element('/model/elec/soma/Ca_conc').B = 1000e3/(2*F*depth*np.pi*sm_diam*sm_len*2)
        # moose.element('/model/elec/soma/Ca_conc').B *= 2
        # moose.element('/model/elec/soma/Ca_conc').B = 0
    except:
        pass

    moose.reinit()
    moose.start(runfor)
    return rdes

rdes = makeModel(runfor = 0.8, stimul = 'Vclamp')
Na_chan = moose.element('/model/elec/soma/Na_chan')
K_SK_chan = moose.element('/model/elec/soma/K_SK_chan')
K_M_chan = moose.element('/model/elec/soma/K_M_chan')
K_DR_chan = moose.element('/model/elec/soma/K_DR_chan')
K_D_chan = moose.element('/model/elec/soma/K_D_chan')
K_BK_chan = moose.element('/model/elec/soma/K_BK_chan')
K_A_chan = moose.element('/model/elec/soma/K_A_chan')
h_chan = moose.element('/model/elec/soma/h_chan')
Ca_T_chan = moose.element('/model/elec/soma/Ca_T_chan')
Ca_N_chan = moose.element('/model/elec/soma/Ca_N_chan')
Ca_L_chan = moose.element('/model/elec/soma/Ca_L_chan')

Gl = 1/moose.element('/model/elec/soma').Rm
Gk = Na_chan.Gk + K_SK_chan.Gk + K_M_chan.Gk + K_DR_chan.Gk + K_D_chan.Gk + K_BK_chan.Gk + K_A_chan.Gk + h_chan.Gk + Ca_T_chan.Gk + Ca_N_chan.Gk + Ca_L_chan.Gk #initial calculation of total active conductance
Ginc = Gl+Gk

print(f'Gin should be {Gin} but it is {Ginc}')
print(f'gl should be between 0.15 and {Gin/sm_area} and it is {Gl/sm_area}')
print(f'Gk should be between 0 and {Gin} and it is {Gk}')
if Ginc<Gin-1e-9 or Ginc>Gin+1e-9:
    if Gk<Gin:
        Gl = Gin - Gk
    elif Gk>Gin:
        Gl = mingl*sm_area
        K_D_chanGbar = (Gin - Gl - (Gk - K_D_chan.Gk))/P_K_D_chan_m70
        if K_D_chanGbar<0:
            K_D_chanGbar = 0
            Gk = -1*(0 - Na_chan.Gk - K_SK_chan.Gk - K_M_chan.Gk - K_DR_chan.Gk - K_BK_chan.Gk - K_A_chan.Gk - h_chan.Gk - Ca_T_chan.Gk - Ca_N_chan.Gk - Ca_L_chan.Gk)
            h_chanGbar = (Gin - Gl - (Gk - h_chan.Gk))/P_h_chan_m70
            if h_chanGbar<0:
                h_chanGbar = 0
                Gk = -1*(0 - Na_chan.Gk - K_SK_chan.Gk - K_M_chan.Gk - K_DR_chan.Gk - K_BK_chan.Gk - K_A_chan.Gk - Ca_T_chan.Gk - Ca_N_chan.Gk - Ca_L_chan.Gk)
                K_M_chanGbar = (Gin - Gl - (Gk - K_M_chan.Gk))/P_K_M_chan_m70
                if K_M_chanGbar<0:
                    print('Recheck model. BK currents are too high')
                    plt.close()
                    sureshoterror
                K_M_chan.Gbar = K_M_chanGbar
            h_chan.Gbar = h_chanGbar
        K_D_chan.Gbar = K_D_chanGbar
    moose.element('/model/elec/soma').Rm = 1/Gl

moose.reinit()
moose.start(0.8)
# # makeModel(runfor)
Gl = 1/moose.element('/model/elec/soma').Rm
Gk = Na_chan.Gk + K_SK_chan.Gk + K_M_chan.Gk + K_DR_chan.Gk + K_D_chan.Gk + K_BK_chan.Gk + K_A_chan.Gk + h_chan.Gk + Ca_T_chan.Gk + Ca_N_chan.Gk + Ca_L_chan.Gk #initial calculation of total active conductance
Ginc = Gl+Gk
Iactive = Na_chan.Ik + K_SK_chan.Ik + K_M_chan.Ik + K_DR_chan.Ik + K_D_chan.Ik + K_BK_chan.Ik + K_A_chan.Ik + h_chan.Ik + Ca_T_chan.Ik + Ca_N_chan.Ik + Ca_L_chan.Ik
Em = -Iactive/Gl + Vrest
moose.element('/model/elec/soma').Em = Em

print('After calculations,')
print(f'Em is {Em}')
print(f'Gin should be {Gin} but it is {Ginc}')
print(f'gl should be between 0.15 and {Gin/sm_area} and it is {Gl/sm_area}')
print(f'Gk should be between 0 and {Gin} and it is {Gk}')

K_D_changbar=K_D_chan.Gbar/sm_area
K_M_changbar=K_M_chan.Gbar/sm_area
K_BK_changbar=K_BK_chan.Gbar/sm_area
h_changbar=h_chan.Gbar/sm_area

moose.delete('/model')
rdes = makeModel(runfor=0.8, stimul='Iclamp', RM=str(sm_area/Gl), Em=str(Em), K_D_changbar=str(K_D_changbar), K_M_changbar=str(K_M_changbar), K_BK_changbar=str(K_BK_changbar), h_changbar=str(h_changbar))
moose.reinit()
moose.start(2)
rdes.display()
