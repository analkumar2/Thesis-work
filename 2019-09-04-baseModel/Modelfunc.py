#exec(open('Modelfunc.py').read())

def Modelfunc(runfor=2, stimul='Vclamp', Injectcurr=150e-12, gl = 0.25, Ca_concCaBasal='0.05e-3', Ca_conctau='0.1', Ca_L_changbar='3', Ca_N_changbar='3', Ca_T_changbar='3', h_changbar='0.25', K_A_changbar='30', K_BK_changbar='8', K_D_changbar='0.05', K_DR_changbar='3', K_M_changbar='1.1', K_SK_changbar='1', Na_changbar='1000'):
    global sm_area
    #Deleting any previous run of the model
    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        moose.delete('/model')
        moose.delete('/library/Ca_conc') #Deletes library Ca_conc so that CaBasal can be changed
    except:
        pass

    #Which stimuli to give
    if stimul == 'Vclamp':
        stimuli = ['soma', '1', '.', 'vclamp', f'{Vleveli} + (t>{preStimTime} && t<{preStimTime+injectTime}) * {Vlevelf-Vleveli}' ]
    elif stimul == 'Iclamp':
        stimuli = ['soma', '1', '.', 'inject', f'(t>={preStimTime} && t<={preStimTime+injectTime}) ? {Injectcurr} : 0']

    #Initial model run to calculate RM and Gbars so that input resistance and Vrest are as needed.
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
            ['soma', 'RM', str(1/gl), 'CM', str(CM), 'initVm', str(Vrest), 'Em', str(Vrest)],
        ],
        chanDistrib = [
            ['Ca_conc', 'soma', 'CaBasal', Ca_concCaBasal, 'tau', Ca_conctau],
            ['Ca_L_chan', 'soma', 'Gbar', Ca_L_changbar],
            ['Ca_N_chan', 'soma', 'Gbar', Ca_N_changbar],
            ['Ca_T_chan', 'soma', 'Gbar', Ca_T_changbar],
            ['h_chan', 'soma', 'Gbar', h_changbar],
            ['K_A_chan', 'soma', 'Gbar', K_A_changbar],
            ['K_BK_chan', 'soma', 'Gbar', K_BK_changbar],
            ['K_D_chan', 'soma', 'Gbar', K_D_changbar],
            ['K_DR_chan', 'soma', 'Gbar', K_DR_changbar],
            ['K_M_chan', 'soma', 'Gbar', K_M_changbar],
            ['K_SK_chan', 'soma', 'Gbar', K_SK_changbar],
            ['Na_chan', 'soma', 'Gbar', Na_changbar],
        ],
        stimList = [
            ['soma', '1', '.', 'vclamp', str(Vrest) ],
        ],
        plotList = [
            # ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
            # ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
            # ['soma', '1', ',', 'inject', 'Injected current MOOSE'],
            # ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
            # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
            # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
        ],
    )

    rdes.buildModel()

    #Changing vclamp parameters
    try:
        moose.element( '/model/elec/soma/vclamp' ).gain = CM*sm_area/elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).tau = 5*elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).ti = elecPlotDt
        moose.element( '/model/elec/soma/vclamp' ).td = 0
    except:
        pass

    #Setting Ca_conc B value
    try:
        moose.element('/model/elec/soma/Ca_conc').B = 1000e3/(2*F*depth*np.pi*sm_diam*sm_len*2)
        # moose.element('/model/elec/soma/Ca_conc').B *= 2
        # moose.element('/model/elec/soma/Ca_conc').B = 0
    except:
        pass

    moose.reinit()
    moose.start(0.8)

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
    Ca_conc = moose.element('/model/elec/soma/Ca_conc')
    soma = moose.element('/model/elec/soma')

    #Calculate Leak conductance, active conductance and input conductance
    Gl = 1/moose.element('/model/elec/soma').Rm
    Gk = Na_chan.Gk + K_SK_chan.Gk + K_M_chan.Gk + K_DR_chan.Gk + K_D_chan.Gk + K_BK_chan.Gk + K_A_chan.Gk + h_chan.Gk + Ca_T_chan.Gk + Ca_N_chan.Gk + Ca_L_chan.Gk #initial calculation of total active conductance
    Ginc = Gl+Gk

    print(f'Gin should be {Gin} but it is {Ginc}')
    print(f'gl should be between 0.15 and {Gin/sm_area} and it is {Gl/sm_area}')
    print(f'Gk should be between 0 and {Gin} and it is {Gk}')
    #First change Gk, then K_D_chanGbar and so on so that input resistance turns out to be around Gin
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
    #Set up proper Em so that Vrest is around the needed value
    Gl = 1/moose.element('/model/elec/soma').Rm
    Gk = Na_chan.Gk + K_SK_chan.Gk + K_M_chan.Gk + K_DR_chan.Gk + K_D_chan.Gk + K_BK_chan.Gk + K_A_chan.Gk + h_chan.Gk + Ca_T_chan.Gk + Ca_N_chan.Gk + Ca_L_chan.Gk #initial calculation of total active conductance
    Ginc = Gl+Gk
    Iactive = Na_chan.Ik + K_SK_chan.Ik + K_M_chan.Ik + K_DR_chan.Ik + K_D_chan.Ik + K_BK_chan.Ik + K_A_chan.Ik + h_chan.Ik + Ca_T_chan.Ik + Ca_N_chan.Ik + Ca_L_chan.Ik
    Em = -Iactive/Gl + Vrest

    print('After calculations,')
    print(f'Em is {Em}')
    print(f'Gin should be {Gin} but it is {Ginc}')
    print(f'gl should be between 0.15 and {Gin/sm_area} and it is {Gl/sm_area}')
    print(f'Gk should be between 0 and {Gin} and it is {Gk}')

    K_D_changbar=str(K_D_chan.Gbar/sm_area)
    K_M_changbar=str(K_M_chan.Gbar/sm_area)
    K_BK_changbar=str(K_BK_chan.Gbar/sm_area)
    h_changbar=str(h_chan.Gbar/sm_area)

    #Actual model run
    moose.delete('/model')
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
            ['soma', 'Rm', str(1/Gl), 'CM', str(CM), 'initVm', str(Vrest), 'Em', str(Em)],
        ],
        chanDistrib = [
            ['Ca_conc', 'soma', 'CaBasal', Ca_concCaBasal, 'tau', Ca_conctau],
            ['Ca_L_chan', 'soma', 'Gbar', Ca_L_changbar],
            ['Ca_N_chan', 'soma', 'Gbar', Ca_N_changbar],
            ['Ca_T_chan', 'soma', 'Gbar', Ca_T_changbar],
            ['h_chan', 'soma', 'Gbar', h_changbar],
            ['K_A_chan', 'soma', 'Gbar', K_A_changbar],
            ['K_BK_chan', 'soma', 'Gbar', K_BK_changbar],
            ['K_D_chan', 'soma', 'Gbar', K_D_changbar],
            ['K_DR_chan', 'soma', 'Gbar', K_DR_changbar],
            ['K_M_chan', 'soma', 'Gbar', K_M_changbar],
            ['K_SK_chan', 'soma', 'Gbar', K_SK_changbar],
            ['Na_chan', 'soma', 'Gbar', Na_changbar],
        ],
        stimList = [stimuli,
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
            ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
            ['soma', '1', '.', 'inject', 'Injected current MOOSE'],
            ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
            # ['soma', '1', 'K_M_chan', 'Ik', 'Channel current MOOSE'],
            # ['soma', '1', 'Na_chan', 'Gk', 'Channel conductance MOOSE'],
        ],
    )
    rdes.buildModel()
    #Setup clock table to record time
    clk = moose.element('/clock')
    plott = moose.Table('/model/graphs/plott')
    moose.connect(plott, 'requestOut', clk, 'getCurrentTime')

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
    Ca_conc = moose.element('/model/elec/soma/Ca_conc')
    soma = moose.element('/model/elec/soma')

    #Outputting final values of all the parameters
    Parametersout = {}
    sm_area = soma.diameter*soma.length*np.pi
    Parametersout['sm_diam'] = soma.diameter
    Parametersout['sm_len'] = soma.length
    Parametersout['RM'] = soma.Rm*sm_area
    Parametersout['CM'] = soma.Cm/(sm_area)
    Parametersout['Em'] = soma.Em
    Parametersout['Ca_concCaBasal'] = Ca_conc.CaBasal
    Parametersout['Ca_conctau'] = Ca_conc.tau
    Parametersout['Ca_concB'] = Ca_conc.B
    Parametersout['Ca_L_changbar'] = Ca_L_chan.Gbar/sm_area
    Parametersout['Ca_N_changbar'] = Ca_N_chan.Gbar/sm_area
    Parametersout['Ca_T_changbar'] = Ca_T_chan.Gbar/sm_area
    Parametersout['h_changbar'] = h_chan.Gbar/sm_area
    Parametersout['K_A_changbar'] = K_A_chan.Gbar/sm_area
    Parametersout['K_BK_changbar'] = K_BK_chan.Gbar/sm_area
    Parametersout['K_D_changbar'] = K_D_chan.Gbar/sm_area
    Parametersout['K_DR_changbar'] = K_DR_chan.Gbar/sm_area
    Parametersout['K_M_changbar'] = K_M_chan.Gbar/sm_area
    Parametersout['K_SK_changbar'] = K_SK_chan.Gbar/sm_area
    Parametersout['Na_changbar'] = Na_chan.Gbar/sm_area

    #Outputting some charateristics of the model
    characteristics = {}
    characteristics['sm_area'] = sm_area
    characteristics['sm_vol'] = sm_area
    characteristics['Input conductance'] = Ginc
    characteristics['Input resistance'] = 1/Ginc
    characteristics['Total active conductance'] = Gk
    characteristics['Resting potential Membrane'] = Vrest

    #Outputting the model outputs
    if stimul=='Vclamp':
        Vmvec=moose.element('/model/graphs/plot0').vector
        Ivec=moose.element('/model/graphs/plot1').vector
        Cavec=moose.element('/model/graphs/plot3').vector
        tvec=plott.vector
    elif stimul=='Iclamp':
        Vmvec=moose.element('/model/graphs/plot0').vector
        Ivec=moose.element('/model/graphs/plot1').vector
        Cavec=moose.element('/model/graphs/plot2').vector
        tvec=plott.vector
    return [Parametersout, characteristics, Vmvec, Ivec, Cavec, tvec]
