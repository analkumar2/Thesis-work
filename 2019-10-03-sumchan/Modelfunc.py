#exec(open('Modelfunc.py').read())


Parameters = {}
Parameters['sm_diam'] = sm_diam
Parameters['sm_len'] = sm_len
Parameters['RM'] = 2.053773343
Parameters['CM'] = CM
Parameters['Em'] = Vrest
Parameters['K_DR_changbar'] = 4.477080862
Parameters['Na_changbar'] = 51.0502324
Parameters['Na_act'] = [-138.9, -5.34, 2.60187360e-05, -7.05694169e+01, 2.64467199e-02, -1.02623323e+02, 1.73923802e-05]
Parameters['Na_inact'] = [250, 12.5, 2.004014672869906e-09, -360.599, 1.866086164497956e-09, -454.5, 0.00047 ]
Parameters['K_DR_act'] = [-114.1, 1.483, 1.41182507e-01, -7.98484941e+01, 4.40570637e+00, -1.14069277e+02, 3.84599317e-13]

def Modelfunc(runfor=2, stimul='Iclamp', Injectcurr=Injectcurr, gl=1/Parameters['RM'], K_DR_changbar=str(Parameters['K_DR_changbar']), Na_changbar=str(Parameters['Na_changbar']), Na_act=Parameters['Na_act'], Na_inact=Parameters['Na_inact'], K_DR_act=Parameters['K_DR_act']):
    sys.stdout = open(os.devnull, 'w') #Supresses any output on terminal
    global sm_area
    #Deleting any previous run of the model
    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        moose.delete('/model')
        moose.delete('/library')
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
            ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
            ['Na_Chan_(Migliore2018).Na_Chan()', 'Na_chan'],
        ],
        passiveDistrib = [
            ['soma', 'RM', str(1/gl), 'CM', str(CM), 'initVm', str(Vrest), 'Em', str(Vrest)],
        ],
        chanDistrib = [
            ['K_DR_chan', 'soma', 'Gbar', K_DR_changbar],
            ['Na_chan', 'soma', 'Gbar', Na_changbar],
        ],
        stimList = [
            ['soma', '1', '.', 'vclamp', str(Vrest) ],
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
            ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
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
    moose.start(0.4)
    # rdes.display()

    Na_chan = moose.element('/model/elec/soma/Na_chan')
    K_DR_chan = moose.element('/model/elec/soma/K_DR_chan')
    soma = moose.element('/model/elec/soma')

    #Calculate Leak conductance, active conductance and input conductance
    Gl = 1/moose.element('/model/elec/soma').Rm
    Gk = Na_chan.Gk + K_DR_chan.Gk #initial calculation of total active conductance
    Ginc = Gl+Gk

    print(f'Gin should be {Gin} but it is {Ginc}')
    print(f'gl should be between 0.15 and {Gin/sm_area} and it is {Gl/sm_area}')
    print(f'Gk should be between 0 and {Gin} and it is {Gk}')
    # First change Gk, then K_D_chanGbar and so on so that input resistance turns out to be around Gin
    if Ginc<Gin-1e-9 or Ginc>Gin+1e-9:
        if Gk<Gin:
            Gl = Gin - Gk
        elif Gk>Gin:
            sureshoterror
        moose.element('/model/elec/soma').Rm = 1/Gl

    moose.reinit()
    moose.start(0.4)
    # rdes.display()
    #Set up proper Em so that Vrest is around the needed value
    Gl = 1/moose.element('/model/elec/soma').Rm
    Gk = Na_chan.Gk + K_DR_chan.Gk #initial calculation of total active conductance
    Ginc = Gl+Gk
    Iactive = Na_chan.Ik + K_DR_chan.Ik
    Em = -Iactive/Gl + Vrest

    print('After calculations,')
    print(f'Em is {Em}')
    print(f'Gin should be {Gin} but it is {Ginc}')
    print(f'gl should be between 0.15 and {Gin/sm_area} and it is {Gl/sm_area}')
    print(f'Gk should be between 0 and {Gin} and it is {Gk}')

    #Actual model run
    moose.delete('/model')
    rdes = rd.rdesigneur(
        elecPlotDt = elecPlotDt,
        cellProto = [
            ['somaProto', 'soma', sm_diam, sm_len],
        ],
        chanProto = [
            ['K_DR_Chan_(Migliore2018).K_DR_Chan()', 'K_DR_chan'],
            ['Na_Chan_(Migliore2018).Na_Chan()', 'Na_chan'],
        ],
        passiveDistrib = [
            ['soma', 'Rm', str(1/Gl), 'CM', str(CM), 'initVm', str(Vrest), 'Em', str(Em)],
        ],
        chanDistrib = [
            ['K_DR_chan', 'soma', 'Gbar', K_DR_changbar],
            ['Na_chan', 'soma', 'Gbar', Na_changbar],
        ],
        stimList = [stimuli,
        ],
        plotList = [
            ['soma', '1', '.', 'Vm', 'Soma Membrane potential MOOSE'],
            ['soma', '1', 'vclamp', 'current', 'Soma holding current MOOSE'],
            ['soma', '1', '.', 'inject', 'Injected current MOOSE'],
            # ['soma', '1', 'Ca_conc', 'Ca', 'soma calcium conc MOOSE'],
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
    # rdes.display()

    Na_chan = moose.element('/model/elec/soma/Na_chan')
    K_DR_chan = moose.element('/model/elec/soma/K_DR_chan')
    soma = moose.element('/model/elec/soma')

    #Outputting final values of all the parameters
    Parametersout = {}
    sm_area = soma.diameter*soma.length*np.pi
    Parametersout['sm_diam'] = soma.diameter
    Parametersout['sm_len'] = soma.length
    Parametersout['RM'] = soma.Rm*sm_area
    Parametersout['CM'] = soma.Cm/(sm_area)
    Parametersout['Em'] = soma.Em
    Parametersout['K_DR_changbar'] = K_DR_chan.Gbar/sm_area
    Parametersout['Na_changbar'] = Na_chan.Gbar/sm_area
    Parametersout['Na_act'] = Na_act
    Parametersout['Na_inact'] = Na_inact
    Parametersout['K_DR_act'] = K_DR_act

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
        # Cavec=moose.element('/model/graphs/plot3').vector
        tvec=plott.vector
    elif stimul=='Iclamp':
        Vmvec=moose.element('/model/graphs/plot0').vector
        Ivec=moose.element('/model/graphs/plot1').vector
        # Cavec=moose.element('/model/graphs/plot2').vector
        tvec=plott.vector

    sys.stdout = sys.__stdout__
    return [Parametersout, characteristics, Vmvec, Ivec, tvec]
