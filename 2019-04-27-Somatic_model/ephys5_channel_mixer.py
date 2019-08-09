#exec(open('Somatic model/ephys5_channel_mixer.py').read())
########################################################################
# This example demonstrates the behaviour of various voltage and calcium-
# gated channels.
# Copyright (C) Upinder S. Bhalla NCBS 2018
# Released under the terms of the GNU Public License V3.
########################################################################
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np
import warnings
import moose
import rdesigneur as rd
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor


########################################
#Deleting any previous run of the model
try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

#Myown
ChP = 'Somatic model/ChannelProtos_Sri2015_base'
exec(open('Somatic model/feature_dict.py').read()) #Get the feature_range
F = 96485.3329
sm_diam = 100e-6
sm_vol = np.pi/4*sm_diam**3
sm_area = np.pi*sm_diam**2
Ca_tau = 0.029
Ca_B = 1/(sm_vol*F*2)
Na_Gbar = 76
K_DR_Gbar = 43.4
K_A_Gbar = 20.7
K_M_Gbar = 3.88294e-02
h_Gbar =  0.11
Ca_T_Gbar = 1.57
Ca_R_Gbar = 0.07
Ca_L_Gbar = 0.51
Ca_N_Gbar = 1.76
K_SK_Gbar = 4.28006e-02
K_BK_Gbar = 5.08559e-03

#######################################
lines = []
tplot = []
axes = []
sliders = []
fields = []

RM = 3.23
RA = 1.0
CM = 0.0083
diameter = sm_diam
runtime = 0.05
elecPlotDt = 0.00005
sliderMode = "Gbar"
currentChanPath = '/model/elec/soma/Na'
currentChanMod = 1.0
preStimTime = 1
injectTime = 0.5
postStimTime = 0.5
runtime = preStimTime + injectTime + postStimTime
#inject = 20e-12
helpText = """
Channel Mixer demo help.\n
This demo illustrates how different classes  of ion channels affect neuronal spiking behaviour. The model consists of a single compartment spherical soma, with Na and K_DR set at a default level to get spiking.  The top graph shows the membrane potential, and the second graph shows the calcium concentration at the soma. The current injection begins at 50 ms, and lasts for 500 ms.
In each graph the solid blue plot indicates the current calculated response of the model to the specified channel modulation and current injection.  The red dotted plot indicates what would have happened had the modulated channel not been there at all.

The default channel density for all channels is 100 Mho/m^2, except K_AHP where the default is 10 Mho/m^2. There are sliders for modulating the levels of each of the channels.  The slider lets you modulate this channel density in a range from zero to 10x the default density.
There are two terms to modulate for the Calcium dynamics: the Ca_tau is the time-course for the pump to expel calcium influx, and the Ca_thickness is the thickness of the calcium shell. The maximum thickness is, of course, the radius of the soma, which is 5 microns. If we have a smaller value than this it means that the calcium flux only comes into a cylindrical shell and therefore gives rise to higher concentrations.

The simulation starts with nominal values for the Na and K_DR channels of 100 Mho/m^2, and all other channels modulated down to zero.

Things to do:
    - Examine what each channel does on its own
    - Examine how to control the calcium influx
    - Once calcium influx is present, examine how it affects the K_Ca and K_AHP channels, and how in turn the cell dynamics are affected by calcium dynamics.
    - Can you get it to fire a burst of APs starting around 400ms?
    - Can you get it to fire just one action potential around 100 ms?
    - Can you get it to fire just one action potential around 400 ms?
    - Can you get it to fire with adaptation, that is, rate gets slower?
    - Can you get it to fire with facilitation, that is, rate gets faster?
    - It is a bit fiddly, as the simulation only runs 500 ms, but can you get it to do a couple of bursts?
"""

class SlideField():
    def __init__(self, name, field = 'modulation', initVal = 1.0, suffix = "Modulation", scale = 1.0, valmin = 0, valmax = 10.0 ):
        self.name = name
        self.path = '/model/elec/soma/' + name
        self.field = field
        assert( moose.exists( self.path ) )
        self.val = initVal
        self.suffix = suffix
        self.scale = scale
        self.valmin = valmin
        self.valmax = valmax

    def setVal( self, val ):
        moose.element( self.path ).setField( self.field, val * self.scale )
        self.val = val
        updateDisplay()

    def setChanMod( self, val ):
        global currentChanPath
        global currentChanMod
        currentChanPath = self.path
        currentChanMod = val
        self.val = val
        updateDisplay()

def features(Vtrace):
    stim_start = preStimTime
    stim_end = preStimTime + injectTime
    features_df = feature_range_df.copy()
    Inputcurr = 150e-12 #in A

    v = np.array(Vtrace) #in mV
    t = np.linspace(0,runtime, len(Vtrace)) #in s
    start_idx = (np.abs(t - stim_start)).argmin()
    end_idx = (np.abs(t - stim_end)).argmin()
    i = np.zeros(len(t))
    i[start_idx:end_idx] = Inputcurr*1e12 #in pA

    sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, filter = 9.9)
    sweep_ext.process_spikes()

    features_df['raw'] = ''
    # E_rest
    features_df.loc[f'E_rest_{Inputcurr*1e12}','raw'] = np.nanmean(v[:int(stim_start*len(v)/runtime)])
    # AP1_amp
    features_df.loc[f'AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #APp_amp
    features_df.loc[f'APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #AP1_width
    features_df.loc[f'AP1_width_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("width")[0]
    #APp_width
    features_df.loc[f'APp_width_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("width")[-2]
    #AP1_thresh
    features_df.loc[f'AP1_thresh_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("threshold_v")[0]
    #APp_thresh
    features_df.loc[f'APp_thresh_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("threshold_v")[-2]
    #AP1_lat
    features_df.loc[f'AP1_lat_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("threshold_t")[0] - stim_start
    #ISI1
    features_df.loc[f'ISI1_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_t")[1] - sweep_ext.spike_feature("peak_t")[0]
    #ISIl
    features_df.loc[f'ISIl_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("peak_t")[-1] - sweep_ext.spike_feature("peak_t")[-2]
    # #ISIavg
    # pt = sweep_ext.spike_feature("peak_t")
    # features_df.loc[f'ISIavg_{Inputcurr*1e12}','raw'] = np.nanmean([s-f for s,f in zip(pt[1:],pt[:-1])])
    features_df.loc[f'ISIavg_{Inputcurr*1e12}','raw'] = 'skip'
    #freq
    features_df.loc[f'freq_{Inputcurr*1e12}','raw'] = len(sweep_ext.spike_feature("peak_t"))/(stim_end - stim_start)
    #Adptn_id = 1-ISI1/ISIl
    features_df.loc[f'Adptn_id_{Inputcurr*1e12}','raw'] = 1 - features_df.loc[f'ISI1_{Inputcurr*1e12}','raw']/features_df.loc[f'ISIl_{Inputcurr*1e12}','raw']
    #fAHP_AP1_amp
    features_df.loc[f'fAHP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("fast_trough_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #fAHP_APp_amp
    features_df.loc[f'fAHP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("fast_trough_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    # #mAHP_AP1_amp
    # features_df.loc[f'mAHP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("slow_trough_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    features_df.loc[f'mAHP_AP1_amp_{Inputcurr*1e12}','raw'] = 'skip'
    #mAHP_APp_amp
    features_df.loc[f'mAHP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("slow_trough_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    # #mAHP_AP1_dur
    # features_df.loc[f'mAHP_AP1_dur_{Inputcurr*1e12}','raw'] = (sweep_ext.spike_feature("slow_trough_t")[0] - sweep_ext.spike_feature("peak_t")[0])/features_df.loc[f'ISI1_{Inputcurr*1e12}','raw']
    features_df.loc[f'mAHP_AP1_dur_{Inputcurr*1e12}','raw'] = 'skip'
    #mAHP_APp_dur = mAHP of second last spike (penultimate)
    features_df.loc[f'mAHP_APp_dur_{Inputcurr*1e12}','raw'] = (sweep_ext.spike_feature("slow_trough_t")[-2] - sweep_ext.spike_feature("peak_t")[-2])/features_df.loc[f'ISIl_{Inputcurr*1e12}','raw']
    # #ADP_AP1_amp
    # features_df.loc[f'ADP_AP1_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("adp_v")[0] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    features_df.loc[f'ADP_AP1_amp_{Inputcurr*1e12}','raw'] = 'skip'
    # #ADP_APp_amp
    # features_df.loc[f'ADP_APp_amp_{Inputcurr*1e12}','raw'] = sweep_ext.spike_feature("adp_v")[-2] - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    features_df.loc[f'ADP_APp_amp_{Inputcurr*1e12}','raw'] = 'skip'
    #mAHP_stimend_amp = within 50ms
    end50_idx = (np.abs(t - stim_end - 50e-3)).argmin()
    features_df.loc[f'mAHP_stimend_amp_{Inputcurr*1e12}','raw'] = np.min(v[end_idx:end50_idx]) - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']
    #sAHP_stimend_amp = within 200ms
    end200_idx = (np.abs(t - stim_end - 200e-3)).argmin()
    features_df.loc[f'sAHP_stimend_amp_{Inputcurr*1e12}','raw'] = np.min(v[end_idx:end200_idx]) - features_df.loc[f'E_rest_{Inputcurr*1e12}','raw']

    features_df = features_df.replace('', np.nan)

    features_df['rescaled'] = ''
    for i in features_df.index:
        if features_df.loc[i,'raw'] == 'skip':
            features_df.loc[i,'rescaled'] = 0
        elif np.isnan(features_df.loc[i,'raw']):
            features_df.loc[i,'rescaled'] = 10 #Penalty for not having the feature altogether
        else:
            Max = features_df.loc[i,'Max']
            Min = features_df.loc[i,'Min']
            features_df.loc[i,'rescaled'] = 2/(Max-Min)*(features_df.loc[i,'raw'] - Min) -1
    features_df = features_df.replace('', np.nan)

    features_df['cost'] = ''
    for i in features_df.index:
        if features_df.loc[i,'rescaled'] !='' and ~np.isnan(features_df.loc[i,'rescaled']):
            if features_df.loc[i,'rescaled'] >1 or features_df.loc[i,'rescaled'] <-1:
                features_df.loc[i,'cost'] = np.abs(features_df.loc[i,'rescaled'])
            else:
                features_df.loc[i,'cost'] = 0
    features_df = features_df.replace('', np.nan)

    return features_df

def makeModel():
    rdes = rd.rdesigneur(
        elecPlotDt = elecPlotDt,
        stealCellFromLibrary = True,
        verbose = False,
        #chanProto = [['make_glu()', 'glu'],['make_GABA()', 'GABA']],
        chanProto = [
            [ChP+'.Na_Chan()', 'Na'],
            [ChP+'.KDR_Chan()', 'K_DR'],
            [ChP+'.KA_Chan()', 'K_A'],
            [ChP+'.KM_Chan()', 'K_M'],
            [ChP+'.h_Chan()', 'h'],
            [ChP+'.CaT_Chan()', 'Ca_T'],
            [ChP+'.CaR_Chan()', 'Ca_R'],
            [ChP+'.CaL_Chan()', 'Ca_L'],
            [ChP+'.CaN_Chan()', 'Ca_N'],
            [ChP+'.Ca_Conc()', 'Ca_conc'],
            [ChP+'.KBK_Chan()', 'K_BK'],
            [ChP+'.KSK_Chan()', 'K_SK'],
        ],
        # cellProto syntax: ['ballAndStick', 'name', somaDia, somaLength, dendDia, dendLength, numDendSegments ]
        # The numerical arguments are all optional
        cellProto =
            [['somaProto', 'cellBase', diameter, diameter]],
        passiveDistrib = [[ '#', 'RM', str(RM), 'CM', str(CM), 'RA', str(RA) ]],
        chanDistrib = [
            ['Ca_conc', 'soma', 'tau', str(Ca_tau) ],
            ['Na', 'soma', 'Gbar', str(Na_Gbar) ],
            ['K_DR', 'soma', 'Gbar', str(K_DR_Gbar) ],
            ['K_A', 'soma', 'Gbar', str(K_A_Gbar) ],
            ['K_M', 'soma', 'Gbar', str(K_M_Gbar) ],
            ['h', 'soma', 'Gbar', str(h_Gbar) ],
            ['Ca_T', 'soma', 'Gbar', str(Ca_T_Gbar) ],
            ['Ca_R', 'soma', 'Gbar', str(Ca_R_Gbar) ],
            ['Ca_L', 'soma', 'Gbar', str(Ca_L_Gbar) ],
            ['Ca_N', 'soma', 'Gbar', str(Ca_N_Gbar) ],
            ['K_BK', 'soma', 'Gbar', str(K_BK_Gbar) ],
            ['K_SK', 'soma', 'Gbar', str(K_SK_Gbar) ],
        ],
        plotList = [['soma', '1','.', 'Vm'], ['soma', '1', 'Ca_conc', 'Ca']]
    )
    moose.element( '/library/Ca_conc' ).CaBasal=0.1e-3
    rdes.buildModel()
    for i in moose.wildcardFind( '/model/elec/soma/#[ISA=ChanBase]' ):
        i.modulation = 1e-6
    moose.element( '/model/elec/soma/Na' ).modulation = 1
    moose.element( '/model/elec/soma/K_DR' ).modulation = 1
    moose.element( '/model/elec/soma/Ca_conc' ).B = Ca_B

def main():
    global fields
    makeModel()
    fields.append( SlideField( 'Na', initVal = 1.0 ) )
    fields.append( SlideField( 'K_DR', initVal = 1.0 ) )
    fields.append( SlideField( 'K_A', initVal = 1.0 ) )
    fields.append( SlideField( 'K_M', initVal = 1.0 ) )
    fields.append( SlideField( 'h', initVal = 1.0 ) )
    fields.append( SlideField( 'Ca_T', initVal = 1.0 ) )
    fields.append( SlideField( 'Ca_R', initVal = 1.0 ) )
    fields.append( SlideField( 'Ca_L', initVal = 1.0 ) )
    fields.append( SlideField( 'Ca_N', initVal = 1.0 ) )
    fields.append( SlideField( 'K_SK', initVal = 1.0 ) )
    fields.append( SlideField( 'K_BK', initVal = 1.0 ) )
    fields.append( SlideField( 'Ca_conc', field = 'tau', initVal = 13.3, suffix = "(ms)" , scale = 0.001, valmax = 500.0 ) )
    fields[-1].name = 'Ca_tau'
    # fields.append( SlideField( 'Ca_conc', field = 'thick', initVal = 1e6*diameter/2.0, suffix = "(microns)", scale = 1e-6, valmax = 1e6*diameter / 2.0) )
    # fields[-1].name = 'Ca_thickness'
    fields.append( SlideField( 'Ca_conc', field = 'B', initVal = Ca_B, suffix = "(m^-3C^-1mol)", scale = 1, valmin = Ca_B*0.01, valmax = Ca_B*100) )
    fields[-1].name = 'Ca_B'
    fields.append( SlideField( '../soma', field = 'Rm', initVal = RM/sm_area*1e-6, suffix = "(Mohm)", scale = 1e6, valmin = 1, valmax = 100000 ) )
    fields[-1].name = 'Rm'
    fields.append( SlideField( '../soma', field = 'Cm', initVal = CM*sm_area*1e12, suffix = "(pF)", scale = 1e-12, valmin = 1, valmax = 1000 ) )
    fields[-1].name = 'Cm'
    fields.append( SlideField( '.', field = 'inject', initVal = 150, suffix = "(pA)", scale = 1e-12, valmax = 300 ) )
    fields[-1].name = 'Inject'
    warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")

    makeDisplay()
    # quit()
    plt.close()

def printSomaVm():
    print("This is somaVm" )

def runMod( mod ):
    moose.element( currentChanPath ).modulation = mod
    #print( 'path = {}, mod = {}'.format( currentChanPath, mod ) )
    moose.reinit()
    moose.start( preStimTime )
    moose.element( '/model/elec/soma' ).inject = fields[-1].val * fields[-1].scale
    moose.start( injectTime )
    moose.element( '/model/elec/soma' ).inject = 0
    moose.start( postStimTime )
    Vm = moose.element( '/model/graphs/plot0' ).vector
    Ca = moose.element( '/model/graphs/plot1' ).vector
    return Vm * 1000, Ca * 1e3

def updateDisplay():
    # Vm0, Ca0 = runMod( 1.0e-9 )
    Vm1, Ca1 = runMod( currentChanMod )
    #print len(Vm0), len(Ca0), len(Vm1), len( Ca1)
    # tplot[0].set_ydata( Vm0 )
    tplot[1].set_ydata( Vm1 )
    # tplot[2].set_ydata( Ca0 )
    tplot[3].set_ydata( Ca1 )

    print_parameters()
    features_df = features(Vm1)
    # print(features_df)
    print(f"Total cost = {np.nansum(features_df['cost'])}")
    # print('#############################################')

def print_parameters():
    soma = moose.element('/model/elec/soma')
    soma_area = np.pi*soma.diameter*soma.length
    print(f'RM={soma.Rm*soma_area}; CM={soma.Cm/soma_area}; RA={soma.Cm*soma_area/soma.length}')
    for i in moose.wildcardFind( '/model/elec/soma/#[ISA=ChanBase]' ):
        print(f'{i.name} Gbar = {i.Gbar/soma_area}')
    Caconcc = moose.element('/model/elec/soma/Ca_conc')
    print(f'Ca_Basal = {Caconcc.CaBasal}; Ca_B = {Caconcc.B}; Ca_tau = {Caconcc.tau}')

def doQuit( event ):
    # quit()
    plt.close()

def makeDisplay():
    global tplot
    global axes
    global sliders

    fig = plt.figure( figsize=(10,12) )
    t = np.arange( 0.0, runtime + elecPlotDt / 2.0, elecPlotDt ) * 1000 #ms
    ax0 = fig.add_subplot(311)
    plt.ylabel( 'Vm (mV)' )
    plt.ylim( -80, +40 )
    plt.xlabel( 'time (ms)' )
    plt.title( "Membrane potential vs. time." )
    ln, = ax0.plot( t, np.zeros(len(t)), 'r:', label='chan absent' )
    tplot.append(ln)
    ln, = ax0.plot( t, np.zeros(len(t)), 'b-', label='chan present' )
    tplot.append(ln)
    plt.legend()
    ####################################
    ax1 = fig.add_subplot(312)
    plt.ylabel( 'Ca (uM)' )
    plt.ylim( 0, 10 )
    plt.xlabel( 'time (ms)' )
    ln, = ax1.plot( t, np.zeros(len(t)), 'r:', label='chan absent' )
    tplot.append(ln)
    ln, = ax1.plot( t, np.zeros(len(t)), 'b-', label='chan present' )
    tplot.append(ln)
    plt.title( "Calcium concentration vs. time." )
    plt.legend()
    helpTextBox = plt.text( -0.14,-0.1, helpText, fontsize = 14, wrap = True, transform=ax1.transAxes, bbox=dict(facecolor='white', alpha=0.9), zorder = 10 )
    ####################################
    ax2 = fig.add_subplot(313)
    plt.axis('off')
    axcolor = 'palegreen'
    axStim = plt.axes( [0.02,0.005, 0.20,0.03], facecolor='green' )
    axReset = plt.axes( [0.25,0.005, 0.30,0.03], facecolor='blue' )
    axQuit = plt.axes( [0.60,0.005, 0.30,0.03], facecolor='blue' )

    for x in np.linspace( 0.05, 0.34, len(fields) ):
        axes.append( plt.axes( [0.25, x, 0.65, 0.015], facecolor=axcolor ) )
    #aInit = Slider( axAinit, 'A init conc', 0, 10, valinit=1.0, valstep=0.2)
    #rax = plt.axes([0.02, 0.05, 0.10, 0.28], facecolor="#EEEFFF")
    #mode = RadioButtons(rax, ('Channels', 'Other Parms'))

    stim = Button( axStim, 'View help', color = 'yellow' )

    reset = Button( axReset, 'Reset', color = 'cyan' )
    q = Button( axQuit, 'Quit', color = 'pink' )

    for i in range( len( axes ) ):
        if i < 2: # Na and K
            valinit = 1.0
        else:
            valinit = 0.0
        sliders.append( Slider( axes[i], fields[i].name+ " " + fields[i].suffix, 0.0, fields[i].valmax, valinit = fields[i].val) )
        if fields[i].field == 'modulation':
            sliders[-1].on_changed( fields[i].setChanMod )
        else:
            sliders[-1].on_changed( fields[i].setVal )

    def resetParms( event ):
        for i in sliders:
            i.reset()

    def toggleHelp( event ):
        if helpTextBox.get_visible():
            helpTextBox.set_visible( False )
            stim.label.set_text( "View Help" )
        else:
            helpTextBox.set_visible( True )
            stim.label.set_text( "Hide Help" )
        #plt.draw()
    #mh = modeHandler()

    #mode.on_clicked( mh.setMode )
    stim.on_clicked( toggleHelp )
    reset.on_clicked( resetParms )
    q.on_clicked( doQuit )

    updateDisplay()

    plt.show()

# Run the 'main' if this script is executed standalone.
if __name__ == '__main__':
        main()
