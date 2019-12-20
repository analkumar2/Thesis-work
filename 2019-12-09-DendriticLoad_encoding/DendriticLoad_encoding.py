# exec(open('DendriticLoad_encoding.py').read())

import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import time

numdendSegments = 10
dendcomptLen = 10e-6
dendcomptDia = 5e-6
dendRM = 3
dendRA = 1
dendCM = 0.75e-2

numsomaSegments = 3
somacomptLen = 10e-6
somacomptDia = 20e-6
somaRM = 3
somaRA = 1
somaCM = 0.75e-2

numAISSegments = 2
AIScomptLen = 2e-6
AIScomptDia = 1e-6
AISRM = 3
AISRA = 1
AISCM = 0.75e-2

numaxonSegments = 10
axoncomptLen = 10e-6
axoncomptDia = 1e-6
axonRM = 3
axonRA = 1
axonCM = 0.75e-2


def makeCelllProto():
    Celll = moose.Neuron( '/library/Celll' )

    prev = rd.buildCompt( Celll, 'dend0', RM = dendRM, RA = dendRA, CM = dendCM, dia = dendcomptDia, x=0, dx=dendcomptLen)
    x = dendcomptLen
    for i in np.arange(1, numdendSegments):
        compt = rd.buildCompt( Celll, 'dend' + str(i), RM=dendRM, RA=dendRA, CM=dendCM, x=x, dx=dendcomptLen, dia=dendcomptDia )
        moose.connect( prev, 'axial', compt, 'raxial' )
        prev = compt
        x += dendcomptLen

    prev = rd.buildCompt( Celll, 'soma0', RM = somaRM, RA = somaRA, CM = somaCM, dia = somacomptDia, x=x, dx=somacomptLen)
    moose.connect( prev, 'axial', compt, 'raxial' )
    x += somacomptLen
    for i in np.arange(1, numsomaSegments):
        compt = rd.buildCompt( Celll, 'soma' + str(i), RM=somaRM, RA=somaRA, CM=somaCM, x=x, dx=somacomptLen, dia=somacomptDia )
        moose.connect( prev, 'axial', compt, 'raxial' )
        prev = compt
        x += somacomptLen

    prev = rd.buildCompt( Celll, 'AIS0', RM = AISRM, RA = AISRA, CM = AISCM, dia = AIScomptDia, x=x, dx=AIScomptLen)
    moose.connect( prev, 'axial', compt, 'raxial' )
    x += AIScomptLen
    for i in np.arange(1, numAISSegments):
        compt = rd.buildCompt( Celll, 'AIS' + str(i), RM=AISRM, RA=AISRA, CM=AISCM, x=x, dx=AIScomptLen, dia=AIScomptDia )
        moose.connect( prev, 'axial', compt, 'raxial' )
        prev = compt
        x += AIScomptLen

    prev = rd.buildCompt( Celll, 'axon0', RM = axonRM, RA = axonRA, CM = axonCM, dia = axoncomptDia, x=x, dx=axoncomptLen)
    moose.connect( prev, 'axial', compt, 'raxial' )
    x += axoncomptLen
    for i in np.arange(1, numaxonSegments):
        compt = rd.buildCompt( Celll, 'axon' + str(i), RM=axonRM, RA=axonRA, CM=axonCM, x=x, dx=axoncomptLen, dia=axoncomptDia )
        moose.connect( prev, 'axial', compt, 'raxial' )
        prev = compt
        x += axoncomptLen

    return Celll

def modelwrapper(Na_soma_Gbar, Na_AIS_Gbar, K_soma_Gbar, K_AIS_Gbar):
    # ttt = time.time()
    try:
        moose.delete('/model')
        # moose.delete('/library')
    except:
        pass

    moose.Neutral( '/library' )
    makeCelllProto()

    rdes = rd.rdesigneur(
            cellProto = [['elec', 'Celll'],],
            chanProto = [['make_HH_Na()', 'Na'],
                        ['make_HH_Na()', 'K'],],
            chanDistrib = [['Na', '#soma#', 'Gbar', f'{Na_soma_Gbar}' ],
                        ['Na', '#AIS#', 'Gbar', f'{Na_AIS_Gbar}' ],
                        ['K', '#soma#', 'Gbar', f'{K_soma_Gbar}' ],
                        ['K', '#AIS#', 'Gbar', f'{K_AIS_Gbar}' ],],
            stimList = [['soma1', '1', '.', 'inject', '(t>0.5 && t<1) * 2e-11' ]],
            plotList = [['soma1', '1', '.', 'Vm', 'Soma1 Membrane potential'],
                        ['AIS20', '1', '.', 'Vm', 'AIS20 Membrane potential'],],
            # moogList = [['#', '1', '.', 'Vm', 'Vm (mV)']],
    )

    rdes.buildModel()
    moose.reinit()
    moose.start(1.5)
    # print(time.time()-ttt)
    return rdes

    # rdes.displayMoogli( 0.001, 2, rotation = 0, azim = 0, elev = 0.0)


for i in range(1000):
    Na_soma_Gbar, Na_AIS_Gbar, K_soma_Gbar, K_AIS_Gbar = np.random.randint(10000, size=4)
    print(i)
    rdes = modelwrapper(Na_soma_Gbar, Na_AIS_Gbar, K_soma_Gbar, K_AIS_Gbar)
    Vvec = moose.element('/model/graphs/plot0').vector
    if np.mean(Vvec[int(len(Vvec)/1.5*0.1):int(len(Vvec)/1.5*0.5)]) < -0.050:
        print('Na_soma_Gbar, Na_AIS_Gbar, K_soma_Gbar, K_AIS_Gbar = ', Na_soma_Gbar, Na_AIS_Gbar, K_soma_Gbar, K_AIS_Gbar)
        rdes.display()
