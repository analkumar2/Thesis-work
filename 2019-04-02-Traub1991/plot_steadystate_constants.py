#Author - Anal Kumar
#disclaimer - The code uses eval() function. Use at your own discretion.


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt
import Channelprotos as Cp #Channelprotos and rdesigneurProtos should be in the same directory as the pwd
try:
    moose.delete('/model')
    moose.delete('/library')
except:
    pass

Cp_list = dir(Cp) #getting a list of all variables and functions in the Channelprotos and rdesigneurProtos

for func in Cp_list: #Channelprotos
    if callable( eval('Cp.%s' %(func)) ): #checking if its a function or another variable
        moose.Neutral('/library') #setting up library
        try:
            eval('Cp.%s(\'%s\')' %(func, func)) # setting up the channel
            Chan = moose.element('/library/' + func)
        except:
            continue
        if Chan.className == 'HHChannel': #use this if the channel setup is HHChanel
            #Xgate
            if Chan.Xpower >=1 and Chan.instant != 1:
                Chanxgate = moose.element(Chan.path + '/gateX')
                min = Chanxgate.min
                max = Chanxgate.max
                minf = Chanxgate.tableA/Chanxgate.tableB
                tau = 1/Chanxgate.tableB
                plt.figure(0)
                plt.plot(np.linspace(min, max, len(minf)), minf, color = 'red', label = 'mInfinity')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('mInfinity')
                plt.title(func + ' mInfinity for its x gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_xgate_mInf.png')
                plt.figure(1)
                plt.plot(np.linspace(min, max, len(tau)), tau, color = 'blue', label = 'Tau')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Time constant (s)')
                plt.title(func + ' time constant for its x gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_xgate_tau.png')
                plt.close('all')
            
            # Ygate
            if Chan.Ypower >=1 and Chan.instant != 2:
                Chanygate = moose.element(Chan.path + '/gateY')
                min = Chanygate.min
                max = Chanygate.max
                minf = Chanygate.tableA/Chanygate.tableB
                tau = 1/Chanygate.tableB
                plt.figure(0)
                plt.plot(np.linspace(min, max, len(minf)), minf, color = 'red', label = 'mInfinity')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('mInfinity')
                plt.title(func + ' mInfinity for its y gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_ygate_mInf.png')
                plt.figure(1)
                plt.plot(np.linspace(min, max, len(tau)), tau, color = 'blue', label = 'Tau')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Time constant (s)')
                plt.title(func + ' time constant for its y gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_ygate_tau.png')
                plt.close('all')
                
            # Zgate
            if Chan.Zpower >=1 and Chan.instant != 4 and Chan.useConcentration == 0:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                minf = Chanzgate.tableA/Chanzgate.tableB
                tau = 1/Chanzgate.tableB
                plt.figure(0)
                plt.plot(np.linspace(min, max, len(minf)), minf, color = 'red', label = 'mInfinity')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('mInfinity')
                plt.title(func + ' mInfinity for its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate_mInf.png')
                plt.figure(1)
                plt.plot(np.linspace(min, max, len(tau)), tau, color = 'blue', label = 'Tau')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Time constant (s)')
                plt.title(func + ' time constant for its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate_tau.png')
                plt.close('all')
            elif Chan.Zpower >=1 and Chan.instant != 4 and Chan.useConcentration == 1:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                minf = Chanzgate.tableA/Chanzgate.tableB
                tau = 1/Chanzgate.tableB
                plt.figure(0)
                plt.plot(np.linspace(min, max, len(minf)), minf, color = 'red', label = 'mInfinity')
                plt.xlabel('Calcium concentration (mol/m^3)')
                plt.ylabel('mInfinity')
                plt.title(func + ' mInfinity for its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate_mInf.png')
                plt.figure(1)
                plt.plot(np.linspace(min, max, len(tau)), tau, color = 'blue', label = 'Tau')
                plt.xlabel('Calcium concentration (mol/m^3)')
                plt.ylabel('Time constant (s)')
                plt.title(func + ' time constant for its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate_tau.png')
                plt.close('all')
        print 'Cp_' + str(func)
        moose.delete('/library')