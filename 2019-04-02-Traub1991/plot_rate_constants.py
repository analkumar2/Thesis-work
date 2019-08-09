#Author - Anal Kumar
#disclaimer - The code uses eval() function. Use at your own discretion.


import moose
import rdesigneur as rd
import numpy as np
import matplotlib.pyplot as plt

import Channelprotos as Cp #Channelprotos and rdesigneurProtos should be in the same directory as the pwd
import rdesigneurProtos as rp

Cp_list = dir(Cp) #getting a list of all variables and functions in the Channelprotos and rdesigneurProtos
rp_list = dir(rp)

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
                alpha = Chanxgate.tableA
                beta = Chanxgate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its x gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_xgate.png')
            elif Chan.Xpower >=1 and Chan.instant == 1:
                Chanxgate = moose.element(Chan.path + '/gateX')
                min = Chanxgate.min
                max = Chanxgate.max
                value = Chanxgate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its x gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_xgate.png')
            
            # Ygate
            if Chan.Ypower >=1 and Chan.instant != 2:
                Chanygate = moose.element(Chan.path + '/gateY')
                min = Chanygate.min
                max = Chanygate.max
                alpha = Chanygate.tableA
                beta = Chanygate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its y gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_ygate.png')
            elif Chan.Ypower >=1 and Chan.instant == 2:
                Chanygate = moose.element(Chan.path + '/gateY')
                min = Chanygate.min
                max = Chanygate.max
                value = Chanygate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its y gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_ygate.png')
                
            # Zgate
            if Chan.Zpower >=1 and Chan.instant != 4 and Chan.useConcentration == 0:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                alpha = Chanzgate.tableA
                beta = Chanzgate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate.png')
            elif Chan.Zpower >=1 and Chan.instant != 4 and Chan.useConcentration == 1:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                alpha = Chanzgate.tableA
                beta = Chanzgate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Calcium concentration (mol/m^3)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate.png')
            elif Chan.Zpower >=1 and Chan.instant == 4 and Chan.useConcentration == 0:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                value = Chanzgate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate.png')
            elif Chan.Zpower >=1 and Chan.instant == 4 and Chan.useConcentration == 1:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                value = Chanzgate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Calcium concentration (mol/m^3)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its z gate')
                plt.legend()
                plt.savefig('./gateparams/Cp_' + func + '_zgate.png')
                
        if Chan.className == 'CaConc':
            pass
        print 'Cp_' + str(func)
        moose.delete('/library')

        
for func in rp_list: #rdesigneurProtos
    if callable(eval('rp.%s' %(func))):
        moose.Neutral('/library') #setting up library
        try:
            eval('rp.%s(\'%s\')' %(func, func)) # setting up the channel
            Chan = moose.element('/library/' + func)
        except:
            continue
        if Chan.className == 'HHChannel': #use this if the channel setup is HHChanel
            #Xgate
            if Chan.Xpower >=1 and Chan.instant != 1:
                Chanxgate = moose.element(Chan.path + '/gateX')
                min = Chanxgate.min
                max = Chanxgate.max
                alpha = Chanxgate.tableA
                beta = Chanxgate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its x gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_xgate.png')
            elif Chan.Xpower >=1 and Chan.instant == 1:
                Chanxgate = moose.element(Chan.path + '/gateX')
                min = Chanxgate.min
                max = Chanxgate.max
                value = Chanxgate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its x gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_xgate.png')
            
            # Ygate
            if Chan.Ypower >=1 and Chan.instant != 2:
                Chanygate = moose.element(Chan.path + '/gateY')
                min = Chanygate.min
                max = Chanygate.max
                alpha = Chanygate.tableA
                beta = Chanygate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its y gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_ygate.png')
            elif Chan.Ypower >=1 and Chan.instant == 2:
                Chanygate = moose.element(Chan.path + '/gateY')
                min = Chanygate.min
                max = Chanygate.max
                value = Chanygate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its y gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_ygate.png')
                
            # Zgate
            if Chan.Zpower >=1 and Chan.instant != 4 and Chan.useConcentration == 0:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                alpha = Chanzgate.tableA
                beta = Chanzgate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its z gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_zgate.png')
            elif Chan.Zpower >=1 and Chan.instant != 4 and Chan.useConcentration == 1:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                alpha = Chanzgate.tableA
                beta = Chanzgate.tableB - alpha
                plt.figure()
                plt.plot(np.linspace(min, max, len(alpha)), alpha, color = 'red', label = 'alpha')
                plt.plot(np.linspace(min, max, len(beta)), beta, color = 'blue', label = 'beta')
                plt.xlabel('Calcium concentration (mol/m^3)')
                plt.ylabel('Rate constant (1/s)')
                plt.title(func + ' rate constants for its z gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_zgate.png')
            elif Chan.Zpower >=1 and Chan.instant == 4 and Chan.useConcentration == 0:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                value = Chanzgate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Membrane potential (V)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its z gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_zgate.png')
            elif Chan.Zpower >=1 and Chan.instant == 4 and Chan.useConcentration == 1:
                Chanzgate = moose.element(Chan.path + '/gateZ')
                min = Chanzgate.min
                max = Chanzgate.max
                value = Chanzgate.tableA
                plt.figure()
                plt.plot(np.linspace(min, max, len(value)), value, color = 'red')
                plt.xlabel('Calcium concentration (mol/m^3)')
                plt.ylabel('Gate value')
                plt.title(func + ' value of its z gate')
                plt.legend()
                plt.savefig('./gateparams/rp_' + func + '_zgate.png')
                
        if Chan.className == 'CaConc':
            pass
        print 'rp_' + str(func)
        moose.delete('/library')
