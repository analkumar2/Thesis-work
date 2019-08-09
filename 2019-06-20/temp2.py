# exec(open('temp2.py').read())

import time
import numpy as np
import moose
import pylab
import rdesigneur as rd

t = time.time()
rdes = rd.rdesigneur(
    chanProto = [['Channels/kdr.xml', 'KDR']],
    chanDistrib = [
        ['KDR', 'soma', 'Gbar', '1200' ]],
    stimList = [['soma', '1', '.', 'inject', '(t>0.1 && t<0.2) * 1e-8' ]],
    plotList = [['soma', '1', '.', 'Vm', 'Membrane potential'],
    ['soma', '1', 'KDR', 'Ik', 'KDR current']]
)

rdes.buildModel()
moose.reinit()
moose.start( 0.3 )
print(time.time()-t)
rdes.display()
