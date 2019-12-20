
import time
import moose
import pylab
import rdesigneur as rd
import numpy as np

rdes = rd.rdesigneur(
    chanProto = [['make_HH_Na()', 'Na'], ['make_HH_K()', 'K']],
    chanDistrib = [
        ['Na', 'soma', 'Gbar', '1200' ],
        ['K', 'soma', 'Gbar', '360' ]],
    stimList = [['soma', '1', '.', 'inject', '(t>1 && t<2) * 1e-8' ]],
    plotList = [['soma', '1', '.', 'Vm', 'Membrane potential']]
)
elecid_ori = rdes.elecid

def wraa(i):
	try:
		moose.delete('/model')
		rdes.elecid = moose.element(elecid_ori)
	except:
		pass

	rdes.buildModel()
	moose.reinit()
	moose.start( 3 )
	print(i)
	print((time.time()-ttt)/i)
	# rdes.display()

ttt = time.time()
iii = np.arange(1, 30001)
[wraa(i) for i in iii]
