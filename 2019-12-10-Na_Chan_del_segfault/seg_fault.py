# exec(open('seg_fault.py').read())

import moose
import pylab
import rdesigneur as rd


# Wrapper function so that the model can be build and run again and again
def rdeswrapper():
    # Deleting any previous run of the model
    try:
        # [moose.delete(x) for x in ['/model', '/library']]
        moose.delete('/model')
    except:
        pass
    ######################################

    rdes = rd.rdesigneur(
        chanProto = [['make_HH_Na()', 'Na'], ['K_A_Chan_(Migliore2018)_ghk.K_A_Chan()', 'K']],
        chanDistrib = [
            ['K', 'soma', 'Gbar', '2000' ],
            ['Na', 'soma', 'Gbar', '100' ],],
        stimList = [['soma', '1', '.', 'inject', '(t>0.1 && t<0.2) * 1e-8' ]],
        plotList = [['soma', '1', '.', 'Vm', 'Membrane potential']]
    )
    rdes.buildModel()
    moose.reinit()
    moose.start( 0.3 )
    rdes.display()
    return rdes

# # Initial run
# print('Initial run')
# rdeswrapper()

# Delete library and run
moose.delete('/library')
print('After libsrary deletion and re-build and re-run')
rdeswrapper()

# Delete Na and run
moose.delete('/library/Na')
print('After libsrary/Na deletion and re-build and re-run')
rdeswrapper()
