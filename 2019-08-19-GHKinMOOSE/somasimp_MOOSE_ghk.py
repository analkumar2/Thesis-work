# exec(open('somasimp_MOOSE_ghk.py').read())

import moose
import numpy as np
import matplotlib.pyplot as plt
import rdesigneur as rd

try:
    # [moose.delete(x) for x in ['/model', '/library']]
    [moose.delete(x) for x in ['/model']]
except:
    pass

ChP = 'simplechan_ghk'

F = 96485.3329
sm_diam = 500e-6
sm_len = 100e-6
sm_vol = np.pi/4*sm_diam**2*sm_len
sm_area = np.pi*sm_diam*sm_len
