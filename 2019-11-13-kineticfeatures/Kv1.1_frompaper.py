#exec(open('Kv1.1_frompaper.py').read())

#As seen in the data sheet (is different in the paper fig 8)

import numpy as np
import matplotlib.pyplot as plt

def mTauFunc(v, vBreak, offset, amp1, amp2, vh1, vh2, slope1, slpo2):
    sigswitch = 1/(1+np.exp((v-vBreak)/3))
    sig1 = sigswitch*amp1/(1+np.exp((v-vh1)/-slope1))
    sig2 = (1-sigswitch)*offset + (amp2-offset)/(1+np.exp((v-vh2)/slope2))
    return sig1+sig2
minf = lambda volt: 1 / (1 + np.exp((volt - -8.386847) / -9.503428))
hinf = lambda volt: 0.251120 + ( 0.785618 / (1 + np.exp((volt - -36.585932) / 13.927393)))
mtau = lambda volt: mTauFunc(volt, -57.076, -401.547, 35050.749, 14595.688, 3286.855, 0.640, 11.107, -40.568, 13.957)
htau = lambda volt: 44.316090 + ( 270.171206 / (1 + np.exp((volt - 16.341487) / 6.397418)))

v = np.linspace(-100,100,2000)
