#exec(open('31659862.py').read())
# This question was asked on StackOverflow. Now solved.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import sys

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
t = np.arange(0.0, 1.0, 0.001)
a0 = 5
f0 = 3
s = a0*np.sin(2*np.pi*f0*t)
l, = plt.plot(t,s, lw=2, color='red')
plt.axis([0, 1, -10, 10])

axcolor = 'lightgoldenrodyellow'

d0      = 0.0
c       = 300000
z0      = d0/c

vmin    = -300.0
vmax    = 3000.0

zmin    = -0.01
zmax    = 2

axfreq  = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axz     = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

svlsr   = Slider(axfreq, 'VLSR', vmin, vmax, valinit=d0, valfmt=u'%1.1f')
sreds   = Slider(axz, 'z', zmin, zmax, valinit=z0, valfmt=u'%1.4f')

def donothing(val):
    pass
def update(val):
    global d0, z0
    delt = svlsr.val/c
    z = sreds.val

    svlsr.observers[svlsrcid] = donothing
    if z!=0.0:
        if z != z0:
            delt = z
            svlsr.set_val(z*c) #set_val causes infinite loop?? Now it doesn't.
    svlsr.observers[svlsrcid] = update

    d0 = delt
    z0    = z
    fac = 1.0 + delt
    l.set_xdata(t*fac)
    fig.canvas.draw_idle()

svlsrcid = svlsr.on_changed(update)
sredscid = sreds.on_changed(update)

plt.show()
