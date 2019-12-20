# exec(open('Na_chan_slider.py').read())
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def Na_Chan_minf(v,vshiftm, slopem):
    #vshiftm and slopem should be dimensionless but think of them to be in mV terms
    mAlpha = (0.182 * (v*1e3- (-38+vshiftm)))/(1-(np.exp(-(v*1e3- (-38+vshiftm))/slopem)))
    mBeta  = (0.124 * (-v*1e3 + (-38+vshiftm)))/(1-(np.exp(-(-v*1e3 + (-38+vshiftm))/slopem)))
    mInf = mAlpha/(mAlpha + mBeta)
    return mInf

fig, ax = plt.subplots()
plt.subplots_adjust(left=0.25, bottom=0.25)
v = np.arange(-0.100, 0.100, 0.0003)
mInf = Na_Chan_minf(v,0,6)
l, = plt.plot(v, mInf, lw=2)
ax.margins(x=0)

axcolor = 'lightgoldenrodyellow'
axslopem = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axshiftm = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)

sslopem = Slider(axslopem, 'slopem', -10, 10, valinit=6, valstep=0.1)
sshiftm = Slider(axshiftm, 'shiftm', -10, 10, valinit=0)


def update(val):
    shiftm = sshiftm.val
    slopem = sslopem.val
    l.set_ydata(Na_Chan_minf(v,shiftm,slopem))
    fig.canvas.draw_idle()


sslopem.on_changed(update)
sshiftm.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    sslopem.reset()
    sshiftm.reset()
button.on_clicked(reset)

plt.show()
