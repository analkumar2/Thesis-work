# exec(open('chan_slider.py').read())
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

def ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, F):
    Inf = 1/(1+np.exp((v-vhalf_inf)/-slope_inf))
    Tau = B*(v+A/B)/(C+np.exp((v+D)/F))
    return [Inf,Tau]

fig, ax = plt.subplots(1,2)
plt.subplots_adjust(left=0.25, bottom=0.50)
Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
v = np.linspace(Vmin,Vmax, Vdivs)
Inf, Tau = ChanGate(v,-0.050,-0.004, 0.2,-3,100,0.05,-0.007 )
infl, = ax[0].plot(v, Inf, lw=2)
taul, = ax[1].plot(v, Tau, lw=2)

exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore2018).py').read())
moose.Neutral('/library')
Na_Chan('Na_Chan')
Mig2018hinf = moose.element('/library/Na_Chan/gateY').tableA/moose.element('/library/Na_Chan/gateY').tableB
Mig2018htau = 1/moose.element('/library/Na_Chan/gateY').tableB
hinfl_ori, = ax[0].plot(v, Mig2018hinf, lw=2)
htaul_ori, = ax[1].plot(v, Mig2018htau, lw=2)

ax[0].margins(x=0)
ax[1].margins(x=0)

axcolor = 'lightgoldenrodyellow'
axvhalf_inf = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
axslope_inf = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
axA = plt.axes([0.25, 0.20, 0.65, 0.03], facecolor=axcolor)
axB = plt.axes([0.25, 0.25, 0.65, 0.03], facecolor=axcolor)
axC = plt.axes([0.25, 0.30, 0.65, 0.03], facecolor=axcolor)
axD = plt.axes([0.25, 0.35, 0.65, 0.03], facecolor=axcolor)
axF = plt.axes([0.25, 0.40, 0.65, 0.03], facecolor=axcolor)

svhalf_inf = Slider(axvhalf_inf, 'vhalf_inf', -0.100, 0, valinit=-0.05, valstep=0.001, valfmt='%1.4f')
sslope_inf = Slider(axslope_inf, 'slope_inf', -0.100, 0, valinit=-0.004, valstep=0.001, valfmt='%1.4f')
sA = Slider(axA, 'A', -100, 100, valinit=-0.015, valstep=0.001, valfmt='%1.4f')
sB = Slider(axB, 'B', -100, 100, valinit=-0.015, valstep=0.001, valfmt='%1.4f')
sC = Slider(axC, 'C', -100, 100, valinit=-0.015, valstep=0.001, valfmt='%1.4f')
sD = Slider(axD, 'D', -100, 100, valinit=-0.015, valstep=0.001, valfmt='%1.4f')
sF = Slider(axF, 'F', -100, 100, valinit=-0.015, valstep=0.001, valfmt='%1.4f')

def update(val):
    vhalf_inf = svhalf_inf.val
    slope_inf = sslope_inf.val
    A = sA.val
    B = sB.val
    C = sC.val
    D = sD.val
    F = sF.val
    I,T = ChanGate(v,vhalf_inf, slope_inf, A, B, C, D, F)
    infl.set_ydata(I)
    taul.set_ydata(T)
    fig.canvas.draw_idle()


svhalf_inf.on_changed(update)
sslope_inf.on_changed(update)
sA.on_changed(update)
sB.on_changed(update)
sC.on_changed(update)
sD.on_changed(update)
sF.on_changed(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')


def reset(event):
    svhalf_inf.reset()
    sslope_inf.reset()
    sA.reset()
    sB.reset()
    sC.reset()
    sD.reset()
    sF.reset()
button.on_clicked(reset)

plt.show()
