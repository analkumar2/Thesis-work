# exec(open('KinCompare.py').read())

exec(open('ChannelSlider.py').read())
Vmin = -0.100
Vmax = 0.100
Vdivs = 3000
v = np.linspace(Vmin,Vmax, Vdivs)

KDRX_Mig_inf = moose.element('/library/K_DR_chan/gateX').tableA/moose.element('/library/K_DR_chan/gateX').tableB
KAX_Mig_inf = moose.element('/library/K_A_chan/gateX').tableA/moose.element('/library/K_A_chan/gateX').tableB
KAY_Mig_inf = moose.element('/library/K_A_chan/gateY').tableA/moose.element('/library/K_A_chan/gateY').tableB

KDRX_Hoff_inf = 1/(1+np.exp((0.013-v)/0.011))
KAX_Hoff_inf = 1/(1+np.exp((0.011-v)/0.018))
KAY_Hoff_inf = 1/(1+np.exp((-0.056-v)/-0.008))

KDRX_Hoff_inf_m = (-0.0035*(v*1e3+30)/(np.exp((v*1e3+30)/-13)-1))/((-0.0035*(v*1e3+30)/(np.exp((v*1e3+30)/-13)-1)) + (0.0035*(v*1e3+30)/(np.exp((v*1e3+30)/13)-1)))
KAX_Hoff_inf_m = (-0.0035*(v*1e3+30)/(np.exp((v*1e3+30)/-13)-1))/((-0.0035*(v*1e3+30)/(np.exp((v*1e3+30)/-13)-1)) + (0.0035*(v*1e3+30)/(np.exp((v*1e3+30)/13)-1)))
KAY_Hoff_inf_m = 1/(1+np.exp((-0.056-v)/-0.008))

plt.plot(v, KDRX_Mig_inf, label='Mig2018')
plt.plot(v, KDRX_Hoff_inf, label='Hoff1997 exp')
plt.plot(v, KDRX_Hoff_inf_m**4, label='Hoff1997 model')
plt.title('KDR activation gate Inf')
plt.legend()
plt.show()
