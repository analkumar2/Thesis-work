import matplotlib.pyplot as plt
import moose

moose.Neutral('/library')

exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore1999).py').read())

Na_Chan('Na_Chan_Mig1999')
gateX_mig1999 = moose.element('/library/Na_Chan_Mig1999/gateX')
gateY_mig1999 = moose.element('/library/Na_Chan_Mig1999/gateY')

exec(open('../../Compilations/Kinetics/Na_Chan_(Migliore2018).py').read())
Na_Chan('Na_Chan_Mig2018')
gateX_mig2018 = moose.element('/library/Na_Chan_Mig2018/gateX')
gateY_mig2018 = moose.element('/library/Na_Chan_Mig2018/gateY')

v = np.linspace(-0.1,0.1,3000)

plt.plot(v, gateX_mig1999.tableA/gateX_mig1999.tableB, label='migliore1999 minf')
plt.plot(v, gateX_mig2018.tableA/gateX_mig2018.tableB, label='migliore2018 minf')
plt.legend()
plt.xlabel('Potential (V)')
plt.ylabel('inf')
plt.title('Mig2018 vs Mig2018')
plt.show()

plt.plot(v, 1/gateX_mig1999.tableB, label='migliore1999 mtau')
plt.plot(v, 1/gateX_mig2018.tableB, label='migliore2018 mtau')
plt.legend()
plt.xlabel('Potential (V)')
plt.ylabel('tau (s)')
plt.title('Mig2018 vs Mig2018')
plt.show()

plt.plot(v, gateY_mig1999.tableA/gateY_mig1999.tableB, label='migliore1999 hinf')
plt.plot(v, gateY_mig2018.tableA/gateY_mig2018.tableB, label='migliore2018 hinf')
plt.legend()
plt.xlabel('Potential (V)')
plt.ylabel('inf')
plt.title('Mig2018 vs Mig2018')
plt.show()

plt.plot(v, 1/gateY_mig1999.tableB, label='migliore1999 htau')
plt.plot(v, 1/gateY_mig2018.tableB, label='migliore2018 htau')
plt.legend()
plt.xlabel('Potential (V)')
plt.ylabel('tau (s)')
plt.title('Mig2018 vs Mig2018')
plt.show()