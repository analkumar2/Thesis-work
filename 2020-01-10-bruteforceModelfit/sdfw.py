import plotexp as pe
import os
for mm in os.listdir('../../WT step input cells'):
	try:
		pe.plotexp('../../WT step input cells' + mm, 150e-12)
	except:
		pass