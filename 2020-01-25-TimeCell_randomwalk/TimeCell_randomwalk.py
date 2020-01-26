import numpy
import pandas as pd
import matplotlib.pyplot as plt
import numpy.random as nrm


class RandomWalkNeuron:

	def __init__(self, idx, prON):
		self.idx = idx
		self.prON = prON

	state = nrm.binomial(1,self.prON)

	def changestate(self):
		self.state = nrm.binomial(1,self.prON)