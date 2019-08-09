# export PYTHONSTARTUP='pythonstartimport.py'

import numpy as np
import matplotlib.pyplot as plt
import csv
import moose
import pandas as pd

# Some common functions
def nearest_index(inp_array, value): #works only for strictly monotonic arrays
    return (np.abs(inp_array-value)).argmin()
