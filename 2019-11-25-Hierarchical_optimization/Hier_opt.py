import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def Hierarchical_opt(function, p0):
    """
    function is the function to be minimized with the parameters to be optimized hierarchicaly to be given in order
    """
    ori_output = function(*p0)
    for i in len(p0):
        def tempfunc()
