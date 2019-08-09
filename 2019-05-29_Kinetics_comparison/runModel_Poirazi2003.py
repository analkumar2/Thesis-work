#exec(open('Kinetics_comparisons/runModel.py').read())

import io
import os
import sys
from neo.io import AxonIO
import quantities as pq
import csv
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
import numpy as np
import warnings
import moose
import pickle
import rdesigneur as rd
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
