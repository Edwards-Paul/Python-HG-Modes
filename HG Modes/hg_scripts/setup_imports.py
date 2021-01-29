##############
#Mass importing file. Maybe not good programming discipline but ok...
##############


# EXTERNAL IMPORTS
## Maths
import numpy as np, matplotlib.pyplot as plt, cmath as cm, scipy
from sympy.parsing import mathematica as mc
from numpy import sin as sin, angle
#from numpy import sqrt as sqrt
from scipy.special import erf as erf, gamma as gamma, comb
import mpmath as mp
import scipy.io
from math import pi, log, exp, sin, cos, atan, radians, degrees,factorial as Factorial
import cmath
from cmath import sqrt as Sqrt
## Times
from time import process_time
import pandas as pd
from pprint import pprint
###progress bar
import time
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer, UnknownLength
## Etc.
import re
from copy import copy as cp
from ast import literal_eval

# CONSTANTS
inf=np.inf
from numpy import pi as pi
from math import e

# INTERNAL IMPORTS
from hg_scripts import paulisa as pl, pl_plot as plplt,plback as plb, tophat_integration_AW as th
#import signals_calc as sig
from hg_scripts.transverse_coord_transform import * #class: item; method: transform_x
from hg_scripts.build_modes import update_modes #method: update_modes
from hg_scripts.beam_from_mat_file import coeff_from_mat, truncated_tophat
from hg_scripts.plot_signals import * #method: plot_dws,plot_lpsT,plot_lpsR
from hg_scripts.mathematica_to_matrix import * #method: parse_math_for_python
from hg_scripts.signals_calc import *