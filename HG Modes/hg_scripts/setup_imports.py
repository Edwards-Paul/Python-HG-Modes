# Mass importing file. Maybe not good programming discipline but ok...

# EXTERNAL IMPORTS

## Maths
import numpy as np, matplotlib.pyplot as plt, cmath as cm, scipy
from sympy.parsing import mathematica as mc
from numpy import sin as sin, pi as pi, angle
#from numpy import sqrt as sqrt
from scipy.special import erf as erf, gamma as gamma, comb
import mpmath as mp
import scipy.io
from math import pi, log, exp, sin, cos, atan, e, radians, degrees,factorial as Factorial
import cmath
from cmath import sqrt as Sqrt

inf=np.inf

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

# INTERNAL IMPORTS
from hg_scripts import paulisa as pl, pl_plot as plplt,plback as plb, tophat_integration_AW_2 as th
#import signals_calc as sig
from hg_scripts.transverse_coord_transform import *
from hg_scripts.build_modes import update_modes
from hg_scripts.beam_from_mat_file import coeff_from_mat