# Calculating and plotting intensity and phase of propagating HG modes
# q param. goes from i(10) -> i(10+100)
# evaluate HG modes again at new q-param. and add together.
# optical params described by waist size (w0) and waist location into q parameter (q(z)=i*Zr+z-z0=q0+z-z0, q0=i*Zr)
# w(z) spot size - radius at which intensity is 1/e**2 max intensity I(0)
# ========================================================================================================

## IMPORTS:

# Math imports
from math import pi, log, exp, sin, cos, atan, sqrt, e, factorial, radians, degrees
from numpy.polynomial.hermite import hermval
from scipy import integrate
import cmath
import numpy as np
# Plot imports
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


# --------------------------------------------------------------------------------------------------------
# OPTICAL PARAMETERS PASSED TO CALCULATION
class Params:

    def __init__(self, wavelength, w0, z0):
        self.wavelength = wavelength
        self.w0 = w0  # Beam waist size, sqrt(zr*wavelength/pi)~1mm
        self.z0 = z0  # waist location, I ~ 1/e**2 ! 0.13
        self.Zr = pi * w0 ** 2 / wavelength  # Rayleigh Range, near field = 2 Zr = 6
        self.q0 = (1j) * self.Zr
        self.k = 2 * pi / (wavelength)  # Wavenumber

    def __str__(self):
        return '\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}'.format('wavelength=', self.wavelength,
                                                       'w0=', self.w0,
                                                       'z0=', self.z0,
                                                       'Zr=', self.Zr,
                                                       'q0=', self.q0,
                                                       'k=', self.k)

    def getWavelength(self):
        return self.wavelength

    def getK(self):
        return self.k

    def getW0(self):
        return self.w0

    def getZ0(self):
        return self.z0

    def getZr(self):
        return self.Zr

    def getQ0(self):
        return self.q0


# --------------------------------------------------------------------------------------------------------
## DEFAULT CONSTANTS :

defaultParams = Params(1064e-9, 1e-3, 0)


# defaultParams = Params(1064e-9, 1e-2, 0)


# --------------------------------------------------------------------------------------------------------
## FUNCTIONS OF OPTICAL PARAMETERS:


## Desc change in the radius of curvature (Rc)
def radius_curvature(z, params):
    # r=z-z0+(Zr**2/(z-z0))

    if z == params.getZ0():
        r = float('inf')
    else:
        r = z - params.getZ0() + (params.getZr() ** 2) / (z - params.getZ0())
    # r = z * (1 + (params.getZr() / z) ** 2)

    return r


## Gouy phase
def gouy_phase(z, order, params):
    Gouy = atan((z - params.getZ0()) / params.getZr())
    return Gouy


## Phase lag
def phase_lag(z, order, params):
    # phaselag = (order + 1) * atan(z / params.getZr())
    phaselag = (order + 1) * gouy_phase(z, order, params)
    return phaselag


## Spot size (LR eq. 9.16)
def w(z, params):
    w = params.getW0() * sqrt(1 + ((z - params.getZ0()) / params.getZr()) ** 2)
    return w


## q-parameter
def q(z, params):
    q = params.getQ0() + z - params.getZ0()
    return q


# --------------------------------------------------------------------------------------------------------
##INITIAL INPUT:
# --------------------------------------------------------------------------------------------------------

##Get user input modes, create 2-d modes array, show array
def modes(*argv):
    # get number of modes passed by user
    NumberModes = (len(argv))

    # create lists for n,m,c of length equal to number of modes
    listN = [None] * NumberModes
    listM = [None] * NumberModes
    listC = [None] * NumberModes

    # parse args into lists (n,m, & c) to set modes
    for i in range(0, len(argv)):
        listN[i], listM[i], listC[i] = argv[i]

    # get a modes 2-d array created from lists of n,m,c of required size
    return (create_modes(listN, listM, listC, NumberModes))


##Create the modes 2-d array
def create_modes(listN, listM, listC, NumberModes):
    # get max n and m to create modes grid with size based on user highest modes
    MaxN = max(listN)
    MaxM = max(listM)

    # get number of rows/cols for modes grid (plus 1 for 0 index..)
    rows = MaxN + 1
    cols = MaxM + 1

    # initialize 2-d array
    modes = [[0 for m in range(cols)] for n in range(rows)]
    # iterate through lists to set modes in grid
    for i in range(0, NumberModes):
        modes[listN[i]][listM[i]] = listC[i]

    return (modes)


##Print modes
def show_modes(modes):
    if not modes:
        print("No modes entered.")

    else:
        rows = len(modes)
        cols = len(modes[0])

        colList = []

        for m in range(cols):
            colList.append(m)
        print("n\m " + str(colList))

        for n in range(rows):
            List = []
            for m in range(cols):
                List.append(modes[n][m])
                if m == cols - 1:
                    print(str(n) + "   " + str(List))


# --------------------------------------------------------------------------------------------------------
##Get X,Y plane for calculation
class Plane:

    def __init__(self, xmin, xmax, xpoints, ymin, ymax, ypoints):
        # minimum and maximum x values with # points
        self.xmin = xmin
        self.xmax = xmax
        self.xpoints = xpoints
        # minimum and maximum y values with # points
        self.ymin = ymin
        self.ymax = ymax
        self.ypoints = ypoints

    def __str__(self):
        return '\n{}{},{}{},{}{},{}{}\n{}{},{}{},{}{},{}{}'.format('xmin=', self.xmin,
                                                                   'xmax=', self.xmax,
                                                                   'xpoints=', self.xpoints,
                                                                   'x step size = ',
                                                                   abs(self.xmax - self.xmin) / self.xpoints,
                                                                   'ymin=', self.ymin,
                                                                   'ymax=', self.ymax,
                                                                   'ypoints=', self.ypoints,
                                                                   'y step size = ',
                                                                   abs(self.ymax - self.ymin) / self.ypoints)

    def getXmin(self):
        return self.xmin

    def getXmax(self):
        return self.xmax

    def getXpoints(self):
        return self.xpoints

    def getYmin(self):
        return self.ymin

    def getYmax(self):
        return self.ymax

    def getYpoints(self):
        return self.ypoints

    # x,y vectors from limits and # points
    def getX(self):
        return np.linspace(self.xmin, self.xmax, self.xpoints)

    def getY(self):
        return np.linspace(self.ymin, self.ymax, self.ypoints)


# --------------------------------------------------------------------------------------------------------
# DEFAULT PLANE OF CALCULATION
defaultPlane = Plane(-2e-2, 2e-2, 1000, -2e-2, 2e-2, 1000)

# --------------------------------------------------------------------------------------------------------
##PLANAR CALCULATIONS OF AMPLITUDE AND PHASE

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def calculate(params, plane, modes, z):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude(params, X, Y, z, modes)

    return (amp)

# --------------------------------------------------------------------------------------------------------
def calculate_tilted(params, plane, modes, z, tilt,shift):

    # map z->z-(x+a)*sin(alpha), x->(x+a)*cos(alpha)+z*sin(alpha)
    # Make coordinate change for tilt

    X, Y = np.meshgrid( (plane.getX()+shift)*cos(tilt), plane.getY())
    Z = z-(X+shift)*sin(tilt)

    amp = amplitude_tilted(params, X, Y, Z, tilt, modes)

    return (amp)
# --------------------------------------------------------------------------------------------------------
# Calculate peak and its location from result
class PeakInt:

    def __init__(self, result):
        intensity = abs(result.amp) ** 2

        map(max, intensity)
        list(map(max, intensity))

        self.peak = max(map(max, intensity))
        self.loc = np.where(intensity == self.peak)

        self.x = self.loc[1] * abs((result.plane.xmax - result.plane.xmin) / result.plane.xpoints) + result.plane.xmin
        self.y = self.loc[0] * abs((result.plane.ymax - result.plane.ymin) / result.plane.ypoints) + result.plane.ymin



class PeakAmp:

    def __init__(self, result):
        amplitude = (result.amp)

        map(max, amplitude)
        list(map(max, amplitude))

        self.peak = max(map(max, amplitude))
        self.loc = np.where(amplitude == self.peak)

        self.x = self.loc[1] * abs((result.plane.xmax - result.plane.xmin) / result.plane.xpoints) + result.plane.xmin
        self.y = self.loc[0] * abs((result.plane.ymax - result.plane.ymin) / result.plane.ypoints) + result.plane.ymin




##AMPLITUDE CALCULATIONS:

# Amplitutude calculation from w0,Zr. (LR eq. 9.26)
def amplitude(params, x, y, z, modes):
    # Unm a sum over modes
    UnmSum = 0
    rows = len(modes)
    cols = len(modes[0])

    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            coeff = 0
            # n array for hermpol
            carrN = [0] * rows
            carrN[n] = modes[n][m]

            # m array for hermpol
            carrM = [0] * cols
            carrM[m] = modes[n][m]

            # Avoid double-counting coefficient
            if (carrN[n] != 0):
                carrN[n] = 1

            order = n + m

            Unm = (1 / sqrt(2 ** (order - 1) * factorial(n) * factorial(m) * pi)) * \
                  (1 / w(z, params)) * e ** ((1j) * (order + 1) * gouy_phase(z, order, params)) * \
                  e ** ((-1j) * ((params.getK() * (x ** 2 + y ** 2) / (2.0 * radius_curvature(z, params)))) -
                        ((x ** 2 + y ** 2) / ((w(z, params)) ** 2))) * \
                  herm_poly(n, x, carrN, z, params) * \
                  herm_poly(m, y, carrM, z, params)

            # Add each to result
            UnmSum += Unm

    return (UnmSum)


#--------------------------------------------------------
# TILTED BEAM AMPLITUDE
# map z->z-x*sin(alpha), x->x*cos(alpha)
def amplitude_tilted(params, x, y, z, tilt, modes):


    # Unm a sum over modes
    UnmSum = 0
    rows = len(modes)
    cols = len(modes[0])

    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            coeff = 0
            # n array for hermpol
            carrN = [0] * rows
            carrN[n] = modes[n][m]

            # m array for hermpol
            carrM = [0] * cols
            carrM[m] = modes[n][m]

            # Avoid double-counting coefficient
            if (carrN[n] != 0):
                carrN[n] = 1

            order = n + m

            Unm = (1 / sqrt(2 ** (order - 1) * factorial(n) * factorial(m) * pi)) * \
                  (1 / w(z, params)) * e ** ((1j) * (order + 1) * gouy_phase(z, order, params)) * \
                  e ** ((-1j) * ((params.getK() * (x ** 2 + y ** 2) / (2.0 * radius_curvature(z, params)))) -
                        ((x ** 2 + y ** 2) / ((w(z, params)) ** 2))) * \
                  herm_poly(n, x, carrN, z, params) * \
                  herm_poly(m, y, carrM, z, params)

            # Add each to result
            UnmSum += Unm

    return (UnmSum)




## Get herm polys from modes
# mode: working mode (as Calculate iterates through n,m grid
# coord: x or y
# carr: coefficient array
def herm_poly(mode, coord, carr, z, params):
    herm = hermval(coord * sqrt(2.0) / w(z, params), carr)
    return herm


## PHASE CALCULATIONS:

def phase(params, x, y, z, modes):
    return (np.angle(amplitude(params, x, y, z, modes)))


# =========================================================================================================
## PRINTING DEFAULTS AND USAGE
def defaults():
    print("DEFAULT PARAMS (PauLisa.defaultParams)\
    \n wavelength =" + str(defaultParams.wavelength) + "\
    m\n waist size(w0) =" + str(defaultParams.w0) + "\
    m\n z0 =" + str(defaultParams.z0) + "\
    m\n Rayleigh Range (Zr) =" + str(defaultParams.Zr) + "m")

    print("\n\nDEFAULT X,Y PLANE (PauLisa.defaultPlane)\
    \n x: " + str(defaultPlane.xmin) + "m to " + str(defaultPlane.xmax) + "m with " + str(defaultPlane.xpoints) + " points.\
    \n y: " + str(defaultPlane.ymin) + "m to " + str(defaultPlane.ymax) + "m with " + str(
        defaultPlane.ypoints) + " points.")

    print("\n\n\nFunction Usage:\
    \nOPTICAL PARAMETERS DEFINITION \
        \n PARAMETERS=PauLisa.Params(wavelength,w0,z0)\
    \n\nPLANE OF PROPAGATION DEFINITION \
        \n PLANE=PauLisa.Plane(xmin,xmax,xpoints,ymin,ymax,ypoints) \
    \n\nMODES DEFNITION AND DISPLAY \
        \n MODESARRAY=PauLisa.Modes((n1,m1,c1),(n2,m2,c2)) \
        \n PauLisa.ShowModes(MODES) \
    \n\nAMPLITUDE CALCULATIONS \
        \n Calculate amplitude over plane: RESULT=PauLisa.Calculate(PARAMS,PLANE,MODES,z) \
    \n Simple calculation from coordinates: PauLisa.Amplitude(PARAMS,x,y,z,MODES) \
    \n\nINTENSITY PLOTTING \
        \n PauLisa.Contour(RESULT, **xlim,**ylim) \
        \n PauLisa.IntensitySliceX(y, *RESULT, **xlim) \
        \n PauLisa.IntensitySliceY(x, *RESULT, **xlim) \
    \n\nPHASE CALCULATION \
        \n PauLisa.Phase(PARAMS,x,y,z,MODES) \
    \n\nPHASE PLOTTING \
        \n PauLisa.PhaseContour(RESULT,**xlim,**ylim) \
        \n PauLisa.PhaseSliceX(y,*RESULT,**xlim)\
        \n PauLisa.PhaseSliceY(x,*RESULT,**xlim)\
    \n\n *VARNAME represents variable number of args of specified type.")

