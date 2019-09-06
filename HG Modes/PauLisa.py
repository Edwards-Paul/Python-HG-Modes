# Calculating and plotting intensity and phase of propagating HG modes
# q param. goes from i(10) -> i(10+100)
# evaluate HG modes again at new q-param. and add together.
# optical params described by waist size (w0) and waist location into q parameter (q(z)=i*Zr+z-z0=q0+z-z0, q0=i*Zr)
# w(z) spot size - radius at which intensity is 1/e**2 max intensity I(0)
# ========================================================================================================

## IMPORTS:

# Math imports
from math import pi, log, exp, sin, cos, atan, sqrt, e, factorial, radians, degrees
from scipy.special import jn
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
        self.z0 = z0 # waist location, I ~ 1/e**2 ! 0.13
        self.Zr = pi * w0 ** 2 / wavelength  # Rayleigh Range, near field = 2 Zr = 6
        self.q0 = (1j) * self.Zr
        self.k = 2 * pi / (wavelength)  # Wavenumber

    def __str__(self):
        return '\n{}{}\n{}{}\n{}{}\n{}{}\n{}{}'.format('wavelength=', self.wavelength ,
                                          'w0=', self.w0 ,
                                          'z0=', self.z0 ,
                                          'Zr=', self.Zr ,
                                          'q0=', self.q0 ,
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
#defaultParams = Params(1064e-9, 1e-2, 0)


# --------------------------------------------------------------------------------------------------------
## FUNCTIONS OF OPTICAL PARAMETERS:


## Desc change in the radius of curvature (Rc)
def Rc(z, params):
    # r=z-z0+(Zr**2/(z-z0))

    if z == 0:
        r = float('inf')
    else:
        r = z - params.getZ0() + ( params.getZr()**2)/(z-params.getZ0())
# r = z * (1 + (params.getZr() / z) ** 2)
    return r


## Gouy phase
def GouyPhase(z, order, params):
    Gouy = atan((z - params.getZ0()) / params.getZr())
    return Gouy

## Phase lag
def PhaseLag(z, order,params):
    #phaselag = (order + 1) * atan(z / params.getZr())
    phaselag = (order + 1) * GouyPhase(z, order, params)
    return phaselag

## Spot size (LR eq. 9.16)
def w(z, params):
    w = params.getW0() * sqrt(1 + ((z-params.getZ0()) / params.getZr()) ** 2)
    return w


## q-parameter
def q(z,params):
    q = params.getQ0() + z - params.getZ0()
    return q


# --------------------------------------------------------------------------------------------------------
##INITIAL INPUT:
# --------------------------------------------------------------------------------------------------------

##Get user input modes, create 2-d modes array, show array
def Modes(*argv):
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
    return (CreateModes(listN, listM, listC, NumberModes))


##Create the modes 2-d array
def CreateModes(listN, listM, listC, NumberModes):
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
def ShowModes(modes):
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
        return '\n{}{},{}{},{}{},{}{}\n{}{},{}{},{}{},{}{}'.format('xmin=',self.xmin,
                                                           'xmax=',self.xmax,
                                                           'xpoints=',self.xpoints,
                                                         'x step size = ', abs(self.xmax-self.xmin)/self.xpoints,
                                                           'ymin=',self.ymin,
                                                           'ymax=',self.ymax,
                                                           'ypoints=',self.ypoints,
                                                         'y step size = ', abs(self.ymax-self.ymin)/self.ypoints)

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

#x,y vectors from limits and # points
    def getX(self):
        return np.arange(self.xmin, self.xmax, abs(self.xmax-self.xmin)/self.xpoints)

    def getY(self):
        return np.arange(self.ymin, self.ymax, abs(self.ymax-self.ymin)/self.ypoints)

# --------------------------------------------------------------------------------------------------------
# DEFAULT PLANE OF CALCULATION
#defaultPlane = Plane(-0.05, 0.05, 1000, -0.05, 0.05, 1000)
defaultPlane = Plane(-2e-2, 2e-2, 1000, -2e-2, 2e-2, 1000)
#defaultPlane = Plane(-5e-6,5e-6,1000,-5e-6,5e-6,1000)
#defaultPlane = Plane(-.2, .2, 1000, -.2, .2, 1000)
# --------------------------------------------------------------------------------------------------------
##PLANAR CALCULATIONS OF AMPLITUDE AND PHASE

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def Calculate(params, plane, modes, z):

    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = Amplitude(params, X, Y, z, modes)
    phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

# --------------------------------------------------------------------------------------------------------
# Calculate peak and its location from result
class peakInt:

    def __init__(self,result):
        intensity = abs(result.amp)**2

        map(max, intensity)
        list(map(max, intensity))

        self.peak = max(map(max, intensity))
        self.loc = np.where(intensity == self.peak)

        self.x = self.loc[1] * abs((result.plane.xmax - result.plane.xmin) / result.plane.xpoints)+result.plane.xmin
        self.y = self.loc[0] * abs((result.plane.ymax - result.plane.ymin) / result.plane.ypoints)+result.plane.ymin

        #print ("x=" + str(self.x) + "y=" + str(self.y) + "peak="+str(self.peak))



class peakAmp:

    def __init__(self,result):
        amplitude = (result.amp)

        map(max, amplitude)
        list(map(max, amplitude))

        self.peak = max(map(max, amplitude))
        self.loc = np.where(amplitude == self.peak)

        self.x = self.loc[1] * abs((result.plane.xmax - result.plane.xmin) / result.plane.xpoints)+result.plane.xmin
        self.y = self.loc[0] * abs((result.plane.ymax - result.plane.ymin) / result.plane.ypoints)+result.plane.ymin

        print ("x=" + str(self.x) + "y=" + str(self.y) + "peak="+str(self.peak))

# --------------------------------------------------------------------------------------------------------
class Result:

    def __init__(self, params, plane, modes, z, amp, phase):
        self.params = params
        self.plane = plane
        self.z = z
        self.amp = amp
        self.modes = modes
        self.phase = phase

    def __str__(self):
        return '{}{}\n\n{}{}\n\n{}{}\n\n{}{}\n\n{}{}\n\n{}{}'.format('PARAMS: ',self.params,
                                                             'PLANE: ',self.plane,
                                                             'MODES: ',self.modes,
                                                             'Z: ',self.z,
                                                             'AMP: ', self.amp,
                                                             'PHASE: ', self.phase)

    def getParams(self):
        return self.params

    def getPlane(self):
        return self.plane

    def getZ(self):
        return self.z

    def getAmp(self):
        return self.amp

    def getPhase(self):
        return self.phase

    def getModes(self):
        return self.modes


# --------------------------------------------------------------------------------------------------------
##AMPLITUDE CALCULATIONS:

# Amplitutude calculation from w0,Zr. (LR eq. 9.26)
def Amplitude(params, x, y, z, modes):
    # Unm a sum over modes
    UnmSum = 0
    rows = len(modes)
    cols = len(modes[0])

    # Iterate through modes
    for n in range(rows):
        for m in range(cols):
            # n array for hermpol
            carrN = [0] * rows
            carrN[n] = modes[n][m]
            #print(carrN)

            # m array for hermpol
            carrM = [0] * cols
            carrM[m] = modes[n][m]
            #print(carrM)

            order = n + m

            Unm = (2 ** (order - 1) * factorial(n) * factorial(m) * pi) ** (-1 / 2) * \
                (1 / w(z, params)) * e ** ((1j) * (order + 1) * GouyPhase(z, order, params)) * \
                  e ** (-(1j) * ( (params.getK() * (x ** 2 + y ** 2)) / (2 * Rc(z, params))) - \
                  ((x ** 2 + y ** 2) / (w(z, params) ** 2))) * \
                  HermPol(n, x, carrN, z, params) * \
                  HermPol(m, y, carrM, z, params)

            # Add each to result
            UnmSum += Unm
            #print(Unm,UnmSum)

    return (UnmSum)


## Get herm polys from modes
# mode: working mode (as Calculate iterates through n,m grid
# coord: x or y
# carr: coefficient array
def HermPol(mode, coord, carr, z, params):
    herm = hermval(coord * sqrt(2) / w(z, params), carr)
    return herm


# --------------------------------------------------------------------------------------------------------
## 2D INTENSITY PLOT:

# These IntensitySlice's require recalculation at x and y plane.
# Plotting calc at, e.g., halfway points in x-y grid for x=0 (center col.) or y=0 (center row) accomplishes the same.
def AmplitudeSliceX(y, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), amp, label = i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break
    plt.title('Amp. along x')
    plt.xlabel('X (m)')
    plt.ylabel('Amplitude')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')

def AmplitudeSliceY(x, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), np.real(amp), label= i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)

    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('amp along y')
    plt.xlabel('Y (m)')
    plt.ylabel('amp')
    plt.legend(loc='upper right')
    plt.grid()
    # plt.savefig('IntCheck.pdf')

def IntensitySliceX(y, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), abs(amp )** 2, label = i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break
    plt.title('Intensity along x')
    plt.xlabel('X (m)')
    plt.ylabel('Intensity')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def IntensitySliceY(x, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), abs(amp ** 2), label= i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)

    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('Intensity along y')
    plt.xlabel('Y (m)')
    plt.ylabel('Intensity')
    plt.legend(loc='upper right')
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def Contour(f, **kwargs):
    fig, ax = plt.subplots(figsize=(7,6))
    cs = plt.contourf(f.plane.getX(), f.plane.getY(), abs(f.getAmp() ** 2))
    # cs = plt.contourf(f[0],f[1],abs(f[3]**2), locator=matplotlib.ticker.LogLocator())
    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])
        # optionally set limits. default is last result's range
    if ('ylim' in kwargs):
        plt.ylim(kwargs['ylim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    ax.yaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), mathText=True)

    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    for i in plt.ylim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-6, -6),mathText=True)
            break


    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    cbar = fig.colorbar(cs)
    plt.title('Intensity')

# formatting axes
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, mathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self,useOffset=offset,useMathText=mathText)
    def _set_orderOfMagnitude(self, nothing):
        self.orderOfMagnitude = self.oom
    def _set_format(self, vmin, vmax):
        self.format = self.fformat
        if self._useMathText:
            self.format = '$%s$' % matplotlib.ticker._mathdefault(self.format)
# --------------------------------------------------------------------------------------------------------
## 3D PLOT:
#
def IntensityPlot(f):
    # Get data

    ## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(f.plane.getX(), f.plane.getY(), abs(f.getAmp()) ** 2, rstride=1, cstride=1,
                    cmap='viridis', edgecolor='none')
    # Labels and render
    plt.xlabel("x")
    plt.ylabel("y")
    ax.set_zlabel("Intensity")
    # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()
#
#
# def IntensityPlot2(z):
#     # Get data
#
#     x = np.arange(-.03, .03, 0.001)
#     y = np.arange(-.03, .03, 0.001)
#     x, y = np.meshgrid(x, y)
#     f2 = np.vectorize(Amplitude2)
#     ## Make the plot, amp. rather than int.
#     fig = plt.figure()
#     ax = Axes3D(fig)
#
#     ax.plot_surface(x, y, abs(f2(0, y, z)) ** 2, rstride=1, cstride=1,
#                     cmap='viridis', edgecolor='none')
#     # Labels and render
#     plt.xlabel("x")
#     plt.ylabel("y")
#     ax.set_zlabel("Intensity")
#     # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
#     plt.show()
#
#
# --------------------------------------------------------------------------------------------------------
## PHASE CALCULATIONS:

def Phase(params, x, y, z, modes):
    return (np.angle(Amplitude(params, x, y, z, modes),deg = True))


# --------------------------------------------------------------------------------------------------------
## Phase Map
# def PhaseMap(f):
#     x = f.getX()
#     y = f.getY()
#     z = f.getZ()
#     modes = f.getModes()
#     x, y = np.meshgrid(x, y)
#     ph = np.vectorize(Phase)
#
#     fig = plt.figure()
#     ax = Axes3D(fig)
#     ax.plot_surface(x, y, ph(x, y, z), rstride=1, cstride=1,
#                     cmap='viridis', edgecolor='none');
#     # Labels and render
#     plt.xlabel("x")
#     plt.ylabel("y");
#     ax.set_zlabel("deg.")
#     # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
#     plt.show()
#
#
# def PhaseMap2(z):
#     # Get data
#     x = np.arange(-.003, .003, 0.00001)
#     y = np.arange(-.003, .003, 0.00001)
#     x, y = np.meshgrid(x, y)
#     f = np.vectorize(Phase2)
#     ## Make the plot, amp. rather than int.
#     fig = plt.figure()
#     ax = Axes3D(fig)
#
#     ax.plot_surface(x, y, f(x, y, z), rstride=1, cstride=1,
#                     cmap='viridis', edgecolor='none');
#     # Labels and render
#     plt.xlabel("x")
#     plt.ylabel("y");
#     ax.set_zlabel("deg.")
#     # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
#     plt.show()


def PhaseContour(f,**kwargs):
    fig, ax = plt.subplots()
    cs = plt.contourf(f.plane.getX(), f.plane.getY(), f.getPhase())
    # cs = plt.contourf(f[0],f[1],abs(f[3]**2), locator=matplotlib.ticker.LogLocator())
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])
        # optionally set limits. default is last result's range
    if ('ylim' in kwargs):
        plt.ylim(kwargs['ylim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    ax.yaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), mathText=True)
    #scale x
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break
    #scale y
    for i in plt.ylim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    cbar = fig.colorbar(cs)
    plt.title('Phase')


def PhaseSliceX(y, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # calc at from z and modes in f (*argv)

    for i in range(0, len(argv)):
        phase = Phase(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), phase, label=i+1)
    # optionally, set limits. default is last result
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    #scale x
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('Phase along x')
    plt.xlabel('X')
    plt.ylabel('deg.')
    plt.legend(loc='upper right')

    plt.grid()

def PhaseSliceY(x, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Calc phase from z and modes in f

    for i in range(0, len(argv)):
        phase = Phase(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), phase, label=i+1)

    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    #scale x
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('Phase along y')
    plt.xlabel('Y (m)')
    plt.ylabel('deg.')
    plt.legend(loc='upper right')

    plt.grid()
    # plt.savefig('IntCheck.pdf')

# =========================================================================================================
## PRINTING DEFAULTS AND USAGE

print("DEFAULT PARAMS (PauLisa.defaultParams)\
\n wavelength =" + str(defaultParams.wavelength) + "\
m\n waist size(w0) =" + str(defaultParams.w0) + "\
m\n z0 ="+ str(defaultParams.z0)+ "\
m\n Rayleigh Range (Zr) ="+ str(defaultParams.Zr)+ "m")

print("\n\nDEFAULT X,Y PLANE (PauLisa.defaultPlane)\
\n x: "+ str(defaultPlane.xmin)+ "m to "+ str(defaultPlane.xmax)+ "m with "+ str(defaultPlane.xpoints)+ " points.\
\n y: "+ str(defaultPlane.ymin)+ "m to "+ str(defaultPlane.ymax)+ "m with "+ str(defaultPlane.ypoints)+ " points.")

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



## from q param

# Amplitude from q-parameter (L.R., Eq. 9.34)
def Amplitude2(params, x, y, z, modes):
    # Unm a sum over modes
    UnmSum = 0
    rows = len(modes)
    cols = len(modes[0])

    for n in range(rows):
        for m in range(cols):
            carrN = [0] * rows
            carrN[n] = modes[n][m]

            #print carrN

            carrM = [0] * cols
            carrM[m] = modes[n][m]
            #print carrM
            order = n + m

            #print n,m

            Unm = (2 / pi) ** 1 / 4 * \
                  np.sqrt(1 / (2 ** n * factorial(n) * params.getW0())) * \
                  np.sqrt(params.getQ0() / q(z,params)) * \
                  ((params.getQ0() * np.conjugate(q(z,params))) / (np.conjugate(params.getQ0()) * q(z,params))) ** n / 2 * \
                  HermPol(n, x, carrN, z, params) * \
                  np.exp((-((1j) * params.getK() * x ** 2) / (2 * q(z,params)))) * \
                  (2 / pi) ** 1 / 4 * \
                  np.sqrt(1 / (2 ** m * factorial(m) * params.getW0())) * \
                  np.sqrt(params.getQ0() / q(z,params)) * \
                  ((params.getQ0() * np.conjugate(q(z,params))) / (np.conjugate(params.getQ0()) * q(z,params))) ** m / 2 * \
                  HermPol(m, y, carrM, z, params) * \
                  np.exp((-((1j) * params.getK() * y ** 2) / (2 * q(z,params))))

            UnmSum += Unm
            #print(Unm, UnmSum)

    return (UnmSum)

def Calculate2(params, plane, modes, z):

    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = Amplitude2(params, X, Y, z, modes)
    phase = Phase2(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

def AmplitudeSliceX2(y, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude2(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), amp, label = i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break
    plt.title('Amp2 along x')
    plt.xlabel('X (m)')
    plt.ylabel('amp')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')

def IntensitySliceX2(y, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude2(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), abs(amp ** 2), label = i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break
    plt.title('Intensity along x')
    plt.xlabel('X (m)')
    plt.ylabel('Intensity')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def IntensitySliceY2(x, *argv, **kwargs):

    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = Amplitude2(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), abs(amp ** 2), label= i+1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)

    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('Intensity along y')
    plt.xlabel('Y (m)')
    plt.ylabel('Intensity')
    plt.legend(loc='upper right')
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def Phase2(params, x, y, z, modes):
    return (np.angle(Amplitude2(params, x, y, z, modes), deg=True))
    #return degrees(cmath.phase(Amplitude2(params, x, y, z, modes)))

def PhaseSliceX2(y, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # calc at from z and modes in f (*argv)

    for i in range(0, len(argv)):
        phase = Phase2(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), phase, label=i+1)
    # optionally, set limits. default is last result
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    #scale x
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('Phase along x')
    plt.xlabel('X')
    plt.ylabel('deg.')
    plt.legend(loc='upper right')

    plt.grid()

def PhaseSliceY2(x, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Calc phase from z and modes in f

    for i in range(0, len(argv)):
        phase = Phase2(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), phase, label=i+1)

    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), mathText=True)
    #scale x
    for i in plt.xlim():
        if  1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), mathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6),mathText=True)
            break

    plt.title('Phase along y')
    plt.xlabel('Y (m)')
    plt.ylabel('deg.')
    plt.legend(loc='upper right')

    plt.grid()
    # plt.savefig('IntCheck.pdf')