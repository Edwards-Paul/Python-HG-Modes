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
from PauLisa import *

# --------------------------------------------------------------------------------------------------------
## 2D INTENSITY PLOT:

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

# These IntensitySlice's require recalculation at x and y plane.
# Plotting calc at, e.g., halfway points in x-y grid for x=0 (center col.) or y=0 (center row) accomplishes the same.
def ampslicex(y, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    if('labels' in kwargs):
        for i in range(0, len(argv)):
            amp = amplitude(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), amp, label=kwargs['labels'][i])
    else:
        for i in range(0, len(argv)):
            amp = amplitude(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), amp, label=i + 1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break
    plt.title('Amp. along x')
    plt.xlabel('X (m)')
    plt.ylabel('Amplitude')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def ampslicey(x, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = amplitude(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), np.real(amp), label=i + 1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)

    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('amp along y')
    plt.xlabel('Y (m)')
    plt.ylabel('amp')
    plt.legend(loc='upper right')
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def intslicex(y, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f

    if('labels' in kwargs):
        for i in range(0, len(argv)):
            amp = amplitude(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), abs(amp) ** 2, label=kwargs['labels'][i])
    else:
        for i in range(0, len(argv)):
            amp = amplitude(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), abs(amp) ** 2, label=i + 1)
    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    if ('ylim' in kwargs):
        plt.ylim(kwargs['ylim'])

#     ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
#     # format x axis to mm or microns depending on order of x limits
#     for i in plt.xlim():
#         if 1e-5 < abs(i) < 1e-2:
#             ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.3f"))
#             ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
#         if abs(i) <= 1e-5:
#             ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.3f"))
#             ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
#             break
    plt.title('Intensity along X ')
    plt.xlabel('X [m]')
    plt.ylabel('Intensity [a.u.]')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def intslicey(x, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = amplitude(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), abs(amp ** 2), label=i + 1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)

    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('Intensity along y')
    plt.xlabel('Y (m)')
    plt.ylabel('Intensity')
    plt.legend(loc='upper right')
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def Contour(f, **kwargs):
    fig, ax = plt.subplots(figsize=(7, 6))
    cs = plt.contourf(f.plane.getX(), f.plane.getY(), abs(f.getAmp() ** 2))
    # cs = plt.contourf(f[0],f[1],abs(f[3]**2), locator=matplotlib.ticker.LogLocator())
    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])
        # optionally set limits. default is last result's range
    if ('ylim' in kwargs):
        plt.ylim(kwargs['ylim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    ax.yaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)

    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.3f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.3f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    for i in plt.ylim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.3f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.3f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    cbar = fig.colorbar(cs)
    plt.title('Intensity')


# formatting axes
class OOMFormatter(matplotlib.ticker.ScalarFormatter):
    def __init__(self, order=0, fformat="%1.1f", offset=True, useMathText=True):
        self.oom = order
        self.fformat = fformat
        matplotlib.ticker.ScalarFormatter.__init__(self, useOffset=offset, useMathText = True)

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


def PhaseContour(f, **kwargs):
    fig, ax = plt.subplots()
    cs = plt.contourf(f.plane.getX(), f.plane.getY(), f.getPhase())
    # cs = plt.contourf(f[0],f[1],abs(f[3]**2), locator=matplotlib.ticker.LogLocator())
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])
        # optionally set limits. default is last result's range
    if ('ylim' in kwargs):
        plt.ylim(kwargs['ylim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    ax.yaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
    # scale x
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break
    # scale y
    for i in plt.ylim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    cbar = fig.colorbar(cs)
    plt.title('Phase')


def phaseslicex(y, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # calc at from z and modes in f (*argv)
    if('labels' in kwargs):
        for i in range(0, len(argv)):
            angle = phase(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), angle, label=kwargs['labels'][i])
    else:
        for i in range(0, len(argv)):
            angle = phase(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), angle, label=i + 1)

    # optionally, set limits. default is last result
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # scale x
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.3f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.3f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('Phase along x')
    plt.xlabel('X [m]')
    plt.ylabel('Angle [rad]')
    plt.legend(loc='upper right')

    plt.grid()


def phaseslicey(x, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Calc phase from z and modes in f

    for i in range(0, len(argv)):
        angle = phase(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), angle, label=i + 1)

    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # scale x
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('Phase along y')
    plt.xlabel('Y [m]')
    plt.ylabel('Angle [rad]')
    plt.legend(loc='upper right')

    plt.grid()
    # plt.savefig('IntCheck.pdf')
    
def ampslicex_q(y, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = amplitude_q(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getX(), amp, label=i + 1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break
    plt.title('Amp2 along x')
    plt.xlabel('X [m]')
    plt.ylabel('Amplitude [a.u.]')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def intslicex_q(y, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    if('labels' in kwargs):
        for i in range(0, len(argv)):
            amp = amplitude_q(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), abs(amp) ** 2, label=kwargs['labels'][i])
    else:
        for i in range(0, len(argv)):
            amp = amplitude_q(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), abs(amp) ** 2, label=i + 1)
    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    if ('ylim' in kwargs):
        plt.ylim(kwargs['ylim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # format x axis to mm or microns depending on order of x limits
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break
    plt.title('Intensity along x')
    plt.xlabel('X [m]')
    plt.ylabel('Intensity [a.u]')
    plt.legend(loc='upper right')
    # ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def intslicey_q(x, *argv, **kwargs):
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)

    # Calc amp from z and modes in f
    for i in range(0, len(argv)):
        amp = amplitude_q(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), abs(amp ** 2), label=i + 1)

    # optionally set limits. default is last result's range
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)

    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('Intensity along y')
    plt.xlabel('Y [m]')
    plt.ylabel('Intensity [a.u]')
    plt.legend(loc='upper right')
    plt.grid()
    # plt.savefig('IntCheck.pdf')


def phase_q(params, x, y, z, modes):
    return (np.angle(amplitude_q(params, x, y, z, modes)))
    # return degrees(cmath.phase(Amplitude2(params, x, y, z, modes)))


def phaseslicex_q(y, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # calc at from z and modes in f (*argv)

    if('labels' in kwargs):
        for i in range(0, len(argv)):
            angle = phase_q(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), angle, label=kwargs['labels'][i])
    else:
        for i in range(0, len(argv)):
            angle = phase_q(argv[i].getParams(), argv[i].plane.getX(), y, argv[i].getZ(), argv[i].getModes())
            plt.plot(argv[i].plane.getX(), angle, label=i + 1)
    # optionally, set limits. default is last result
    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.3f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # scale x
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('Phase along x')
    plt.xlabel('X [m]')
    plt.ylabel('Angle [rad]')
    plt.legend(loc='upper right')

    plt.grid()


def phaseslicey_q(x, *argv, **kwargs):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # Calc phase from z and modes in f

    for i in range(0, len(argv)):
        angle = phase_q(argv[i].getParams(), x, argv[i].plane.getY(), argv[i].getZ(), argv[i].getModes())
        plt.plot(argv[i].plane.getY(), angle, label=i + 1)

    if ('xlim' in kwargs):
        plt.xlim(kwargs['xlim'])

    ax.xaxis.set_major_formatter(OOMFormatter(0, "%1.2f"))
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
    # scale x
    for i in plt.xlim():
        if 1e-5 < abs(i) < 1e-2:
            ax.xaxis.set_major_formatter(OOMFormatter(-3, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-3, -3), useMathText=True)
        if abs(i) <= 1e-5:
            ax.xaxis.set_major_formatter(OOMFormatter(-6, "%1.2f"))
            ax.ticklabel_format(axis='x', style='sci', scilimits=(-6, -6), useMathText=True)
            break

    plt.title('Phase along y')
    plt.xlabel('Y [m]')
    plt.ylabel('Angle [rad]')
    plt.legend(loc='upper right')

    plt.grid()
    # plt.savefig('IntCheck.pdf')