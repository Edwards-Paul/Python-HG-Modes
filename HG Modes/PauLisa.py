## Calculating and plotting intensity profile propagated large relative to (Zr=10, far-field 100m)
#--------------------------------------------------------------------------------------------------------
#q param. goes from i(10) -> i(10+100) 
#evaluate HG modes again at new q-param. and add together. 
#basis described by waist size (w0) and waist location into q parameter (q(z)=i*Zr+z-z0=q0+z-z0, q0=i*Zr)
#w(z) spot size - radius at which intensity is 1/e**2 max intensity I(0)
#========================================================================================================

## IMPORTS:

# Math importss
from math import pi, log, exp, sin, cos, atan, sqrt, e, factorial, radians, degrees
from scipy.special import jn
from numpy.polynomial.hermite import hermval
from scipy import integrate
from array import *
import operator
import cmath
import numpy as np
# Plot imports
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#--------------------------------------------------------------------------------------------------------
## CONSTANTS :


wavelength = 1064e-9     #Wavelength = 1064 nm
k = 2*pi /( wavelength ) #Wavenumber
W0 = 0.001 #Beam waist size, Virgo uses w0=0.02m intensity patterns for reference
z0=0
z=100
#Zr=pi*W0*W0/ wavelength 
Zr=10 #Rayleigh range (approx. length of near-field region), im. part of q-param.
q0=(1j) * Zr


#--------------------------------------------------------------------------------------------------------
## FUNCTIONS OF CHARACTERISTIC VARIABLES:


## Desc change in the radius of curvature (Rc)
def Rc(z) :
    #r=z-z0+(Zr**2/(z-z0))
    
    if z==0:
        r=float('inf')
    else :
        r=z*(1 + ( Zr / z )**2)
    return (r)



## Gouy phase
def GouyPhase(z,order) :
    PhaseLag=(order+1)*atan(z/Zr)
    Gouy=atan((z-z0)/Zr)
    return(PhaseLag)

## Spot size
def w(z):
    w=W0*sqrt(1+(z/Zr)**2)
    return(w)

## q-parameter
def q(z):
    q=q0+z-z0
    return(q)

#--------------------------------------------------------------------------------------------------------
##INITIAL INPUT:


##Get user input modes, create 2-d modes array, show array
def Modes(*argv):
    
    #get number of modes passed by user
    NumberModes = (len(argv))
    
    #create lists for n,m,c of length equal to number of modes
    listN = [None]*NumberModes
    listM = [None]*NumberModes
    listC = [None]*NumberModes
    
    #parse args into lists (n,m, & c) to set modes
    for i in range(0,len(argv)):
        listN[i],listM[i],listC[i] = argv[i]
    
    
    #get a modes 2-d array created from lists of n,m,c of required size
    return(CreateModes(listN,listM,listC,NumberModes))

##Create the modes 2-d array
def CreateModes(listN,listM,listC,NumberModes):
   
   #get max n and m to create modes grid with size based on user highest modes
    MaxN = max(listN)
    MaxM = max(listM)
    
   #get number of rows/cols for modes grid (plus 1 for 0 index..)
    rows=MaxN+1
    cols=MaxM+1
    
   #initialize 2-d array
    modes = [[0 for m in range(cols)] for n in range(rows)]
   #iterate through lists to set modes in grid
    for i in range (0,NumberModes):    
        modes[listN[i]][listM[i]]=listC[i]
   
    return(modes)
    
##Print modes
def ShowModes(modes):  
    
    
    if not modes:
        print("No modes entered.")
    
    else:
        rows = len(modes)   
        cols = len(modes[0])
        
        colList=[]

        for m in range(cols):
             colList.append(m)
        print("n\m " + str(colList))

        for n in range(rows):
            List=[]     
            for m in range(cols):
                List.append(modes[n][m])
                if m==cols-1:
                    print(str(n) +"   " + str(List))
            

#--------------------------------------------------------------------------------------------------------
##SET BASIS

#Without q
def Basis(wZero):
    global w0
    
    w0 = wZero
    print("w0 = ",w0)

#Using q-parameter
def qBasis(qZero):
    global q0
    
    q0 = qZero
    print("q0 = ",q0)

    
#--------------------------------------------------------------------------------------------------------
##PLANAR CALCULATIONS OF AMPLITUDE AND PHASE

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def Calculate(xmin,xmax,ymin,ymax,z,modes):
    
    plane=[xmin,xmax,ymin,ymax,z]
    
    x = np.arange(xmin, xmax, (xmax-xmin)/1000)
    y = np.arange(ymin, ymax, (ymax-ymin)/1000)
    X, Y = np.meshgrid(x, y)
    Z = Amplitude(X,Y,z,modes)
    return([x,y,Z,modes]) 
#return([x,y,Amplitude((X,Y),z,modes)])
    
#--------------------------------------------------------------------------------------------------------
##AMPLITUDE CALCULATIONS:

# Amplitutude calculation from w0,zR basis
def Amplitude(x,y,z,modes) :
    #Unm a sum over modes
    UnmSum=0
    rows = len(modes)   
    cols = len(modes[0])
    
    #Iterate through modes
    for n in range(rows):
        for m in range(cols):           
            carrN=[0] * rows
            carrN[n]=modes[n][m]

            carrM=[0] * cols
            carrM[m]=modes[n][m]

            order=n+m

            Unm = (2 ** (order - 1) * factorial(n) * factorial(m) * pi) ** (-1 / 2) *\
            (1 / w(z)) * e ** ((1j) * (order + 1) * GouyPhase(z,order)) *\
            e ** (-(1j) * (k * (x ** 2 + y ** 2) / (2 * Rc(z))) - ((x ** 2 + y ** 2) / (w(z) ** 2))) *\
            HermPol(n, x, carrN) *\
            HermPol(m, y, carrM)

     #Add each to result
            UnmSum+=Unm
          
    return(UnmSum)


# Amplitude from q-parameter basis
def Amplitude2(x,y,z) :

    #Unm a sum over modes
    UnmSum=0
   
    for n in range(rows):
        for m in range(cols):
            carrN=[0] * rows
            carrN[n]=modes[n][m]

            carrM=[0] * cols
            carrM[m]=modes[n][m]

            order=n+m
            
            Unm = (2/pi)**1/4 * \
            cmath.sqrt( 1 / (2**n * factorial(n) * W0) ) * \
            cmath.sqrt( q0 / q(z) ) * \
            ( ( q0 * np.conjugate(q(z)) ) / ( np.conjugate(q0) * q(z) ) )**n/2 * \
            HermPol(n, x, carrN) * \
            cmath.exp( (-( (1j) * k * x**2 )/( 2 * q(z))) ) * \
            (2/pi)**1/4 * \
            cmath.sqrt( 1 / (2**m * factorial(m) * W0) ) * \
            cmath.sqrt(q0 / q(z) ) * \
            ( ( q0 * np.conjugate(q(z) ) )/ ( np.conjugate(q0) * q(z) ) )**m/2 * \
            HermPol(m, y, carrM) * \
            cmath.exp( (-( (1j) * k * y**2 )/( 2 * q(z))) )

            UnmSum+=Unm
                    
    return((UnmSum))

#Get herm polys from modes
def HermPol(mode, coord, carr):
   
    herm = hermval(coord*sqrt(2)/w(z),carr)
    return herm

#-------------------------------------------------------------------------------------------------------
## PHASE CALCULATIONS:

def Phase(x,y,z):
    return degrees(cmath.phase(Amplitude(x,y,z)))

def Phase2(x,y,z):
    return degrees(cmath.phase(Amplitude2(x,y,z)))

#--------------------------------------------------------------------------------------------------------
## 2D INTENSITY PLOT:

def IntensitySlice(f):
    fig=plt.figure()
    x=f[0]
    #at y = 0
    amp=f[2]
    z=amp[len(amp/2)]
    plt.semilogy(x,abs(z**2))
    plt.xlabel('X')
    plt.ylabel('Intensity')
    plt.grid()
    
    #plt.savefig('IntCheck.pdf')

def IntensitySliceX(f):
    fig=plt.figure()
    x=f[0]
    #at y = 0
    amp=f[2]
    z=amp[len(amp[0])/2-1]
    plt.semilogy(x,abs(z**2))
    plt.xlabel('X')
    plt.ylabel('Intensity')
    plt.grid()
    
    #plt.savefig('IntCheck.pdf')
    

def IntensitySliceY(f):
    fig=plt.figure()
    y=f[1]
    #at y = 0
    amp=f[2]
    z=amp[len(amp)/2-1]
    plt.semilogy(y,abs(z**2))
    plt.xlabel('Y')
    plt.ylabel('Intensity')
    plt.grid()
    
    #plt.savefig('IntCheck.pdf')
    
def Contour(f):    
    h = plt.contourf(f[0],f[1],abs(f[2]**2))    
    
#--------------------------------------------------------------------------------------------------------
## 3D PLOT:

def IntensityPlot(z):
# Get data

    x = np.arange(-.03, .03, 0.001)
    y = np.arange(-.03, .03, 0.001)
    x, y = np.meshgrid(x, y)
    f = np.vectorize(Amplitude)
## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(x, y, abs(f(0,y,z))**2, rstride=1, cstride=1,
               cmap='viridis', edgecolor='none');
# Labels and render
    plt.xlabel("x")
    plt.ylabel("y");
    ax.set_zlabel("Intensity")
   # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()
        

def IntensityPlot2(z):
# Get data

    x = np.arange(-.03, .03, 0.001)
    y = np.arange(-.03, .03, 0.001)
    x, y = np.meshgrid(x, y)
    f2 = np.vectorize(Amplitude2)
## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(x, y, abs(f2(0,y,z))**2, rstride=1, cstride=1,
               cmap='viridis', edgecolor='none');
# Labels and render
    plt.xlabel("x")
    plt.ylabel("y");
    ax.set_zlabel("Intensity")
   # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()

#--------------------------------------------------------------------------------------------------------
## 2D PLOT

def AmpPlot(z):
   
# Get data
    x = np.arange(-.03, .03, 0.001)
    y = np.arange(-.03, .03, 0.001)
    x, y = np.meshgrid(x, y)
    f = np.vectorize(Amplitude)
## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(x, y, f, rstride=1, cstride=1,
               cmap='viridis', edgecolor='none');
# Labels and render
    plt.xlabel("x")
    plt.ylabel("y");
    ax.set_zlabel("Amplitude")
   # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()
    
def AmpPlot2(z):
   
# Get data
    x = np.arange(-.03, .03, 0.001)
    y = np.arange(-.03, .03, 0.001)
    x, y = np.meshgrid(x, y)
    f = np.vectorize(Amplitude)
## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(x, y, f, rstride=1, cstride=1,
               cmap='viridis', edgecolor='none');
# Labels and render
    plt.xlabel("x")
    plt.ylabel("y");
    ax.set_zlabel("Amplitude")
   # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()
    
#--------------------------------------------------------------------------------------------------------
## Phase Map
def PhaseMap(z):
   
# Get data
    x = np.arange(-.003, .003, 0.00001)
    y = np.arange(-.003, .003, 0.00001)
    x, y = np.meshgrid(x, y)
    f = np.vectorize(Phase)
## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(x, y, f(x,y,z), rstride=1, cstride=1,
               cmap='viridis', edgecolor='none');
# Labels and render
    plt.xlabel("x")
    plt.ylabel("y");
    ax.set_zlabel("deg.")
   # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()

def PhaseMap2(z):
   
# Get data
    x = np.arange(-.003, .003, 0.00001)
    y = np.arange(-.003, .003, 0.00001)
    x, y = np.meshgrid(x, y)
    f = np.vectorize(Phase2)
## Make the plot, amp. rather than int.
    fig = plt.figure()
    ax = Axes3D(fig)

    ax.plot_surface(x, y, f(x,y,z), rstride=1, cstride=1,
               cmap='viridis', edgecolor='none');
# Labels and render
    plt.xlabel("x")
    plt.ylabel("y");
    ax.set_zlabel("deg.")
   # plt.title('Intensity Profile HG Modes n,m=' + str(carrN) + "," + str(carrM))
    plt.show()    
    
#=========================================================================================================

print("Usage:\nMODESARRAY=PauLisa.Modes((n1,m1,c1),(n2,m2,c2))\nPauLisa.ShowModes(MODESARRAY)\nAMPLITUDES=PauLisa.Calculate(xmin,xmax,-ymin,ymax,z,MODESARRAY)")