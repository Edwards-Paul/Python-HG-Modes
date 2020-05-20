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

    if (z == params.getZ0() ):
        r = float('inf')
    else:
        r = z - params.getZ0() + (params.getZr() ** 2) / (z - params.getZ0())
    # r = z * (1 + (params.getZr() / z) ** 2)

    return r


## Gouy phase
def gouy_phase(z, params):
    gouy = atan((z - params.getZ0()) / params.getZr())
    return gouy


## Phase lag
def phase_lag(z, order, params):
    # phaselag = (order + 1) * atan(z / params.getZr())
    phase = (order + 1) * gouy_phase(z, order, params)
    return phase


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

##Get user input modes, from array of 3-element arrays
def modes(modes_array):
    # get number of modes passed by user
    
    NumberModes = (len(modes_array))
    # create lists for n,m,c of length equal to number of modes
    listN = [None] * NumberModes
    listM = [None] * NumberModes
    listC = [None] * NumberModes

    # parse args into lists (n,m, & c) to set modes
    for i in range(0, len(modes_array)):      
        listN[i], listM[i], listC[i] = modes_array[i][0],modes_array[i][1],modes_array[i][2]

    # get a modes 2-d array created from lists of n,m,c of required size
    return (create_modes(listN, listM, listC, NumberModes))

def return_array(arr):
    arr1 = arr
    return(arr)

##Create modes from 1D array (as in sim)############################################
def modes_from_1d_array(modes_array):

    NumberModes = int(len(modes_array))
    listModesN = [None] * NumberModes
    listModesM = [None] * NumberModes
    listModesC = [None] * NumberModes
    
    for i in range(NumberModes):
        A=i
        N= N_f(A)
        listModesN[i] = int(m(N,A))
        listModesM[i] = int(n(N,A))
        listModesC[i] = modes_array[i]
    
    modes_list =create_modes(listModesN, listModesM, listModesC, NumberModes) 
    return (modes_array)

##Create modes from 1D array (as in sim)############################################
def modes_from_1d_array_LV(modes_array_real,modes_array_imag):

    NumberModes = int(len(modes_array_real))
    listModesN = [None] * NumberModes
    listModesM = [None] * NumberModes
    listModesC = [None] * NumberModes
    
    for i in range(NumberModes):
        A=i
        N= N_f(A)
        listModesN[i] = int(m(N,A))
        listModesM[i] = int(n(N,A))
        C_real = modes_array_real[i]
        C_imag = modes_array_imag[i]
        listModesC[i] = C_real+ (1)*C_imag
    
    modes_list =create_modes(listModesN, listModesM, listModesC, NumberModes) 
    return (np.array(modes_list))

def n(N,A):
    res = (N+1)*(N+2)/2 - (A+1)
    return(res)

def m(N,A):
    res = A - (N*(N+1)/2)
    return(res)

def N_f(A):   
    res = np.floor((np.sqrt(8*A+1)-1)/2)  
    return(res)


##################################################################################

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

##Create the modes 2-d array
def create_modes_order18(listN, listM, listC, NumberModes):
    # get max n and m to create modes grid with size based on user highest modes
    MaxN = max(listN)
    MaxM = max(listM)

    # get number of rows/cols for modes grid (plus 1 for 0 index..)
    rows = 18 + 1
    cols = 18 + 1

    # initialize 2-d array
    modes = [[0 for m in range(cols)] for n in range(rows)]
    
    # iterate through lists to set modes in grid
    for i in range(0, NumberModes):
        if (listN[i] + listM[i]) < 19:
            modes[listN[i]][listM[i]] = listC[i]

    return (modes)

##Create the modes 2-d array up to mode order N from tophat coeff.
def create_modes_orderN(listN, listM, listC, NumberModes,N):
    # get max n and m to create modes grid with size based on user highest modes
    MaxN = max(listN)
    MaxM = max(listM)

    # get number of rows/cols for modes grid (plus 1 for 0 index..)
    rows = N + 1
    cols = N + 1

    # initialize 2-d array
    modes = [[0 for m in range(cols)] for n in range(rows)]
    
    # iterate through lists to set modes in grid
    for i in range(0, NumberModes):
        if (listN[i] + listM[i]) < (N+1):
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

# where alpha an array
def iterative_scatter_case2(z,params,a,alpha,modes):
    modes_list = []
    
    for i in range(len(alpha)):
        mode = scatter_case2(z,params,a,alpha[i],modes)
        modes_list.append(mode)
    return(modes_list)

    #eta from doc, but must implement (n+m+1). def separately because mode dep.
# def eta(n,m,z,z_r)
#     res = (n+m+1)*(z_r/(z**2+z_r**2))
#     return(res)

# Apply general mode scatter formulas for x-dependence (uni-directional offset coord.), RETURN modes
def scatter_case2_first_order(z,params,a,alpha,modes):
    
    
    k = params.getK()
    w0 = params.getW0()
    zr = params.getZr()
    w_z = w(z,params)
    gouy = gouy_phase(z,params) 
    r_c = radius_curvature(z, params)
    phi = e**((1j)*gouy)
    
    #iterate through modes
    rows = len(modes)
    cols = len(modes[0])
    
    number_modes = rows*cols
    
    #build new modes (up to 2 orders larger for quad. dep.)
    new_modes = [[0 for m in range(cols+2)] for n in range(rows+2)]
    

    
    # Add u(n+2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                sc_plus_2 = x_plus_2(w0,z,zr,n)       #scattering factor from x terms (factor with off_plus_2)        
                
                                #eta from doc, but must implement (n+m+1)
                eta = (n+m+1)*(zr/(z**2+zr**2))

                # this from x^2*u_{n,m}
                off_2 = ( -(1j)*2*k*a*alpha/(w_z**2) 
                           +(1j)*a*alpha*eta/(w_z**2)       
                           + (a*alpha*2*k**2/r_c)
                           - (a*alpha*k*eta/ r_c)              
                ) 
                
                #Add at n+2
                new_modes[n+2][m] += (c_nm*(sc_plus_2*off_2))
                
    # Add u(n+1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                sc_plus_1 = x_plus_1(w0,z,zr,n)       #scattering factor from x terms        

                
                    #eta from doc, but must implement (n+m+1)
                eta = (n+m+1)*(zr/(z**2+zr**2))


                # this from x*u_{n,m}
                off_1 = (
                            (1j)*k*alpha
                           - (1j) *alpha * eta
                           - (2*a/w_z**2)
                           -(a*(1j)*k/r_c)
                )
                
                #Add at n+1
                new_modes[n+1][m] += ( c_nm*(sc_plus_1*off_1) )
                
    # Add u(n,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                
                sc_plus_1 = x_plus_1(w0,z,zr,n)       #scattering factor from x terms        
               
                sc_zero_2 = x_zero_2(w0,z,zr,n)
                
                #eta from doc, but must implement (n+m+1)
                eta = (n+m+1)*(zr/(z**2+zr**2))

                # this from x^2*u_{n,m}
                off_2 = ( -(1j)*2*k*a*alpha/(w_z**2) 
                           +(1j)*a*alpha*eta/(w_z**2)       
                           + (a*alpha*2*k**2/r_c)
                           - (a*alpha*k*eta/ r_c)              
                ) 
                
                                # this from x*u_{n-1,m}
                off_1_nminus1 = (
                        ( (1j) *k*alpha*2*a*sqrt(n) * phi/w_z )
                        -
                        ( (1j) * a*alpha * eta * 2 *sqrt(n) * phi / w_z)

                )
                
                #Add at n
                new_modes[n][m] += ( c_nm*(1 + (sc_plus_1*off_1_nminus1) + (sc_zero_2*off_2) ) )
                
    # Add u(n-1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>0) ): #neglect n-1
                c_nm = modes[n][m]
                
                sc_minus_1 = x_minus_1(w0,z,zr,n)       #scattering factor from x terms        
                
                #from herm poly
                sc = ( 2*a*sqrt(n)*e**( (1j)*gouy)/w_z )
                
                    #eta from doc, but must implement (n+m+1)
                eta = (n+m+1)*(zr/(z**2+zr**2))


                # this from x*u_{n,m}
                off_1 = (
                            (1j)*k*alpha
                           - (1j) *alpha * eta
                           - (2*a/w_z**2)
                           -(a*(1j)*k/r_c)
                )
                
                #Add at n
                new_modes[n-1][m] += ( c_nm*(sc + (sc_minus_1*off_1) ) )
           
                
    # Add u(n-2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>1) ): #neglect n=0,1
                c_nm = modes[n][m]
                
                sc_minus_1 = x_minus_1(w0,z,zr,n)       #scattering factor from x terms        

                sc_minus_2 = x_minus_2(w0,z,zr,n)       #scattering factor from x terms        
                
                    #eta from doc, but must implement (n+m+1)
                eta = (n+m+1)*(zr/(z**2+zr**2))

                # this from x^2*u_{n,m}
                off_2 = ( -(1j)*2*k*a*alpha/(w_z**2) 
                           +(1j)*a*alpha*eta/(w_z**2)       
                           + (a*alpha*2*k**2/r_c)
                           - (a*alpha*k*eta/ r_c)              
                ) 
                
                                # this from x*u_{n-1,m}
                off_1_nminus1 = (
                        ( (1j) *k*alpha*2*a*sqrt(n) * phi/w_z )
                        -
                        ( (1j) * a*alpha * eta * 2 *sqrt(n) * phi / w_z)

                )
                
                #Add at n
                new_modes[n-2][m] += ( c_nm*( (sc_minus_1*off_1_nminus1) + (sc_minus_2*off_2) ) )
    
    #account for symmetry, 0,m_final = m_final,0
    if (rows%2==0 and cols%2==0 and rows==cols):
        new_modes[0][rows+2]=new_modes[rows+2][0]
    
    return(new_modes)

# Apply general mode scatter formulas for x-dependence (uni-directional offset coord.), RETURN modes
def scatter_first_order(z,params,a,alpha,modes):
    
    #on rotation
    #z_sub=z
    #z = z*np.cos(alpha)-(a)*alpha
    #x = (a)+z_sub*alpha
    #z = z*np.cos(alpha)-(x_sub+a)*alpha
    
    k = params.getK()
    w0 = params.getW0()
    zr = params.getZr()
    w_z = w(z,params)
    gouy = gouy_phase(z,params) 
    
    #iterate through modes
    rows = len(modes)
    cols = len(modes[0])
    
    number_modes = rows*cols
    
    #build new modes (up to 2 orders larger for quad. dep.)
    new_modes = [[0 for m in range(cols+2)] for n in range(rows+2)]
    

    
    # Add u(n+2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                sc_plus_2 = x_plus_2(w0,z,zr,n)       #scattering factor from x terms (factor with off_plus_2)        
                off_plus_2 = -( (1j)*2*k*a*alpha/(w_z**2) )  #factor from shift/tilt offset terms(factor with sc_plus_2)
                
                #Add at n+2
                new_modes[n+2][m] += (c_nm*(sc_plus_2*off_plus_2) )
                
    # Add u(n+1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                sc_plus_1 = x_plus_1(w0,z,zr,n)       #scattering factor from x terms        
                off_plus_1 = ( (1j)*k*alpha - (2*a)/w_z**2 )  #factor from shift/tilt offset terms
                
                #Add at n+1
                new_modes[n+1][m] += ( c_nm*(sc_plus_1*off_plus_1) )

    # Add u(n,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                
                sc_plus_1 = x_plus_1(w0,z,zr,n)       #scattering factor from x terms        
                off_plus_1 = ( (1j)*k*alpha*(2*a*sqrt(n))/w_z * e**((1j)*gouy) )  #factor from shift/tilt offset terms
                
                sc_zero_2 = x_zero_2(w0,z,zr,n)
                off_zero_2 = -( (1j)*2*k*a*alpha/w_z**2 )
                
                new_terms = -(2*a*alpha*z/w_z**2 - (1j)*k*z*alpha**2*a/w_z**2)
                
                #Add at n
                new_modes[n][m] += ( c_nm*(1 + new_terms + (sc_plus_1*off_plus_1) + (sc_zero_2*off_zero_2) ) )
                
    # Add u(n-1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>0) ): #neglect n-1
                c_nm = modes[n][m]
                
                sc_minus_1 = x_minus_1(w0,z,zr,n)       #scattering factor from x terms        
                off_minus_1 = ( (1j)*k*alpha - 2*a/w_z**2)  #factor from shift/tilt offset terms

                sc = ( 2*a*sqrt(n)*e**( (1j)*gouy)/w_z )
                
                #Add at n
                new_modes[n-1][m] += ( c_nm*(sc + (sc_minus_1*off_minus_1) ) )
               
                
    # Add u(n-2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>1) ): #neglect n=0,1
                c_nm = modes[n][m]
                
                sc_minus_1 = x_minus_1(w0,z,zr,n)       #scattering factor from x terms        
                off_minus_1 = ( (1j)*k*alpha*a*2*sqrt(n)*e**((1j)*gouy)/w_z )  #factor from shift/tilt offset terms

                sc_minus_2 = x_minus_2(w0,z,zr,n)       #scattering factor from x terms        
                off_minus_2 = -( (1j)*2*k*a*alpha/w_z**2 )  #factor from shift/tilt offset terms
                
                #Add at n
                new_modes[n-2][m] += ( c_nm*( (sc_minus_1*off_minus_1) + (sc_minus_2*off_minus_2) ) )
    
    #account for symmetry, 0,m_final = m_final,0
    if (rows%2==0 and cols%2==0 and rows==cols):
        new_modes[0][rows+2]=new_modes[rows+2][0]
    
    return(new_modes)

# --------------------------------------------------------------------------------------------------------
#x dep
def x_plus_1(w0,z,zr,n):    
    factor = (w0/2)*( ( 1-(1j)*(z/zr) )*np.sqrt(n+1))
    return(factor)

#x dep
def x_minus_1(w0,z,zr,n):
    factor = (w0/2)*( np.sqrt(n)*(1+(1j)*(z/zr)) )
    return(factor)

#x**2 dep
def x_plus_2(w0,z,zr,n):
    factor = (w0/2)**2 * ( (1-(1j)*z/zr)**2 * np.sqrt((n+1)*(n+2)) )
    return(factor)

#x**2 dep
def x_zero_2(w0,z,zr,n):
    factor = (w0/2)**2 *( (2*n+1) * (1+(z/zr)**2 ) )
    return(factor)

#x**2 dep
def x_minus_2(w0,z,zr,n):
    factor = (w0/2)**2 * ( np.sqrt(n*(n-1)) * (1+(1j)*(z/zr))**2 )
    return(factor)

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
# defaultPlane = Plane(-0.05, 0.05, 1000, -0.05, 0.05, 1000)
defaultPlane = Plane(-2e-2, 2e-2, 1000, -2e-2, 2e-2, 1000)


# defaultPlane = Plane(-5e-6,5e-6,1000,-5e-6,5e-6,1000)
# defaultPlane = Plane(-.2, .2, 1000, -.2, .2, 1000)
# --------------------------------------------------------------------------------------------------------
##PLANAR CALCULATIONS OF AMPLITUDE AND PHASE

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def calculate(params, plane, modes, z):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude(params, X, Y, z, modes)
    
    #get phase from amp.
    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def calculate_shift_tilt(params, plane, modes, z,a,alpha):
    #x and y are arrays (xmin,xmax)
    X, Y = plane.getX(), plane.getY()
    
    rows = len(X)
    cols = len(Y)
    
    amp_array =  [[0 for m in range(cols)] for n in range(rows)]
    
    for i in range(rows):
        for j in range(cols):
            amp_array[i][j] = amplitude_shift_tilt(params, X[i], Y[j], z, modes,a,alpha)
        
    #amp = amplitude_case2(params, X, Y, z, modes,a,alpha)

    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

def calculate_x(params, plane, modes, z):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude_x(params, X, Y, z, modes)
    
    #get phase from amp.
    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

def calculate_x2(params, plane, modes, z):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude_x2(params, X, Y, z, modes)
    
    #get phase from amp.
    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def calculate_case2(params, plane, modes, z,a,alpha):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude_case2(params, X, Y, z, modes,a,alpha)

    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

# Calculate Amplitude and Phase from user x,y range and z. Parse for plots.
def calculate_case2_no_scatter(params, plane, modes, z,a,alpha):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude_case2_no_scatter(params, X, Y, z, modes,a,alpha)

    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)

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
        return '{}{}\n\n{}{}\n\n{}{}\n\n{}{}\n\n{}{}\n\n{}{}'.format('PARAMS: ', self.params,
                                                                     'PLANE: ', self.plane,
                                                                     'MODES: ', self.modes,
                                                                     'Z: ', self.z,
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
def amplitude(params, x, y, z, modes):
    # Unm, U spatial sum over HG modes n,m
    Unm_Sum = 0
    Unm = 0
    rows = len(modes)
    cols = len(modes[0])

    k = params.getK()
    r_c = radius_curvature(z, params)
    w_z = w(z,params)
    gouy = gouy_phase(z,params)

    
    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            #coefficient for mode (n,m) from modes array
            coeff = modes[n][m]

            # neglect 0 coefficient modes
            if(coeff != 0):
                
                # n array for hermpol
                carrN = [0] * rows
                carrN[n] = modes[n][m]
                
                # m array for hermpol
                carrM = [0] * cols
                carrM[m] = modes[n][m]
                
                # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
                if (carrN[n] != 0):
                    carrN[n] = 1
                if (carrM[m] != 0):
                    carrM[m] = 1
                    
                order = n + m

                Unm = (
                    coeff*
                    HG(n,m,params, x, y, z, modes,rows,cols)
                )
                
                # Add each to result
                Unm_Sum += Unm

    return (Unm_Sum)

# Amplitutude calculation from w0,Zr. (LR eq. 9.26)
def amplitude_x(params, x, y, z, modes):
    # Unm, U spatial sum over HG modes n,m
    Unm_Sum = 0
    Unm = 0
    rows = len(modes)
    cols = len(modes[0])

    k = params.getK()
    r_c = radius_curvature(z, params)
    w_z = w(z,params)
    gouy = gouy_phase(z,params)

    
    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            #coefficient for mode (n,m) from modes array
            coeff = modes[n][m]

            # neglect 0 coefficient modes
            if(coeff != 0):
                
                # n array for hermpol
                carrN = [0] * rows
                carrN[n] = modes[n][m]
                
                # m array for hermpol
                carrM = [0] * cols
                carrM[m] = modes[n][m]
                
                # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
                if (carrN[n] != 0):
                    carrN[n] = 1
                if (carrM[m] != 0):
                    carrM[m] = 1
                    
                order = n + m

                Unm = (
                    x*coeff*
                    HG(n,m,params, x, y, z, modes,rows,cols)
                )
                
                # Add each to result
                Unm_Sum += Unm

    return (Unm_Sum)

# Amplitutude calculation from w0,Zr. (LR eq. 9.26)
def amplitude_x2(params, x, y, z, modes):
    # Unm, U spatial sum over HG modes n,m
    Unm_Sum = 0
    Unm = 0
    rows = len(modes)
    cols = len(modes[0])

    k = params.getK()
    r_c = radius_curvature(z, params)
    w_z = w(z,params)
    gouy = gouy_phase(z,params)

    
    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            #coefficient for mode (n,m) from modes array
            coeff = modes[n][m]

            # neglect 0 coefficient modes
            if(coeff != 0):
                
                # n array for hermpol
                carrN = [0] * rows
                carrN[n] = modes[n][m]
                
                # m array for hermpol
                carrM = [0] * cols
                carrM[m] = modes[n][m]
                
                # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
                if (carrN[n] != 0):
                    carrN[n] = 1
                if (carrM[m] != 0):
                    carrM[m] = 1
                    
                order = n + m

                Unm = (
                    (x**2)*coeff*
                    HG(n,m,params, x, y, z, modes,rows,cols)
                )
                
                # Add each to result
                Unm_Sum += Unm

    return (Unm_Sum)

## Get herm polys from modes
# mode: working mode (as Calculate iterates through n,m grid
# coord: x(n) or y(m) coordinate space
# carr: coefficient array
def herm_poly(mode, coord, coeff_array, w_z):
    herm_argument=(coord * sqrt(2) / w_z)
    herm = hermval(herm_argument, coeff_array)    
    return herm

# mode: n or m
# coord: x or y
# calculate hermite polynomials H_n(sqrt(2)x/w(z)) or H_m(sqrt(2)y/w(z)) 

    
#only for n
def herm_poly_case2(mode, coord,a, coeff_array, z, w_z):
    herm_argument=( (coord - a) * sqrt(2) / w_z)
    herm = hermval(herm_argument, coeff_array)
    return herm

# # Amplitutude calculation from w0,Zr. (LR eq. 9.26)
# def amplitude_case2(params, x, y, z, modes,a,alpha):
#     # Unm, U spatial sum over HG modes n,m
#     Unm_Sum = 0
#     rows = len(modes)
#     cols = len(modes[0])

#     k = params.getK()
#     r_c = radius_curvature(z, params)
#     w_z = w(z,params)
#     gouy = gouy_phase(z,params)
    
#     # Iterate through modes
#     for n in range(rows):
#         for m in range(cols):

#             #coefficient for mode (n,m) from modes array
#             coeff = modes[n][m]

#             # neglect 0 coefficient modes
#             if(coeff != 0):
                
#                 # n array for hermpol
#                 carrN = [0] * rows
#                 carrN[n] = modes[n][m]
                
#                 # m array for hermpol
#                 carrM = [0] * cols
#                 carrM[m] = modes[n][m]
                
#                 # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
#                 if (carrN[n] != 0):
#                     carrN[n] = 1
#                 if (carrM[m] != 0):
#                     carrM[m] = 1
                    
#                 order = n + m

#                 Unm = (
#                     coeff*
#                     (
#                         1 / sqrt( 2**(order-1)*factorial(n)*factorial(m)*pi )
#                     ) 
#                     * (1 / w_z ) 
#                     * e**( 
#                             (1j)*(order+1)*gouy 
#                          ) 
#                     * e**
#                     ( 
#                         (
#                             -(1j)*(k * ( (x-a)**2 + y**2) ) 
#                             / 
#                             (2 * r_c)
#                         ) 
#                         - 
#                         (
#                             ( (x-a) **2 + y**2) / w_z**2
#                         )
#                     ) 
#                     * herm_poly_case2(n, x,a, carrN, z, w_z) 
#                     * herm_poly(m, y, carrM, w_z) 
#                     * e**(-(1j)*k*(z-x*np.sin(alpha) )) 
#                 )
                
#                 # Add each to result
#                 Unm_Sum += Unm

#     return (Unm_Sum)

# Amplitutude calculation from w0,Zr. (LR eq. 9.26)
def amplitude_case2(params, x, y, z, modes,a,alpha):
    # Unm, U spatial sum over HG modes n,m
    

    Unm_Sum = 0
    rows = len(modes)
    cols = len(modes[0])
    
    z_r = params.getZr()
    k = params.getK()
    r_c = radius_curvature(z, params)
    #r_c = z+z_r**2/z
    w_z = w(z,params)
    #w_z = params.w0**2*sqrt(1-(z/params.Zr)**2)
    gouy = gouy_phase(z,params)
    
        #shift and tilt transformation on coord's
    x_sub = x
    z_sub = z
#     x = (x+a)*np.cos(alpha) + z*np.sin(alpha)
#     z = z*np.cos(alpha)-(x_sub+a)*alpha
    x = (x_sub+a + z*alpha)
    z = (z-x_sub*alpha)  

    
    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            #coefficient for mode (n,m) from modes array
            coeff = modes[n][m]

            # neglect 0 coefficient modes
            if(coeff != 0):
                
                # n array for hermpol
                carrN = [0] * rows
                carrN[n] = modes[n][m]
                
                # m array for hermpol
                carrM = [0] * cols
                carrM[m] = modes[n][m]
                
                # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
                if (carrN[n] != 0):
                    carrN[n] = 1
                if (carrM[m] != 0):
                    carrM[m] = 1
                    
                order = (n + m)

                Unm = (
                    coeff*
                    (
                        1 / sqrt( 2**(order-1)*factorial(n)*factorial(m)*pi )
                    ) 
                    * (1 / w_z ) 
                    * e**( 
                            (1j)*(order+1)*gouy 
                         ) 
                    * e**
                    ( 
                        (
                            -(1j)*(k * ( x_sub**2 + y**2) ) 
                            / 
                            (2 * r_c)
                        ) 
                        - 
                        (
                            ( x_sub**2 + y**2) / w_z**2
                        )
                    ) 
                    * herm_poly(n, x_sub,carrN, w_z) 
                    * herm_poly(m, y, carrM, w_z) 
                    * e**(-(1j)*k*(z)) 
                )
                
                # Add each to result
                Unm_Sum += Unm

    return (Unm_Sum)

# Amplitutude calculation from w0,Zr. (LR eq. 9.26)
def amplitude_shift_tilt(params, x, y, z, modes,a,alpha):
    # Unm, U spatial sum over HG modes n,m

    #shift and tilt transformation on coord's
#    x_sub = x
#    z_sub = z
#     x = ( (x_sub+a)*np.cos(alpha) + z_sub*np.sin(alpha) )
#     z = ( z_sub*np.cos(alpha)-(x_sub+a)*np.sin(alpha) )
#    x = (x_sub + z*alpha)
#    z = (z-(x_sub+a)*alpha)  

    Unm_Sum = 0
    rows = len(modes)
    cols = len(modes[0])
    
    z_r = params.getZr()
    k = params.getK()
    r_c = radius_curvature(z, params)
    w_z = w(z,params)
    gouy = gouy_phase(z,params)
    

    
    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            #coefficient for mode (n,m) from modes array
            coeff = modes[n][m]

            # neglect 0 coefficient modes
            if(coeff != 0):
                
                # n array for hermpol
                carrN = [0] * rows
                carrN[n] = modes[n][m]
                
                # m array for hermpol
                carrM = [0] * cols
                carrM[m] = modes[n][m]
                
                # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
                if (carrN[n] != 0):
                    carrN[n] = 1
                if (carrM[m] != 0):
                    carrM[m] = 1
                    
                order = (n + m)

                Unm = (
                    coeff*
                    (
                        1 / sqrt( 2**(order-1)*factorial(n)*factorial(m)*pi )
                    ) 
                    * (1 / w_z ) 
                    * e**( 
                            (1j)*(order+1)*gouy 
                         ) 
                    * e**
                    ( 
                        (
                            -(1j)*(k * ( x**2 + y**2) ) 
                            / 
                            (2 * r_c)
                        ) 
                        - 
                        (
                            ( x**2 + y**2) / w_z**2
                        )
                    ) 
                    * herm_poly(n, x, carrN, w_z) 
                    * herm_poly(m, y, carrM, w_z) 
                    * e**(-(1j)*k*z) 
                )
                
                # Add each to result
                Unm_Sum += Unm

    return (Unm_Sum)

#Amplitutude calculation from w0,Zr. (LR eq. 9.26)
#shift is x+a
def amplitude_case2_no_scatter(params, x, y, z, modes,a,alpha):
    # Unm, U spatial sum over HG modes n,m
    Unm_Sum = 0
    Unm_nminus1 = 0
    Unm_second_order = 0
    rows = len(modes)
    cols = len(modes[0])

    k = params.getK()
    r_c = radius_curvature(z, params)
    w_z = w(z,params)
    w_0 = w(0,params)
    gouy = gouy_phase(z,params)
    
    # Iterate through modes
    for n in range(rows):
        for m in range(cols):

            #coefficient for mode (n,m) from modes array
            coeff = modes[n][m]

            # neglect 0 coefficient modes
            if(coeff != 0):
                
                # Arrays for hermpoly
                carrN = [0] * rows
                carrM = [0] * cols
                
                # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes. Factor coeff. after.
                carrN[n] = 1
                carrM[m] = 1
                    
                Unm = (coeff*
                    HG(n,m,params, x, y, z, modes,rows,cols)
                    * (1 
                       - x*2*a/(w_z**2)
                       + x*(1j)*k*alpha 
                       - (x**2)*(1j)*2*k*a*alpha/(w_z**2) 
                      )
                )
                
                # n array for hermpol
                if(n>0):
                    Unm_nminus1 = (coeff*
                        (a*2*sqrt(n)/w_z)
                        * e**((1j)*gouy)
                        *(1+(1j)*k*alpha*x)
                        *HG( (n-1),m,params, x, y, z, modes,rows,cols)
                    )
                
                
                first=0
                second=0
                if(n>0):
                    first = HG( (n-1),m,params, x, y, z, modes,rows,cols)
                    
                if(n>1):
                    second = (HG( (n-2),m,params, x, y, z, modes,rows,cols)*
                    (4/w_z**2)
                                                    * sqrt(n*(n-1))
                                                    *e**(2*(1j)*gouy) 
                             )
                
                Unm_second_order = (
                                    coeff/2*
                                        (
                                             a**2*(1+(1j)*k*alpha*x)
                                            *(
                                                 (4*(x**2)-2*(w_z**2))/(w_z**4)
                                                    *HG(n,m,params, x, y, z, modes,rows,cols)   
                                                 - 8*x/(w_z**3)*(sqrt(n))*e**((1j)*gouy)
                                                    *first
                                                 +  second
                                            )
                                        )                                   
                                    )
                
                #Add each to result
                Unm_Sum += Unm + Unm_nminus1 + Unm_second_order

    return (Unm_Sum)

# def amplitude_case2_no_scatter(params, x, y, z, modes,a,alpha):
#     # Unm, U spatial sum over HG modes n,m
#     Unm_Sum = 0
#     Unm_nminus1 = 0
#     rows = len(modes)
#     cols = len(modes[0])

#     k = params.getK()
#     r_c = radius_curvature(z, params)
#     w_z = w(z,params)
#     w_0 = w(0,params)
#     gouy = gouy_phase(z,params)
    
#     # Iterate through modes
#     for n in range(rows):
#         for m in range(cols):

#             #coefficient for mode (n,m) from modes array
#             coeff = modes[n][m]

#             # neglect 0 coefficient modes
#             if(coeff != 0):
                
#                 # n array for hermpol
#                 carrN = [0] * rows
#                 carrN[n] = modes[n][m]
                
#                 carrN_minus1 = [0] * rows
#                 carrN_minus1[n] = modes[n][m]
                
#                 # m array for hermpol
#                 carrM = [0] * cols
#                 carrM[m] = modes[n][m]
                
#                 # Only calculate Herm pol. sum w/ coeff. 1 for non-zeroes
#                 if (carrN[n] != 0):
#                     carrN[n] = 1
#                     carrN_minus1[(n-1)]=1
#                 if (carrM[m] != 0):
#                     carrM[m] = 1
                    
#                 order = n + m

#                 Unm = (coeff*
#                         (
#                         (
#                             1 / sqrt( 2**(order-1)*factorial(n)*factorial(m)*pi )
#                         ) 
#                         * (1 / w_z ) 
#                         * e**( 
#                                 (1j)*(order+1)*gouy 
#                              ) 
#                         * e**
#                         ( 
#                             (
#                                 -(1j)*(k * (x**2 + y**2) ) 
#                                 / 
#                                 (2 * r_c)
#                             ) 
#                             - 
#                             (
#                                 (x**2 + y**2) / w_z**2
#                             )
#                         ) 
#                         * herm_poly(n, x, carrN, w_z) 
#                         * herm_poly(m, y, carrM, w_z) 
#                         * e**(-(1j)*k*z)
#                         * e**((1j)*gouy)
#                     )
#                     * (1 
#                        - x*2*a/(w_z**2)
#                        + x*(1j)*k*alpha 
#                        - (x**2)*(1j)*2*k*a*alpha/(w_z**2) 
#                       )
#                 )
                
#                 carrN[n]=0
#                 # n array for hermpol
#                 if(n>0):
#                     Unm_nminus1 = (coeff*
#                         (a*2*sqrt(n)/w_z)
#                         * e**((1j)*gouy)
#                         *(1+(1j)*k*alpha*x)
#                         *(
#         (
#             1 / sqrt( 2**(order-2)*factorial(n-1)*factorial(m)*pi )
#         ) 
#         * (1 / w_z ) 
#         * e**( 
#                 (1j)*(order)*gouy 
#              ) 
#         * e**
#         ( 
#             (
#                 -(1j)*(k * (x**2 + y**2) ) 
#                 / 
#                 (2 * r_c)
#             ) 
#             - 
#             (
#                 (x**2 + y**2) / w_z**2
#             )
#         ) 
#         * herm_poly(n-1, x, carrN_minus1, w_z) 
#         * herm_poly(m, y, carrM, w_z) 
#         * e**(-(1j)*k*z)
#         * e**((1j)*gouy)
#     )
#                     )
#                     print(Unm_nminus1)
#                 # Add each to result
#                 Unm_Sum += Unm + Unm_nminus1

#     return (Unm_Sum)


#Calculate HG modes
def HG(n,m,params, x, y, z, modes,rows,cols):
    # Unm, U spatial sum over HG modes n,m

    k = params.getK()
    r_c = radius_curvature(z, params)
    w_z = w(z,params)
    w_0 = w(0,params)
    gouy = gouy_phase(z,params)
   
    Unm = 0
    
    # n array for hermpol
    carrN = [0] * rows

    # m array for hermpol
    carrM = [0] * cols

    # HG is called for non-zero coeff. The coeff. is factored in by calling fxn.
    carrN[n] = 1
    carrM[m] = 1

    order = n + m

    Unm = (
        (
            1 / sqrt( 2**(order-1)*factorial(n)*factorial(m)*pi )
        ) 
        * (1 / w_z ) 
        * e**( 
                (1j)*(order+1)*gouy 
             ) 
        * e**
        ( 
            (
                -(1j)*(k * (x**2 + y**2) ) 
                / 
                (2 * r_c)
            ) 
            - 
            (
                (x**2 + y**2) / w_z**2
            )
        ) 
        * herm_poly(n, x, carrN, w_z) 
        * herm_poly(m, y, carrM, w_z) 
        * e**(-(1j)*k*z)
        * e**((1j)*gouy)
    )
               
    return (Unm)

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


## from q param

# Amplitude from q-parameter (L.R., Eq. 9.34)
def amplitude_q(params, x, y, z, modes):
    # Unm a sum over modes
    UnmSum = 0
    rows = len(modes)
    cols = len(modes[0])

    for n in range(rows):
        for m in range(cols):
            carrN = [0] * rows
            carrN[n] = modes[n][m]

            # print carrN

            carrM = [0] * cols
            carrM[m] = modes[n][m]

            order = n + m

            # Avoid double-counting coefficient
            if (carrN[n] != 0):
                carrN[n] = 1

            Unm = (2 / pi) ** (0.25) * \
                  np.sqrt(1.0 / (2 ** n * factorial(n) * params.getW0())) * \
                  np.sqrt(params.getQ0() / q(z, params)) * \
                  ((params.getQ0() * np.conjugate(q(z, params))) / (np.conjugate(params.getQ0()) * q(z, params))) ** (
                              n / 2.0) * \
                  herm_poly(n, x, carrN, z, params) * \
                  np.exp((-((1j) * params.getK() * x ** 2) / (2 * q(z, params)))) * \
                  (2 / pi) ** (0.25) * \
                  np.sqrt(1.0 / (2 ** m * factorial(m) * params.getW0())) * \
                  np.sqrt(params.getQ0() / q(z, params)) * \
                  ((params.getQ0() * np.conjugate(q(z, params))) / (np.conjugate(params.getQ0()) * q(z, params))) ** (
                              m / 2.0) * \
                  herm_poly(m, y, carrM, z, params) * \
                  np.exp((-((1j) * params.getK() * y ** 2) / (2 * q(z, params))))

            UnmSum += Unm
            # print(Unm, UnmSum)

    return (UnmSum)


def calculate_q(params, plane, modes, z):
    X, Y = np.meshgrid(plane.getX(), plane.getY())
    amp = amplitude_q(params, X, Y, z, modes)

    result_row = len(amp)
    result_col = len(amp[0])
    phase = np.zeros((result_row,result_col), dtype=float)
    for r in range(result_row):
        for c in range(result_col):
            phase[r,c] = (np.angle(amp[r][c]))
    #phase = Phase2(params, X, Y, z, modes)

    f = Result(params, plane, modes, z, amp, phase)

    return (f)



