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
from PauLisa import *

# SECOND ORDER TILT, first in shift
# Apply general mode scatter formulas for x-dependence (uni-directional offset coord.), RETURN modes
def rotate_RX_2(z,params,a,alpha,modes):
    
    b = alpha
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
    new_modes = [[0 for m in range(cols+3)] for n in range(rows+3)]
    
    J_0 = (
            1
        -1j*k*a*b
        +1j*b**2*k*z/2
        +(1/w_z**2) *(2*a*b*z-b**2*z**2)
    )
    
    J_1 = (
            1j*k*b
        +(1/w_z**2)*(2*a
                     -2*b*z
                   +4*1j*a*(b**2)*k*z
                   +1j*a*(b**2)*k*z)
    )
    
    J_2 = (
        -(b**2)*(k**2)/2
        +(1/w_z**2)*(
            2*1j*a*b*k
            -2*1j*(b**2)*k*z
            +b**2)
        +(1/w_z**4)*(
            -4*a*b*z
            +2*(b**2)*(z**2)
        )
    )
    
    J_3 = (
        -a*(b**2)*(k**2)/(w_z**2) 
        - 4*1j*a*(b**2)*k*z/(w_z**4)
    )
        
    K_0 = (
           -a
        +a*(b**2)/2
        -1/2*1j*a*(b**2)*k*z
        +b*z
        -1j*a*(b**2)*k*z
        +3*a*(b**2)*(z**2)/(w_z**2)
    )
    
    K_1 = (
        -1j*a*b*k
        +4*a*b*z/(w_z**2)
        -(b**2)/2
        +1j*(b**2)*k*z
        -2*(b**2)*(z**2)/(w_z**2)
    )    

    K_2 = (
        +(1/2)*a*(b**2)*(k**2)
        -2*a*(b**2)/(w_z**2)
        +4*1j*a*(b**2)*k*z/(w_z**2)
        -6*a*(b**2)*(z**2)/(w_z**4)
    )   
    
    L_0 = (
            -2*a*b*z
        +(b**2)*(z**2)
    )
    
    L_1 = (
        a*(b**2)
        - 2*1j*a*(b**2)*k*z
        +6*a*(b**2)*(z**2)/(w_z**2)
    )
    
    # Add u(n+3,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]   
                
                x_3 = x_plus_3(w_z,phi,n)
                #Add at n+3
                new_modes[n+3][m] += (
                    c_nm*(
                        x_3*J_3
                    )
                )
    
    # Add u(n+2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                
                x_2 = x_plus_2(w0,z,zr,n)
                
                #Add at n+2
                new_modes[n+2][m] += (
                    c_nm*(
                        x_2*J_2
                    )
                )
                
    # Add u(n+1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm= modes[n][m]
                
                x_31 = x_plus_31(w_z,phi,n)
                x_1 = x_plus_1(w0,z,zr,n)
                x_2 = x_plus_2(w0,z,zr,n-1)
                
                K = 2*sqrt(n)*phi/w_z
                
                #Add at n+1
                new_modes[n+1][m] += (
                    c_nm*(
                        x_31*J_3
                        + x_1*J_1
                        + x_2*K*K_2
                    )
                )
                
    # Add u(n,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]

                x_1 = x_plus_1(w0,z,zr,n-1)
                x_2_0 = x_zero_2(w0,z,zr,n)
                
                K = 2*sqrt(n)*phi/w_z
        
                 ##### ADDED NEW TERMS FROM OFFSET COEFF
                #Add at n
                new_modes[n][m] += (
                    c_nm*(
                        J_0
                        + x_2_0*J_2
                        + x_1*K*K_1
                    )
                )
                
    # Add u(n-1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>0) ): #neglect n-1
                c_nm = modes[n][m]
                
                x_31 = x_minus_31(w_z,phi,n)
                x_1_m = x_minus_1(w0,z,zr,n) #(n-1)
                x_1_p = x_plus_1(w0,z,zr,n-2) 
                x_2_0 = x_zero_2(w0,z,zr,n-1)
                
                K = 2*sqrt(n)*phi/w_z
                L = 2*sqrt(n*(n-1))*phi**2/w_z**2 
                
                
                #Add at n-1
                new_modes[n-1][m] += (
                    c_nm*(
                        x_31*J_3
                        + x_1_m*J_1
                        + x_2_0*K*K_2
                        + K*K_0
                        + x_1_p*L*L_1
                    ) 
                )
           
                
    # Add u(n-2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>1) ): #neglect n=0,1
                c_nm = modes[n][m]
                
                x_2 = x_minus_2(w0,z,zr,n)
                x_1 = x_minus_1(w0,z,zr,n-1)
                
                K = 2*sqrt(n)*phi/w_z
                L = 2*sqrt(n*(n-1))*phi**2/w_z**2 
                
                #Add at n-2
                new_modes[n-2][m] += (
                    c_nm*(
                        x_2*J_2
                        + x_1*K*K_1
                        + L*L_0
                    ) 
                )
                
    # Add u(n-3,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>2) ): #neglect n=0,1
                c_nm = modes[n][m]
                
                x_3 = x_minus_3(w_z,phi,n)
                x_2 = x_minus_2(w0,z,zr,n-1)
                x_1 = x_minus_1(w0,z,zr,n-2)
                
                K = 2*sqrt(n)*phi/w_z
                L = 2*sqrt(n*(n-1))*phi**2/w_z**2 
                
                #Add at n-3
                new_modes[n-3][m] += (
                    c_nm*(
                        x_3*J_3
                        +x_2*K*K_2
                        +x_1*L*L_1
                    ) 
                )                
    
    #account for symmetry, 0,m_final = m_final,0
    if (rows%2==0 and cols%2==0 and rows==cols):
        new_modes[0][rows+2]=new_modes[rows+2][0]
    
    return(new_modes)


# Apply general mode scatter formulas for x-dependence (uni-directional offset coord.), RETURN modes
def rotate_RX(z,params,a,alpha,modes):
    
    
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
    
    B_n_1 = (
        -a*2/w_z**2
        -a*1j*k/r_c
        +1j*k*alpha
        -alpha*2*z/w_z**2
        -alpha*1j*k*z/r_c
    )
    
    B_n_2 = (
        -alpha*a*1j*2*k/w_z**2
        +a*alpha*k**2/r_c
        +a*alpha*4*z/w_z**4
        +(a*alpha*4*1j*k*z)/(w_z**2 *r_c)
        -a*alpha*k**2 *z/(r_c**2)
    )
    
    # Add u(n+2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]              
                #Add at n+2
                new_modes[n+2][m] += (
                    c_nm*(x_plus_2(w0,z,zr,n)*B_n_2)
                )
                
    # Add u(n+1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm= modes[n][m]
                #Add at n+1
                new_modes[n+1][m] += (
                    c_nm*(x_plus_1(w0,z,zr,n) * B_n_1)
                )
                
    # Add u(n,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                
                #removed 1+ on each term, these are at B_nminus1_0
                B_nminus1_1 = (
                    (a*2*sqrt(n)/w_z)*phi
                    *
                    ( 1j*k*alpha
                     - 2*z*alpha/w_z**2
                     - 1j*k*z*alpha/r_c
                    )
                    +
                    (alpha*z*2*sqrt(n)/w_z)*phi
                    *
                    (
                     - 2*a/w_z**2
                     - 1j*k*a/r_c
                    )
                )
                
                new_terms = (
                      (1j)*k*a*alpha
                    - 2*z*a*alpha/w_z**2
                    - (1j)*k*z*a*alpha/r_c
                )
                 ##### ADDED NEW TERMS FROM OFFSET COEFF
                #Add at n
                new_modes[n][m] += (
                    c_nm*(1 + new_terms + x_plus_1(w0,z,zr,n-1)*B_nminus1_1 + x_zero_2(w0,z,zr,n)*B_n_2 )
                )
                
    # Add u(n-1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>0) ): #neglect n-1
                c_nm = modes[n][m]
                
                B_nminus1_0 = (
                    (a*2*sqrt(n)/w_z)*phi
                    +
                    (alpha*z*2*sqrt(n)/w_z)*phi
                )
        
                #Add at n
                new_modes[n-1][m] += (
                    c_nm*(B_nminus1_0 + x_minus_1(w0,z,zr,n)*B_n_1 ) 
                )
           
                
    # Add u(n-2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>1) ): #neglect n=0,1
                c_nm = modes[n][m]
                
                B_nminus1_1 = (
                (a*2*sqrt(n)/w_z)*phi
                *
                (
                    1j*k*alpha
                 - 2*z*alpha/w_z**2
                 - 1j*k*z*alpha/r_c
                )
                                        +
                (alpha*z*2*sqrt(n)/w_z)*phi
                *
                (
                 - 2*a/w_z**2
                 - 1j*k*a/r_c
                )
                )
                
                B_nminus2_0 = (
                    (z*alpha*2*sqrt(n-1)*phi/w_z) * (a*2*sqrt(n)*phi/w_z) 
                )
       
                #Add at n
                new_modes[n-2][m] += (
                    c_nm*(B_nminus2_0 + x_minus_1(w0,z,zr,n-1)*B_nminus1_1 + x_minus_2(w0,z,zr,n)*B_n_2 ) 
                )
    
    #account for symmetry, 0,m_final = m_final,0
    if (rows%2==0 and cols%2==0 and rows==cols):
        new_modes[0][rows+2]=new_modes[rows+2][0]
    
    return(new_modes)

# Apply general mode scatter formulas for x-dependence (uni-directional offset coord.), RETURN modes
def rotate_LO(z,params,a,alpha,modes):
    
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
    
    B_n_1 = (
        -a*2/w_z**2
        -a*1j*k/r_c
        +1j*k*alpha
        -alpha*2*z/w_z**2
        -alpha*1j*k*z/r_c
    )
    
    B_n_2 = (
        -alpha*a*1j*2*k/w_z**2
        +a*alpha*k**2/r_c
        +a*alpha*4*z/w_z**4
        +(a*alpha*4*1j*k*z)/(w_z**2 *r_c)
        -a*alpha*k**2 *z/(r_c**2)
    )
    
    
    
    # Add u(n+2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]              
                #Add at n+2
                new_modes[n+2][m] += (
                    c_nm*(x_plus_2(w0,z,zr,n)*B_n_2)
                )
                
    # Add u(n+1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm= modes[n][m]
                #Add at n+1
                new_modes[n+1][m] += (
                    c_nm*(x_plus_1(w0,z,zr,n) * B_n_1)
                )
                
    # Add u(n,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0):
                c_nm = modes[n][m]
                
                #removed 1+ on each term, these are at B_nminus1_0
                B_nminus1_1 = (
                    (a*2*sqrt(n)/w_z)*phi
                    *
                    (
                     + 1j*k*alpha
                     - 2*z*alpha/w_z**2
                     - 1j*k*z*alpha/r_c
                    )
                    +
                    (alpha*z*2*sqrt(n)/w_z)*phi
                    *
                    (
                     - 2*a/w_z**2
                     - 1j*k*a/r_c
                    )
                )
                
                new_terms = (
                    0
                )
                 ##### ADDED NEW TERMS FROM OFFSET COEFF
                #Add at n
                new_modes[n][m] += (
                    c_nm*(1 + new_terms + x_plus_1(w0,z,zr,n-1)*B_nminus1_1 + x_zero_2(w0,z,zr,n)*B_n_2 )
                )
                
    # Add u(n-1,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>0) ): #neglect n-1
                c_nm = modes[n][m]
                
                B_nminus1_0 = (
                    (a*2*sqrt(n)/w_z)*phi
                    +
                    (alpha*z*2*sqrt(n)/w_z)*phi
                )
        
                #Add at n
                new_modes[n-1][m] += (
                    c_nm*(B_nminus1_0 + x_minus_1(w0,z,zr,n)*B_n_1 ) 
                )
           
                
    # Add u(n-2,m)
    for n in range(rows):
        for m in range(cols):            
            #ignore zero coeff.
            if (modes[n][m]!=0 and (n>1) ): #neglect n=0,1
                c_nm = modes[n][m]
                
                B_nminus1_1 = (
                (a*2*sqrt(n)/w_z)*phi
                *
                (
                    1j*k*alpha
                 - 2*z*alpha/w_z**2
                 - 1j*k*z*alpha/r_c
                )
                                        +
                (alpha*z*2*sqrt(n)/w_z)*phi
                *
                (
                 - 2*a/w_z**2
                 - 1j*k*a/r_c
                )
                )
                
                B_nminus2_0 = (
                    (z*alpha*2*sqrt(n-1)*phi/w_z) * (a*2*sqrt(n)*phi/w_z) 
                )
       
                #Add at n
                new_modes[n-2][m] += (
                    c_nm*(B_nminus2_0 + x_minus_1(w0,z,zr,n-1)*B_nminus1_1 + x_minus_2(w0,z,zr,n)*B_n_2 ) 
                )
    
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

#x**3 dep
def x_plus_3(w_z,phi,n):
    factor = (w_z**3) * sqrt((n+3)*(n+2)*(n+1)) *phi**(-3)/8
    return(factor)

def x_plus_31(w_z,phi,n):
    factor = 3*(n+1)**(3/2)*(w_z**3)* phi**(-1)/8 
    return(factor)

def x_minus_31(w_z,phi,n):
    factor = 3*n**(3/2)*(w_z**3)*phi/8
    return(factor)

def x_minus_3(w_z,phi,n):
    factor = sqrt(n*(n-1)*(n-2))*(w_z**3)*(phi**3)/8
    return(factor)

