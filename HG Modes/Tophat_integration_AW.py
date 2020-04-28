# ========================================================================================================

# IMPORTS:

## Math imports
from math import pi, log, exp, sin, cos, atan, e, radians, degrees
from scipy.special import gamma as gamma
from scipy.special import erf as erf
from math import factorial as fact
import cmath
from cmath import sqrt as sqrt
import numpy as np


# --------------------------------------------------------------------------------------------------------

#sigma fxn
def sigma(w_1,w_2,zR_1,zR_2,d_1,Z):
    res = (
        1/(w_1**2 *(1+(1j)*d_1/zR_1))
        +
        1/((w_2**2 *(1-(1j)*Z/zR_2)))
    )
    return(res)

#even m G(x;m) fxn
def G_even(x,m):
    
    stop = int((m/2-1)+1) #stop sum,+1 for index array
    start = 0  #start M=0
    sum_M = sum([
           ( fact(M+1)/fact(2*(M+1))
        *(2*x)**(2*M+1) )
           for M in range(0, stop,1)])
     
    res = ( gamma((m+1)/2) * (erf(x)/2 - e**(-x**2)/pi * sum_M) )
    return(res)

#odd m G(x;m) fxn
def G_odd(x,m):
    
    stop = int((m-1)/2 + 1) #stop sum, +1 for index array
    start = 0  #start M=0
    sum_M = sum([
            ( x**(2*M)/fact(M) )
           for M in range(0, stop,1)])
     
    res = ( - fact((m-1)/2) * e**(-x**2)/2 * sum_M)   
    return(res)


# --------------------------------------------------------------------------------------------------------

#a = n, b = m (modes); w = waist, zR=rayleigh,k=wave#,d=z-z0,Z=z,Y=pd dimensions
def overlap(a,b,w_1,w_2,zR_1,zR_2,k_1,k_2,d_1,W_2,x_1,x_2,Z,Y):
    
    s = sigma(w_1,w_2,zR_1,zR_2,d_1,Z) #calculate sigma once
    
    
    prefactor_num = ( 2**(a+b+1) * sqrt(fact(a)*fact(b))*e**(-(1j)*(k_2*Z-k_1*d_1)) )
    
    prefactor_denom = pi*sqrt(s**(2+a+b))*w_1*w_2**(1+a+b)*(1+(1j)*d_1/zR_1)*(1-(1j)*Z/zR_2)**(1+a+b)
    
    prefactor = prefactor_num/prefactor_denom
    
    stop_A = int(np.floor(a/2))+1
    stop_B = int(np.floor(b/2))+1
    sum_terms = (
        sum( (
            sum_first_term(s,A,B,a,b,x_1,x_2,W_2)
            * sum_second_term(s,Y,b,B,x_2,x_1,A,a)
        )
            for A in range(0,stop_A,1) 
            for B in range(0,stop_B,1) )
    )
    
    res = prefactor * sum_terms
    return(res)
    
def sum_first_term(s,A,B,a,b,x_1,x_2,W_2):

    m = a-2*A
    
    sum_ab = 0 
    if m % 2 == 0:
        sum_ab = (
            (
                (-s/8)**(A+B) 
                 * W_2**(2*(A+B)) 
                 * gamma( (b+1)/2-B) 
            )
            /
            (
                fact(A)*fact(B)*fact(a-2*A)*fact(b-2*B)
            )
        )
    
    #even G(x;m) fxn
    if m % 2 == 0:
        g=(G_even(sqrt(s)*x_2,a-2*A)
            -
           G_even(sqrt(s)*x_1,a-2*A)
        )
        
    else:
        g=(G_odd(sqrt(s)*x_2,a-2*A)
            -
           G_odd(sqrt(s)*x_1,a-2*A)
        )
    
    
    return(sum_ab*g)
    
def sum_second_term(s,Y,b,B,x_2,x_1,A,a):
    res = (
        erf(cmath.sqrt(s)*Y)
        -
        2*e**(-s*Y**2)/sqrt(pi)
        *
        sum_b(b,B,s,Y,x_2,x_1,A,a)
    )
    
    return(res)

def sum_b(b,B,s,Y,x_2,x_1,A,a):
    start = 0 #M=0
    stop = int((b/2-(B+1))) +1 #end(w/ array)
    
    m = a-2*A
    

    res = (
           sum([
            fact(M+1)/fact(2*(M+1))
               * (2*sqrt(s)*Y)**(2*M+1)
           for M in range(0, stop,1)])     
    )
   
    
    return(res)