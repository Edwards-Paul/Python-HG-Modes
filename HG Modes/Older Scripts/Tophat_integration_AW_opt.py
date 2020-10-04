# ========================================================================================================

# IMPORTS:

## Math imports
import math as math
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
    M = int(m/2-1)
    sum_M = 0
    
    if (M>=0):
        sum_M = ( fact(M+1)/fact(2*(M+1))
            *(2*x)**(2*M+1) )

    res = ( gamma((m+1)/2) * (erf(x)/2 - (e**(-x**2))/sqrt(pi) * sum_M) )

    #print('m',m,'x',x,'sum_M',sum_M,'res',res)
    
    return(res)

#odd m G(x;m) fxn
def G_odd(x,m):

    M  = int((m-1)/2)
    sum_M = 0
    
    if (M>=0):
        sum_M = ( (x**(2*M))/fact(M) )
     
    res = ( - fact((m-1)/2) * (e**(-x**2))/2 * sum_M)   
    return(res)


# --------------------------------------------------------------------------------------------------------

#a = n, b = m (modes); w = waist, zR=rayleigh,k=wave#,d=z-z0,Z=z,Y=pd dimensions
def overlap(a,b,w_1,w_2,zR_1,zR_2,k_1,k_2,d_1,W_2,x_1,x_2,Z,Y):
    
    if(b%2==0):
        s = sigma(w_1,w_2,zR_1,zR_2,d_1,Z) #calculate sigma once
        #print('s',s)

        prefactor_num = ( 2**(a+b+1) * sqrt(fact(a)*fact(b))*e**(-(1j)*(k_2*Z-k_1*d_1)) )

        prefactor_denom = pi*sqrt(s**(2+a+b))*w_1*w_2**(1+a+b)*(1+(1j)*d_1/zR_1)*(1-(1j)*Z/zR_2)**(1+a+b)

        prefactor = prefactor_num/prefactor_denom

        A = int(np.floor(a/2))
        B = int(np.floor(b/2))
        print('a','b',a,b,'A','B',A,B)
        sum_terms = (
                sum_first_term(s,A,B,a,b,x_1,x_2,W_2)
                * sum_second_term(s,Y,b,B,x_2,x_1)
             )
        
        
        res = prefactor * sum_terms

        return(res)
    else:
        return(0)
    
#second line, sum A->fl(a/2), B->fl(b/2). Call 3rd line(b,B) and 4th line (a,A)    
def sum_first_term(s,A,B,a,b,x_1,x_2,W_2):

    m = int(a-2*A)
    
    sum_ab = 0 
    g=0
    
    print('a','b',a,b,'A','B',A,B)
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
    if (m % 2 == 0) :
       # print(s,sqrt(s)*x_2,sqrt(s))
        g=(G_even(sqrt(s)*x_2,m)
            -
           G_even(sqrt(s)*x_1,m)
        )
        
    if (m % 2 != 0):
        g=(G_odd(sqrt(s)*x_2,m)
            -
           G_odd(sqrt(s)*x_1,m)
        )
    
    print('g',g,'a',a,'A',A,'b',b,'B',B)
    
    #print(sum_ab,g,(-s/8)**(A+B))
    
    return(sum_ab*g)

#line 3
def sum_second_term(s,Y,b,B,x_2,x_1):
    
    res = (
        erf(cmath.sqrt(s)*Y)
        -
        (2*e**(-s*Y**2)/sqrt(pi)
        *
        sum_b(b,B,s,Y,x_2,x_1)) #call summation on b/2-(B+1)
    )
    #print(res)
    return(res)

#sum M=0 -> b/2-(B+1)
def sum_b(b,B,s,Y,x_2,x_1):
    start = 0 #M=0
    stop = int((b/2-(B+1))) +1 #end(w/ array)
    
    res = (
           sum([
            fact(M+1)/fact(2*(M+1))
               * (2*sqrt(s)*Y)**(2*M+1)
           for M in range(0, stop,1)])     
    )
   
    #print('b',b,'B',B,res)
    return(res)    
#     M = int((b/2-(B+1)))  #end(w/ array)
#     res = 0
    
#     if (M>=0):
#         res = (
#                 fact(M+1)/fact(2*(M+1))
#                    * (2*sqrt(s)*Y)**(2*M+1)    
#         )
   
#     #print('b',b,'B',B,res)
#     return(res)