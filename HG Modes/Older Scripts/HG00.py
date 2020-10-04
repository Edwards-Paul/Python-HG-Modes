import PauLisa as pl, numpy as np, matplotlib.pyplot as plt
from scipy.special import erfi as erfi
from scipy.special import erf as erf

pi=np.pi
lam=1
alpha=1
a=1
w=1

  
    #function    
def dws(alpha,a,w,lam,gap):
    size = (int(len(alpha)))
    dws_result=[0]*size
    for i in range(size):
        dws_result[i]=(phi_r(alpha[i],a,w,lam,gap)-phi_l(alpha[i],a,w,lam,gap))
    return(dws_result)

def lps(alpha,a,w,lam,gap):
    size = (int(len(alpha)))
    lps_result=[0]*size
    for i in range(size):
        lps_result[i]=(0.50*(phi_r(alpha[i],a,w,lam,gap)+phi_l(alpha[i],a,w,lam,gap) ))
        
    return(lps_result)

def phi_r(alpha,a,w,lam,gap):
    r= gap
    return(np.arctan 
            (
                (
                pi*alpha/lam
                *
                (
                a
                    *
                    (
                        -erf(np.sqrt(2)*r/w)
                        +
                        1
                        +
                        (2*np.sqrt(2)*r*np.exp(-2*r**2/(w**2)))
                        /
                        (np.sqrt(pi)*w)
                    )
                +
                w*np.sqrt(2/pi)
                    *
                    (
                        np.exp(-2*r**2/w**2)
                    )
                )
            )
                /
             (
                 -erf(np.sqrt(2)*r/w)
                 +
                 1
                 +
                 np.sqrt(2/pi)*(a/w)
                 *
                 np.exp(-2*r**2/w**2)
             )      
            )
    )

def phi_l(alpha,a,w,lam,gap):
    l= (-gap)
    return(np.arctan 
            (
                (
                pi*alpha/lam
                *
                (
                a
                    *
                    (
                        erf(np.sqrt(2)*l/w)
                        +
                        1
                        -
                        (2*np.sqrt(2)*l*np.exp(-2*l**2/(w**2)))
                        /
                        (np.sqrt(pi)*w)
                    )
                -
                w*np.sqrt(2/pi)
                    *
                    (
                        np.exp(-2*l**2/w**2)
                    )
                )
            )
                /
             (
                 erf(np.sqrt(2)*l/w)
                 +
                 1
                 -
                 np.sqrt(2/pi)*(a/w)
                 *
                 np.exp(-2*l**2/w**2)
             )      
            )
    )
