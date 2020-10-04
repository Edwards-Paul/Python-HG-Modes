# Takes in mat file with array of coeffients and returns a 2D list of HG modes

from hg_scripts import paulisa as pl
import scipy, numpy as np, matplotlib.pyplot as plt
from math import pi, log, exp, sin, cos, atan, e, radians, degrees,factorial as Factorial


def coeff_from_mat(filename,waist,plotOrNot):
    # mat = scipy.io.loadmat('Top_hat_for_paul.mat')
    mat = scipy.io.loadmat(filename)

    coef=mat['coeftopUnitInt'].ravel()
    
    # coef=mat['coeftop'].ravel()

    def N_f(A):

        res = np.floor((np.sqrt(8*A+1)-1)/2)

        #res = (np.sqrt(9+8*A)-3)/2
        return(res)

    def m(N,A):
        res = (N+1)*(N+2)/2 - (A+1)
        return(res)

    #def n(N,A):
    #    m= (N+1)*(N+2)/2 - (A+1)
    #    res = N-m
    #    return(res)
    def n(N,A):
        res = A - (N*(N+1)/2)
        return(res)

    NumberModes = int(len(coef))
    listModesN = [None] * NumberModes
    listModesM = [None] * NumberModes
    listModesC = [None] * NumberModes

    area = pi*waist**2

    #for i in range(len(coef)):
    #    A=i
    #    N= N_f(A)
    #
    #    if (m(N,A)%2 == 0) & (n(N,A)%2 == 0):
    #        print(coef[A], '\t\t\t' , m(N,A), ',' , n(N,A))

    for i in range(NumberModes):
        A=i
        N= N_f(A)
        listModesN[i] = int(m(N,A))
        listModesM[i] = int(n(N,A))
        listModesC[i] = coef[i]
    
    top_modes = pl.create_modes(listModesM,listModesN,listModesC,NumberModes)
    
    if(plotOrNot):
        print(mat['readmepaul'])
        
        for x in mat:
            print(x)
        
        params = pl.Params(1064e-9,.33e-3,0)
        plane = pl.Plane(-2e-3,2e-3,101,-2e-3,2e-3,101)
        # temp_modes=rotate_RX(10e-3,params,10e-3,0,modes)
        # pl.show_modes(modes)
        calc=pl.calculate(params,plane,top_modes,0)
        #plot
        fig, ax = plt.subplots(figsize=(12, 12))
        cs = plt.contourf(calc.plane.getX(), calc.plane.getY(), abs(calc.getAmp() ** 2))

        ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0), useMathText=True)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)

        plt.xlim([-1e-3,1e-3])
        plt.ylim([-1e-3,1e-3])

        cbar = fig.colorbar(cs)

        plt.grid()
        
    return(top_modes)