# This module to evaluate coupling coefficients and return HG mode matrix of (n x m)

from hg_scripts import paulisa as pl
from hg_scripts.transverse_coord_transform import *
from math import pi, log, exp, sin, cos, atan, e, radians, degrees,factorial as Factorial
import cmath
from cmath import sqrt as Sqrt

#represents the sums
def update_modes (z,params,a,b,modes,sols_matrix):
    #get params
    
    K = params.getK()
    w0 = params.getW0()
    zr = params.getZr()
    w = pl.w(z,params)
    gouy = pl.gouy_phase(z,params) 
    R = pl.radius_curvature(z,params)
    
    #build new modes (up to 2 orders larger for quad. dep.)
    rows = len(modes)
    cols = len(modes[0])  
    number_modes = rows*cols
    exp_order = len(sols_matrix)
    new_modes = [[0 for m in range(cols+exp_order)] for n in range(rows+exp_order)]
    
    #placeholders
    x=1
    j=1j
    p=e**(1j*gouy)
    
    #calculate all coupling factors a->e\d
       
    for n in range(len(modes)):
        for m in range(len(modes[n])):
            if (modes[n][m]!=0): #break if c_nm = 0
                c_nm = modes[n][m] #assign c_nm

                for x_order in range(len(sols_matrix[0])):
                    for p_order in range(len(sols_matrix)):

                        if(sols_matrix[x_order][p_order]!='' and (p_order<=n) ):
                            n_start = n-p_order
                            #append each element in x,p matrix to coupling list
                            coupling = eval(sols_matrix[x_order][p_order])*(p**(p_order))
                            #print(p_order,x_order,sols_matrix[x_order][p_order],coupling)

                            #do x transformation 
                            #start at n - [order of p], which x transformation depends on
#                             x_order=0
#                             p_order=0
                            if(x_order>0):
                                q = transform_x(p_order,x_order,a,w,gouy,n,w0,z,zr)
                                #print(n_start,x_order)
                                #empty the q
                                while(q):
                                    item = cp(q.pop())
                                    X = item.coeff #just x coupling
                                    N = item.N #final n order
                                    #print("N,x_order,p_order,n,m",N,x_order,p_order,n,m)
                                    if(N>=0):
                                        new_modes[N][m]+= c_nm*coupling*X
                                        
                            #N is either n or n-p_order, no x-dependence
                            else:
                                N= n-p_order
                                #print("N,x_order,p_order,n,m",N,x_order,p_order,n,m)
                                new_modes[N][m]+= c_nm*coupling                            
# return(new_modes,ind_list)
    return(new_modes)