#!/usr/bin/env python
# coding: utf-8


from hg_scripts import paulisa as pl, pl_plot as plplt,plback as plb, tophat_integration_AW as th
#from Signals_Rc_2 import *
from hg_scripts import build_modes

from numpy import pi as pi
from numpy import angle
from math import e
from numpy import sqrt as sqrt
from scipy.special import erf as erf
import mpmath as mp
import scipy.io

from math import pi, log, exp, sin, cos, atan, e, radians, degrees
from scipy.special import gamma as gamma
from scipy.special import erf as erf
from math import factorial as fact
import cmath
from cmath import sqrt as sqrt
import numpy as np

from time import process_time

inf=np.inf

import pandas as pd
from pprint import pprint

from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    AdaptiveETA, FileTransferSpeed, FormatLabel, Percentage, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer, UnknownLength


# ## Integration for signals



#loop through and sum over modes overlap where a,b -> n,m
def iterate_modes_overlap(w_1,w_2,zR_1,zR_2,k_1,k_2,d_1,W_2,x_1,x_2,Z,Y,modes):   
    
    N = len(modes)
    M = len(modes[0])
    #res_arr = [[0 for i in range(N)] for j in range(M)]
    res = 0

    for a in range(N):
        for b in range(M):            
            #ignore zero coeff.
            if (modes[a][b]!=0):
                c_nm = modes[a][b]
                #result is the sum of all overlap coefficients each with a factor of mode coefficient c_nm
                res += ( c_nm * th.overlap(a,b,w_1,w_2,zR_1,zR_2,k_1,k_2,d_1,W_2,x_1,x_2,Z,Y) )
                
    return(res)



def signals(v,modes,sols_matrix):
    #points determine size of arrays (equivalent to # data points plotted)
    points = v.points
    time_start = process_time()

    time_ave = 0

    #modes_arr = [0]*points
    cl = [0]*points
    cr = [0]*points
    phi_l = [0]*points
    phi_r = [0]*points
    dws = [0]*points
    lps = [0]*points
    total_lps = [0]*points
    
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=points).start() #Start a progress bar for times
    
    for p in range (points):
        time_ave_start = process_time()
        #shift x+zsina
#         new_modes = updated_modes
        new_modes= build_modes.update_modes(v.z,v.params,v.a,v.alpha_arr[p],modes,sols_matrix)

        #create arrays of overlap coefficients left and right
        cl[p] = iterate_modes_overlap(v.w_1,v.w_2,v.zR_1,v.zR_2,v.k_1,v.k_2,v.d_1,v.W_2,v.x_1L,v.x_2L,v.Z,v.Y,new_modes) # left overlap (-2e-3,0)

        cr[p] = iterate_modes_overlap(v.w_1,v.w_2,v.zR_1,v.zR_2,v.k_1,v.k_2,v.d_1,v.W_2,v.x_1R,v.x_2R,v.Z,v.Y,new_modes) # right overlap (0,2e-3)
        time_ave += (process_time()-time_ave_start)

        #create arrays of phases for left and right coeff.
        phi_r[p]=angle(cr[p])
        phi_l[p]=angle(cl[p])
        
        #create arrays of dws &lps for phases in phase arrays
        dws[p] = (phi_r[p]-phi_l[p])
        lps[p] = 0.5*(phi_r[p]+phi_l[p])/v.k_1 *1e9 # Regular LPS in nm.
#         total_lps[p] = (phi_r[p]+phi_l[p])/v.k_1*1e3
        total_lps[p] = angle(cr[p]+cl[p])/v.k_1 *1e9 # Total LPS in nm.
        pbar.update(p+1)
        
    pbar.finish()
    
    #total time for all DWS & LPS points (not including scattering)
    time_elapsed = (process_time() - time_start)
    print(time_elapsed,'s')
    
    
    return(dws,lps,total_lps)


class Vars:
#modes,alpha with points
    def __init__(self, 
                 lam=1064e-9,
                 z_PD=0,z=0,
                 z_m=0,Z=0,
                 z_LO=0, d_1=0,
                 Y=20e-3,
                 x_1R=10e-6,x_2R=20e-3,
                 x_1L=-20e-3,x_2L=-10e-6,
                 w_1=1e-3, w_2=1e-3,#change this? was waist
                 zR_1=pi*1e-3**2/1064e-9,zR_2=pi*1e-3**2/1064e-9, #zr_2 was pi*waist
                 k_1=5905249.348852994,k_2=5905249.348852994,
                 params=pl.Params(1064e-9,1e-3,0),#1e-3 was waist
                 W_2=1,
                 points=101,
                 a=0e-6,alpha_arr=np.linspace(-500e-6,500e-6,101),modes_arr=[0]*101):
        self.lam = lam


        self.z_PD = z_PD #PD location
        self.z = z #prop distance assumed PD

        self.z_m = z_m #meas beam 
        self.Z = Z #distance PD-m

        self.z_LO = z_LO #Local ref. beam
        self.d_1 = d_1 #distance PD-LO


        self.Y = Y #Y int bound

        self.x_1R = x_1R
        self.x_2R = x_2R 

        self.x_1L = x_1L #assumed symmetric with right side
        self.x_2L = x_2L

        self.w_1 = w_1 #waist LO 
        self.w_2 = w_2 #waist MS

        self.zR_1 = zR_1 #rayleigh LO
        self.zR_2 = zR_2


        self.k_1 = k_1 #wavenum LO
        self.k_2 = k_2

        self.params = params #for building tophat coefficients

        self.W_2 = pl.w(z,pl.Params(1064e-9,1e-3,0)) #tophat beam rad , w(z) or 1e-3?; 1e-3 was waist



        #integration and misalignment
        self.points = points

        self.a = a
        self.alpha_arr = alpha_arr

        self.modes_arr = modes_arr

class Gen_Vars:
#modes,alpha with points
    def __init__(self,file_no,lo_size,rx_size,long_off,lat_off,gaps):
        #constants
        self.k_1 = 5905249.348852994 #wavenum LO
        self.k_2 = 5905249.348852994 #wavenum RX
        self.lam = 1064e-9
        self.z_m = 0 #meas beam 
        self.z_LO = 0 #Local ref. beam
        self.points = 61 #number of data points in tilt

   
        #longitudinal offset terms (RX and LO beam waists assumed at z=0)
        self.z_PD = long_off #PD location
        self.z = long_off #propagation distance assumed PD     
        self.Z = long_off #distance PD-m
        self.d_1 = long_off #distance PD-LO
        
        #PD dimensions
        self.Y = lo_size #Y int bound
        self.x_1R = gaps #start at right gap (pos.)
        self.x_2R = lo_size #end at right PD half (pos.)
        self.x_1L = -lo_size #start at left gap (neg.)
        self.x_2L = -gaps #end at left PD half (neg.)

        self.w_1 = lo_size #waist LO 
        self.w_2 = rx_size #waist RX

        self.zR_1 = pi*lo_size**2/1064e-9 #rayleigh LO
        self.zR_2 = pi*rx_size**2/1064e-9 #rayleigh RX


        self.params = pl.Params(1064e-9,rx_size,long_off) #for building tophat coefficients, z=0 or dist?

        self.W_2 = pl.w(long_off,pl.Params(1064e-9,rx_size,long_off)) #tophat beam rad , w(z) or 1e-3?

        self.file_no = file_no

        #integration and misalignment
        self.a = lat_off
        self.alpha_arr = np.linspace(-150e-6,150e-6,61) #rotation angles
        self.modes_arr = [0]*61


