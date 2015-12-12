import numpy as np
import pylab as pl
import pyfits as py
from stellarpop import tools, distances
import itertools

dist = distances.Distance()


def MassK(K,z,age='4.000'):
    sed = tools.getSED('BC_Z=1.0_age='+age+'gyr')
    kfilt = tools.filterfromfile('Kp_NIRC2')
    K_modrest = tools.ABFM(kfilt,sed,0.0) 
    K_modobs = tools.ABFM(kfilt,sed,z)
    K_rest = K + K_modrest - K_modobs
    mass_K = 0.4*(K_rest - K_modrest)
    K_sun = 5.19
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    K_restabs = K_rest - DM
    lum = -0.4*(K_restabs - K_sun)
    return lum, K_rest

def MassR(R,z,age='4.000'):
    sed = tools.getSED('BC_Z=1.0_age='+age+'gyr')
    rfilt = tools.filterfromfile('F814W_ACS')  
    R_modrest = tools.ABFM(rfilt,sed,0.0)
    R_modobs = tools.ABFM(rfilt,sed,z)
    R_rest = R + R_modrest - R_modobs
    mass_R = 0.4*(R_rest - R_modrest)
    R_sun = 4.57
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    R_restabs = R_rest - DM
    lum = -0.4*(R_restabs - R_sun)
    return lum, R_rest

def MassB1(B,z,age='4.000'):
    sed = tools.getSED('BC_Z=1.0_age='+age+'gyr')
    bfilt = tools.filterfromfile('F606W_ACS')
    B_modrest = tools.ABFM(bfilt,sed,0.0) 
    B_modobs = tools.ABFM(bfilt,sed,z)
    B_rest = B + B_modrest - B_modobs
    mass_B = 0.4*(B_rest - B_modrest)
    B_sun = 4.74
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    B_restabs = B_rest - DM
    lum = -0.4*(B_restabs - B_sun)
    return lum, B_rest

def MassB2(B,z,age='4.000'):
    sed = tools.getSED('BC_Z=1.0_age='+age+'gyr')
    bfilt = tools.filterfromfile('F555W_ACS')
    B_modrest = tools.ABFM(bfilt,sed,0.0) 
    B_modobs = tools.ABFM(bfilt,sed,z)
    B_rest = B + B_modrest - B_modobs
    mass_B = 0.4*(B_rest - B_modrest)
    B_sun = 4.83
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    B_restabs = B_rest - DM
    lum = -0.4*(B_restabs - B_sun)
    return lum, B_rest

def MassB1toB(B,z,age='4.000'):
    sed = tools.getSED('BC_Z=1.0_age='+age+'gyr')
    bfilt = tools.filterfromfile('F606W_ACS')
    restfilt = tools.filterfromfile('F435W_ACS')
    B_modrest = tools.ABFM(restfilt,sed,0.0) 
    B_modobs = tools.ABFM(bfilt,sed,z)
    B_rest = B + B_modrest - B_modobs
    mass_B = 0.4*(B_rest - B_modrest)
    B_sun = 5.372 # ref: EzGal
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    B_restabs = B_rest - DM
    lum = -0.4*(B_restabs - B_sun)
    return lum, B_rest

def MassB2toB(B,z,age='4.000'):
    sed = tools.getSED('BC_Z=1.0_age='+age+'gyr')
    bfilt = tools.filterfromfile('F555W_ACS')
    restfilt = tools.filterfromfile('F435W_ACS')
    B_modrest = tools.ABFM(restfilt,sed,0.0) 
    B_modobs = tools.ABFM(bfilt,sed,z)
    B_rest = B + B_modrest - B_modobs
    mass_B = 0.4*(B_rest - B_modrest)
    B_sun = 5.372
    DL = dist.luminosity_distance(z)
    DM = 5.*np.log10(DL*1e6) - 5.
    B_restabs = B_rest - DM
    lum = -0.4*(B_restabs - B_sun)
    return lum, B_rest
